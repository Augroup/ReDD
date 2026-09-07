import torch
from torch import nn
import torch.nn.functional as F
def init_module_weights(model):
    for m in model.modules():
        if isinstance(m, nn.Conv1d):
            nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
            if m.bias is not None:
                nn.init.zeros_(m.bias)
        elif isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)
            if m.bias is not None:
                nn.init.zeros_(m.bias)
        elif isinstance(m, nn.LSTM):
            for name, param in m.named_parameters():
                if 'weight_ih' in name:
                    nn.init.xavier_uniform_(param.data)
                elif 'weight_hh' in name:
                    nn.init.orthogonal_(param.data)
                elif 'bias' in name:
                    param.data.fill_(0)
                    n = param.size(0)
                    param.data[n//4:n//2].fill_(1.0) # Forget gate bias
            
# --------- Residual CNN Block ---------
class ResidualCNNBlock(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size=5, dropout=0.1):
        super().__init__()
        self.conv = nn.Conv1d(in_channels, out_channels, kernel_size, padding=kernel_size // 2)
        self.bn = nn.BatchNorm1d(out_channels)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(dropout)
        self.shortcut = nn.Conv1d(in_channels, out_channels, 1) if in_channels != out_channels else nn.Identity()
    def forward(self, x):
        residual = self.shortcut(x)
        out = self.conv(x)
        out = self.bn(out)
        out = self.relu(out)
        out = self.dropout(out)
        return out + residual

class DeepCNN(nn.Module):
    def __init__(self, in_channels, layers=[32, 64, 128], kernel_size=5, dropout=0.1):
        super().__init__()
        cnn_layers = []
        cur_in = in_channels
        for out_ch in layers:
            cnn_layers.append(ResidualCNNBlock(cur_in, out_ch, kernel_size, dropout))
            cur_in = out_ch
        self.cnn = nn.Sequential(*cnn_layers)
    def forward(self, x):
        x = x.transpose(1, 2)
        x = self.cnn(x)
        x = x.transpose(1, 2)
        return x

# --------- BiLSTM Block ---------
class BidirectionalLSTM(nn.Module):
    def __init__(self, input_dim, hidden_dim=256, num_layers=2, dropout=0.2):
        super().__init__()
        self.bilstm = nn.LSTM(input_dim, hidden_dim, num_layers=num_layers,
                              bidirectional=True, dropout=dropout, batch_first=True)
    def forward(self, x):
        x, _ = self.bilstm(x)
        return x  # [batch, seq, 2*hidden_dim]

# --------- Positional Encoding ---------
class PositionalEncoding(nn.Module):
    def __init__(self, dim, max_len=500):
        super().__init__()
        pe = torch.zeros(max_len, dim)
        position = torch.arange(0, max_len, dtype=torch.float32).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, dim, 2).float() * -(torch.log(torch.tensor(10000.0)) / dim))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer("pe", pe)
    def forward(self, x):
        return x + self.pe[:, :x.size(1)]

# --------- Transformer Encoder ---------
class TransformerEncoderBlock(nn.Module):
    def __init__(self, input_dim, embed_dim, num_heads=8, num_layers=4, dropout=0.1, max_len=1025):
        super().__init__()
        self.input_proj = nn.Linear(input_dim, embed_dim)
        self.pos_encoding = PositionalEncoding(embed_dim, max_len)
        encoder_layer = nn.TransformerEncoderLayer(embed_dim, num_heads,
                                                  dim_feedforward=embed_dim * 4,
                                                  dropout=dropout, batch_first=True)
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers)
        self.norm = nn.LayerNorm(embed_dim)
    def forward(self, x):
        x = self.input_proj(x)
        x = self.pos_encoding(x)
        x = self.encoder(x)
        x = self.norm(x)
        return x
class AttentionDecoder(nn.Module):
    def __init__(self, embed_dim, num_heads, num_classes, dropout=0.1, is_mod_head=False):
        super().__init__()
        self.mha = nn.MultiheadAttention(embed_dim, num_heads, batch_first=True)
        self.batch_norm = nn.BatchNorm1d(embed_dim)
        self.dropout = nn.Dropout(dropout)
        
        # CHANGE: For the mod head, we output 2 values (alpha and beta)
        # PROPOSED (Refined Reasoning)
        
        self.is_mod_head = is_mod_head
        
        if self.is_mod_head:
#             self.final_fc = nn.Sequential(
#                 nn.Linear(embed_dim, 512),
#                 nn.ReLU(),
#                 nn.Dropout(dropout), # Crucial for 1024-dim models
#                 nn.Linear(512, 256),
#                 nn.ReLU(),
#                 nn.Linear(256, num_classes), # Final Alpha/Beta output
#             )
            self.final_fc = nn.Linear(embed_dim, num_classes)
            self.precise_fc = nn.Linear(embed_dim, num_classes*10)
        else:
            self.final_fc = nn.Linear(embed_dim, num_classes)
        
    def forward(self, x):
        x, _ = self.mha(x, x, x)
        x = self.batch_norm(x.transpose(1, 2)).transpose(1, 2)
        if self.is_mod_head:
            x = x.mean(axis=1) 
            x = self.dropout(x)
            return torch.cat([self.final_fc(x),self.precise_fc(x)],dim=1)
        x = self.dropout(x)
        return self.final_fc(x)
# --------- Main Model ---------
class ReDDModel(nn.Module):
    def __init__(self, featuredim=5, window_size=17, vocab_size=5, mod_classes=1,
                CNN_dims=[8, 16, 32],CNN_kernel_size=5,LSTM_dim=64,
                 transformer_embed=128, transformer_heads=8,transformer_layers=3, 
                 attention_heads=8,
                 dropout_CNN=0.2,dropout_LSTM=0.2,dropout_transformer=0.1,dropout_attention=0.1):
        super().__init__()
        # Feature extractor
        self.cnn = DeepCNN(in_channels=featuredim,
                           layers=CNN_dims, kernel_size=5, dropout=dropout_CNN)
        self.bilstm = BidirectionalLSTM(input_dim=CNN_dims[-1], hidden_dim=LSTM_dim, num_layers=2, dropout=dropout_LSTM)
        self.transformer = TransformerEncoderBlock(input_dim=LSTM_dim*2, embed_dim=transformer_embed,
                                                   num_heads=transformer_heads, num_layers=transformer_layers,
                                                   dropout=dropout_transformer, max_len=window_size)

        # Decoders
        self.attn_ref = AttentionDecoder(embed_dim=transformer_embed, num_heads=attention_heads,
                                         num_classes=vocab_size, dropout=dropout_attention)
        self.attn_call = AttentionDecoder(embed_dim=transformer_embed, num_heads=attention_heads,
                                          num_classes=vocab_size, dropout=dropout_attention)
        self.attn_mod = AttentionDecoder(transformer_embed, attention_heads, 1, 
                                         dropout_attention, is_mod_head=True)
        
    def forward(self, x):
        # x: [batch, seq_len, featuredim]
        x = self.cnn(x)                        # [batch, seq_len, 128]
        x = self.bilstm(x)                     # [batch, seq_len, 512]
        x = self.transformer(x)                # [batch, seq_len, transformer_embed]
        preds_ref = self.attn_ref(x)           # [batch, seq_len, vocab_size]
        preds_call = self.attn_call(x)         # [batch, seq_len, vocab_size]
        preds_mod = self.attn_mod(x)           # [batch, seq_len, mod_classes]
        return preds_ref, preds_call, preds_mod

    def init_weights(self):
        init_module_weights(self)