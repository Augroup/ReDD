import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', type=str, required=True,help="Input")
parser.add_argument('-o','--output', type=str, required=True,help="Output")
parser.add_argument('-c','--count', type=int, default=100,required=False,help="Number of reads each isoform will keep")
args = parser.parse_args()
df=pd.read_csv(args.input,sep='\t',header=None)
df.columns = ['ID','Isoform']
def sample_or_keep(rows,count):
    if len(rows) > count:
        return rows.sample(n=count)
    else:
        return rows
df.groupby('Isoform').apply(lambda rows:sample_or_keep(rows,args.count))['ID'].to_csv(args.output,sep='\t',index=False,header=False)