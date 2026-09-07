#ifndef EVENTALIGN_REDD_INCL
#define EVENTALIGN_REDD_INCL

#include "config.hpp"
#include "pore_model.hpp"
#include "intervals.hpp"
#include <algorithm>
#include <vector>
#include <string>

#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>

#include <cmath>
#include <cstddef>
#include <limits>
#include <unordered_map>


//template <KmerLen K, typename KmerType=typename std::conditional<(K < 8), u16, u32>::type>
//struct Eventalign : public DataFrame<i64, KmerType, i32, i32, float> {
//    using Super = DataFrame<i64, KmerType, i32, i32, float>;
//
//    template <size_t I>
//    using ColType = Super::typename ColType<I>;// = typename Super::ColType<I>;
//
//    static constexpr typename Super::NameArray names = {"pac", "model_kmer", "event_index", "standardized_level"}; 
//    ColType<0> &pac = std::get<0>(Super::data_); 
//    ColType<1> &model_kmer = std::get<1>(Super::data_);                      
//    ColType<2> &event_index = std::get<2>(Super::data_);                      
//    ColType<3> &standardized_level = std::get<3>(Super::data_);
//    using Super::DataFrame;                              
//};                
//constexpr Eventalign::NameArray Eventalign::names;
typedef struct {
    std::vector<std::vector<std::vector<float>>> X;
    std::vector<std::vector<int>> y_ref;
    std::vector<std::vector<int>> y_call; 
    std::vector<std::string> info;
    std::vector<std::vector<std::vector<float>>> X_candidate;
    std::vector<std::vector<int>> y_ref_candidate;
    std::vector<std::vector<int>> y_call_candidate; 
    std::vector<std::string> info_candidate;
} redd_data_t;
static const char* complement_dna = "TGCA";
static const int rank_dna[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static inline int get_encoding(char base) {
    if (base == 'A') { 
        return 1;
    } else if (base == 'C') {
        return 2;
    } else if (base == 'G') {
        return 3;
    } else if (base == 'T') {
        return 4;
    } else if (base == 'N') {
        return 0;
    } else {
        return 0;
    }
}
typedef struct {
    u_int64_t ref_pos;
    char ref_base;
    char read_base;
    size_t event_idx_start;
    size_t event_idx_end;
    float mean;
    float stdev;
    float skewness;
    float kurtosis;
    size_t length;
} redd_feature_t;
// Helper function to compute median for a float array
float get_median_pa(float* arr, int64_t size) {
    // Create a copy to avoid modifying the original array
    float* temp = new float[size];
    memcpy(temp, arr, size * sizeof(float));
    std::sort(temp, temp + size);

    float med;
    if (size % 2 == 0)
        med = (temp[size/2 - 1] + temp[size/2]) / 2.0f;
    else
        med = temp[size/2];

    delete[] temp;
    return med;
}
// Compute MAD for a raw array of floats
float get_median_mad_pa(float* arr, int64_t size) {
    float med = get_median_pa(arr, size);
    // Compute absolute deviations
    float* deviations = new float[size];
    for (size_t i = 0; i < size; ++i)
        deviations[i] = std::fabs(arr[i] - med);

    float mad_value = get_median_pa(deviations, size);
    delete[] deviations;
    return mad_value;
}
//void write_eventalign(Eventalign<K> df, ProcessedRead read, i32 read_idx, bool fwd) {
std::vector<float> get_scaled_samples_redd(float *samples, uint64_t start_idx, uint64_t end_idx,float raw_signal_median,float raw_signal_mad){
    std::vector<float> out;
    for(uint64_t i = start_idx; i < end_idx; ++i) {
        double s = samples[i];
        double scaled_s = s - raw_signal_median;
        scaled_s /= raw_signal_mad;
        // scaled_s = s;
        if (scaled_s < -5.0){
            scaled_s = -5.0;
        }
        if (scaled_s > 5.0){
            scaled_s = 5.0;
        }
        out.push_back(scaled_s);
    }
    return out;
}
// Compute stat
void get_raw_signals_feature(const std::vector<float>& data,float& mean, size_t& length,
    float& stdev, float& skewness, float& kurtosis) {
    length = data.size();
    if (length == 0) {
        // mean = stdev = skewness = kurtosis = std::numeric_limits<float>::quiet_NaN();
        mean = stdev = skewness = kurtosis = 0;
        return;
    }
    // Step 1: Compute mean
    double sum = 0.0;
    for(size_t i = 0; i < length; ++i)
        sum += data[i];
    mean = static_cast<float>(sum / length);

    // Step 2: Central moments
    double m2 = 0.0, m3 = 0.0, m4 = 0.0;
    for(size_t i = 0; i < length; ++i) {
        double diff = data[i] - mean;
        m2 += diff * diff;
        m3 += diff * diff * diff;
        m4 += diff * diff * diff * diff;
    }

    // Variance with Bessel's correction
    double var = (length > 1) ? m2 / (length - 1) : 0.0;
    stdev = static_cast<float>(std::sqrt(var));

    // Standardized moments (sample estimators)
    if(length > 2 && stdev > 0) {
        skewness = static_cast<float>(
            (length * m3) / ((length - 1) * (length - 2) * std::pow(stdev, 3))
        );
    } else {
        skewness = 0;
    }

    if(length > 3 && stdev > 0) {
        double n = (double)length;
        double denom = (n - 1) * (n - 2) * (n - 3);
        double g2 = (n * (n + 1) * m4 - 3 * m2 * m2 * (n - 1))
                    / (denom * std::pow(stdev, 4));
        kurtosis = static_cast<float>(g2);
    } else {
        kurtosis = 0;
    }

}
template <typename ModelType>
redd_data_t write_eventalign_redd_new(Alignment<ModelType> &aln, bool write_name, bool signal_index, py::array_t<float> signal_np,py::dict _redd_candidate_ratio_map,int redd_window_size) {
    // init_redd_hdf5_file(redd_window_size,output_file);
    auto signal = PyArray<float>(signal_np);
    auto coord = aln.seq.coord;
    auto &model = aln.seq.model;
    auto sample_rate = model.PRMS.sample_rate;

    float * signal_arr = &signal[0];
    int64_t len_raw_signal = signal.size();
    // for (size_t i = 0; i < len_raw_signal; ++i){
    //     signal_arr[i] = signal_arr[i] * 23.1105 + 92.5619;
    // }
    float raw_signal_median = get_median_pa(signal_arr,len_raw_signal);
    float raw_signal_mad = get_median_mad_pa(signal_arr,len_raw_signal);
    // int n_collapse = 1;

    size_t event_idx_start = -1;
    size_t event_idx_end = -1;
    char strand = coord.fwd() ? '+' : '-';
    auto read_name = aln.read_id;
    auto ref_name = coord.name;


    int half_redd_window_size = redd_window_size/2;
    std::unordered_map<std::string, std::unordered_map<u_int64_t, float>> redd_candidate_ratio_map = _redd_candidate_ratio_map.cast<std::unordered_map<std::string, std::unordered_map<u_int64_t, float>>>();



    std::vector<redd_feature_t> redd_feature_vec;
    redd_data_t single_read_redd_data;

    

    for (int i = 0-redd_window_size; i < (int)aln.dtw.size()+redd_window_size; i++) {
        if ((i < 0) || (i >= aln.dtw.size())){
            // n_collapse = 1;
            redd_feature_t redd_feature = {0,'N','N',0,0,0,0,0,0,0};
            redd_feature_vec.push_back(redd_feature);
        } else {
            auto &samps = aln.dtw.samples.coords[i];
            if (!samps.is_valid()) {
                continue;
            }
            auto ref = aln.seq.mpos[i];
            if (ref < 0) ref = -ref-1;
            u_int64_t ref_position = ref - model.PRMS.shift;
            auto kmer = aln.seq.kmer[i], model_kmer = kmer;
            if (model.PRMS.reverse) model_kmer = model.kmer_rev(kmer);
            auto ref_kmer = model_kmer;
            if (!coord.fwd()) ref_kmer = model.kmer_revcomp(ref_kmer);
            std::string ref_kmer_seq = model.kmer_to_str(ref_kmer);
            std::string model_kmer_seq = model.kmer_to_str(model_kmer);

            event_idx_start = i;
            event_idx_end = i;
            uint64_t start_idx = samps.start+1; //inclusive
            uint64_t end_idx = samps.end; //non-inclusive
            std::vector<float> samples = get_scaled_samples_redd(signal_arr, start_idx, end_idx,raw_signal_median,raw_signal_mad);
            // n_collapse = 1;
            // while (i + n_collapse < aln.dtw.size() && aln.seq.mpos[i] ==  aln.seq.mpos[i+n_collapse]){
            //     auto &new_samps = aln.dtw.samples.coords[i+n_collapse];
            //     uint64_t new_start_idx = new_samps.start+1; //inclusive
            //     uint64_t new_end_idx = new_samps.end; //non-inclusive
            //     std::vector<float> new_samples = get_scaled_samples_redd(signal_arr, new_start_idx, new_end_idx,raw_signal_median,raw_signal_mad);
            //     samples.insert(samples.end(), new_samples.begin(), new_samples.end());
            //     // if(strcmp(ea.model_kmer,alignments[i+n_collapse].model_kmer)!=0){ //TODO: NNNN kmers must be handled
            //     //     fprintf(stderr, "model kmer does not match! %s vs %s\n",ea.model_kmer,alignments[i+n_collapse].model_kmer);
            //     // }
            //     event_idx_end = i+n_collapse;
            //     n_collapse++;
                
            // }
            float mean = 0.0;
            float stdev = 0.0;
            float skewness = 0.0;
            float kurtosis = 0.0;
            size_t length = 0;
            get_raw_signals_feature(samples,mean, length,stdev, skewness, kurtosis);
            redd_feature_t redd_feature;
            if (coord.fwd()){
                redd_feature = {ref_position,ref_kmer_seq[0],model_kmer_seq[0],event_idx_start,event_idx_end,mean,stdev,skewness,kurtosis,length};
            } else {
                // if map to the RC strand, make the ref base RC instead of read
                redd_feature = {ref_position,complement_dna[rank_dna[(int)ref_kmer_seq[0]]],model_kmer_seq[model_kmer_seq.length()-1],event_idx_start,event_idx_end,mean,stdev,skewness,kurtosis,length};
            }
            if (!redd_feature_vec.empty()){
                // check whether consecutive with previous event
                redd_feature_t prev_base_redd_feature = redd_feature_vec.back();
                if ((prev_base_redd_feature.event_idx_end + 1  != redd_feature.event_idx_start && coord.fwd()) && (prev_base_redd_feature.event_idx_end - 1  != redd_feature.event_idx_start && !coord.fwd())){
                    redd_feature_vec.clear();
                }     
            };
            redd_feature_vec.push_back(redd_feature);
            }

            if (redd_feature_vec.size() >= redd_window_size){
                if (redd_feature_vec[redd_feature_vec.size()-1-half_redd_window_size].ref_base == 'A'){
                    redd_feature_t center_feature = redd_feature_vec[redd_feature_vec.size()-1-half_redd_window_size];
                    float ratio = -1;
                    if (!redd_candidate_ratio_map.empty()){
                        if (redd_candidate_ratio_map.find(ref_name) != redd_candidate_ratio_map.end()){

                            if (redd_candidate_ratio_map[ref_name].find(center_feature.ref_pos) != redd_candidate_ratio_map[ref_name].end()){
                                ratio = redd_candidate_ratio_map[ref_name][center_feature.ref_pos];
                            }
                        }
                    }
                    bool is_candidate = false;
                    std::string info_str;
                    std::stringstream ss;
                    if (ratio != -1){
                        is_candidate = true;
                        ss << "I:" << ratio <<"\t" << read_name << "\t" << ref_name << "\t" << center_feature.ref_pos + 1 << "\t" << strand;
                    } else {
                        is_candidate = false;
                        ss << "I\t" << read_name << "\t" << ref_name << "\t" << center_feature.ref_pos + 1 << "\t" << strand;
                    }
                    info_str = ss.str();
                    if (is_candidate){
                        single_read_redd_data.X_candidate.push_back(std::vector<std::vector<float>>());
                        single_read_redd_data.y_call_candidate.push_back(std::vector<int>());
                        single_read_redd_data.y_ref_candidate.push_back(std::vector<int>());
                        single_read_redd_data.info_candidate.push_back(info_str);
                    } else {
                        single_read_redd_data.X.push_back(std::vector<std::vector<float>>());
                        single_read_redd_data.y_call.push_back(std::vector<int>());
                        single_read_redd_data.y_ref.push_back(std::vector<int>());
                        single_read_redd_data.info.push_back(info_str);
                    }

                    for (int feature_index = redd_feature_vec.size() - redd_window_size; feature_index < redd_feature_vec.size();feature_index+=1 ){

                        redd_feature_t feature = redd_feature_vec[feature_index];
                        if (is_candidate){
                            single_read_redd_data.X_candidate.back().push_back({feature.mean,feature.stdev,(float)feature.length,feature.skewness,feature.kurtosis});
                            single_read_redd_data.y_ref_candidate.back().push_back(get_encoding(feature.ref_base));
                            single_read_redd_data.y_call_candidate.back().push_back(get_encoding(feature.read_base));
                        } else {
                            single_read_redd_data.X.back().push_back({feature.mean,feature.stdev,(float)feature.length,feature.skewness,feature.kurtosis});
                            single_read_redd_data.y_ref.back().push_back(get_encoding(feature.ref_base));
                            single_read_redd_data.y_call.back().push_back(get_encoding(feature.read_base));

                        }

                    };


                }
            }
            if (redd_feature_vec.size() > 2 * redd_window_size){
                redd_feature_vec.erase(redd_feature_vec.begin(), redd_feature_vec.begin() + redd_window_size);
            }
    }

        
        // //if (aln.dtw.current[i] == ValArray<float>::NA) continue;
        // //auto &evt = read.events[i];
        // auto kmer = aln.seq.kmer[i], model_kmer = kmer;
        // if (model.PRMS.reverse) model_kmer = model.kmer_rev(kmer);
        // auto ref_kmer = model_kmer;
        // if (!coord.fwd()) ref_kmer = model.kmer_revcomp(ref_kmer);

        // auto ref = aln.seq.mpos[i];
        // if (ref < 0) ref = -ref-1;


        // ss << coord.name << "\t"
        //    << ref - model.PRMS.shift << "\t"
        //    << model.kmer_to_str(ref_kmer) << "\t";

        // if (write_name) {
        //    ss << aln.read_id;
        // } else {
        //    ss << aln.id;
        // }
        // ss << "\tt" << "\t"
        //    << i << "\t"
        //    << model.current.norm_to_pa(aln.dtw.current[i]) << "\t"
        //    << (model.current.norm_to_pa_sd(aln.dtw.current_sd[i])) << "\t"
        //    << (samps.length() / sample_rate) << "\t"
        //    << model.kmer_to_str(model_kmer) << "\t"
        //    << model.current.norm_to_pa(model.current.mean[kmer]) << "\t"
        //    << model.current.norm_to_pa_sd(model.current.stdv[kmer]) << "\t"
        //    << aln.dtw.current[i];

        // if (signal_index) {
        //    ss << "\t" << samps.start << "\t" << samps.end; //<< "\n";
        // }

        // if (signal.size() > 0) {
        //     ss << "\t" << signal[samps.start];
        //     for (size_t j = samps.start+1; j < samps.end; j++) {
        //         ss << "," << signal[j];
        //     }
        // }

        // ss << "\n";
    // }
    // ss << "A\n";
    // return ss.str();
    // if (redd_candidate_ratio_map.empty()){
    //     append_arr_to_dataset(redd_window_size,hdf5_output_file_prefix+".hdf5",single_read_redd_data.X,single_read_redd_data.y_ref,single_read_redd_data.y_call, single_read_redd_data.info);
    // } else {
    //     append_arr_to_dataset(redd_window_size,hdf5_output_file_prefix+".candidate.hdf5",single_read_redd_data.X_candidate,single_read_redd_data.y_ref_candidate,single_read_redd_data.y_call_candidate, single_read_redd_data.info_candidate);
    //     append_arr_to_dataset(redd_window_size,hdf5_output_file_prefix+".noncandidate.hdf5",single_read_redd_data.X,single_read_redd_data.y_ref,single_read_redd_data.y_call, single_read_redd_data.info);
    // }
    return single_read_redd_data;
};

#endif
