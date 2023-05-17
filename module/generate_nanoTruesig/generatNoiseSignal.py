import sys
import os

if 'module' not in sys.path:
    sys.path.append('module')

from poremodel_util import *
import numpy as np
import time
import datetime
import argparse


###########################################################################################################################
#--------------- step 0: official kmer pore model ---------#
#-> load plain text file for official pore model
#   data structure is a dict with 4096 dims (6mer 4^6)
#   data.keys()   :   kmer
#   data.values() :   (mean, vari)
#-> the file 'template_median68pA.model' could be downloaded from
#   https://github.com/nanoporetech/kmer_models/blob/master/r9.4_180mv_450bps_6mer/template_median68pA.model
def load_official_poremodel(input_file):
    model_data = np.genfromtxt(input_file, delimiter="\t", dtype=None, comments='#', names=True)
    model_dict = dict([(x[0].decode('utf-8'), (x[1], x[2])) for x in model_data])
    return model_dict

def sequence_official_poremodel(sequence, kmer_poremodel):
    k=len(list(kmer_poremodel.keys())[0])
    length=len(sequence)
    # check sequence length
    if length < k:
        # Assign mean and std value
        kmer_means=list()
        kmer_stdvs=list()
        [kmer_means.extend((float(90.2083),)*1) for i in range(length)]
        [kmer_stdvs.extend((float(2.0),)*1) for i in range(length)]
    else:
        # Divide sequence into kmers
        kmers = [sequence[i:i + k] for i in range(0, length - k + 1)]
        # Assign mean and std value
        kmer_means, kmer_stdvs = zip(*[kmer_poremodel[kmer] for kmer in kmers])
        kmer_means=list(kmer_means)
        kmer_stdvs=list(kmer_stdvs)
        # Append tail
        #[kmer_means.extend((float(90.2083),)*1) for i in range(k-1)]
        #[kmer_stdvs.extend((float(2.0),)*1) for i in range(k-1)]
    # return
    kmer_means = np.array(kmer_means)
    kmer_stdvs = np.array(kmer_stdvs)
    return kmer_means,kmer_stdvs


#----------- main program: sequence to raw signal --------------#
# default parameters:
#     repeat_alpha=0.1
#     repeat_more=1
#     event_std=1.0
#     filter_freq=850
#     noise_std=1.5
def sequence_to_true_signal(input_part, output_folder, perfect=False, p_len=1,
    repeat_alpha=0.1, repeat_more=1, event_std=1.0, filter_freq=950, noise_std=1.0, sigroot='signal',seed=0):
    #--- unzip input args ---#
    sequence = input_part[0]
    seq_name = input_part[1]
    #--- get kmer signal ----#
    kmer_poremodel = load_official_poremodel('module/generate_nanoTruesig/model_data/template_median68pA.model')
    mean_result, std_result = sequence_official_poremodel(sequence, kmer_poremodel)
    #--- kmer simulator -----#
    if perfect:
        final_result, final_ali = repeat_k_time(p_len, mean_result)
    else:
        #-> 1. repeat N times
        indep_result, final_ali, event_idx = repeat_n_time(repeat_alpha, mean_result,
            repeat_more, seed=seed)
        np.random.seed(seed)
        event_std = np.random.uniform(-1*event_std*std_result[event_idx], event_std*std_result[event_idx])
        final_result = mean_result[event_idx] + event_std
        #-> 2. low pass filter
        if filter_freq>0:
            h,h_start,N = low_pass_filter(4000.0, filter_freq, 40.0)
            final_result = np.convolve(final_result,h)[h_start+1:-(N-h_start-1)+1]
        #-> 3. add gauss noise
        if noise_std>0:
            final_result = final_result + add_noise(noise_std,
                len(final_result), seed=seed)
    #--- make integer -------#
    final_result = np.array(final_result)
    final_result = np.array(list(map(int, 5.7*final_result+14)))
    #print(final_result)
    write_output(final_result, output_folder + '/' + sigroot + '_{}.txt'.format(seq_name))

def get_parameters():
    """This script generates noiseless nanopore signals based on DNA sequences"""

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True,
                        help='input file, it is a fasta file.')

    parser.add_argument('-o', type=str, required=True,
                        help='output folder, the output signal in here.')

    parser.add_argument('-introduction', type=str, required=False,
                        help='This script generates noiseless nanopore signals based on DNA sequences')

    args = parser.parse_args()

    return args


def main():
    print('script name: %s' % sys.argv[0])
    print('script start time: %s\n' % datetime.datetime.now())
    time_start = time.time()
    args = get_parameters()
    seq_list = get_seq_list(args.i)
    id_list = get_id_list(args.i)
    zip_id_seq = list(zip(seq_list, id_list))
    isExists_out = os.path.exists(args.o)
    if not isExists_out:
        os.makedirs(args.o)

    for seq in zip_id_seq:
        sequence_to_true_signal(seq, output_folder=args.o, sigroot='timeSeries')
        #sequence_to_true_signal(in_list[0], output_folder, filter_freq=950, noise_std=1.0)

    time_end = time.time()
    time_sum = time_end - time_start

    print('Note:')
    print('    sequence fasta path: %s\n' % args.i)
    print('    output path: %s\n' % args.o)
    print('    run time: %f' % time_sum)

    print('script start time: %s' % datetime.datetime.now())

###########################################################################################################################


if __name__ == "__main__":
    main()
