
import os
from tqdm import tqdm
import sys
import time
import datetime
import argparse

"""description: filter all qualified sequences under a certain sequence length."""

def initialization_parameters():

    parser = argparse.ArgumentParser()
    parser.add_argument('-seq_len', type=int, required=True,
                        help='Specifies the length of the sequence.')

    parser.add_argument('-output', type=str, required=True,
                        help='The output file for storing qualified sequences.')

    args = parser.parse_args()
    return args


def cal_GC_content(base_seq):
    rate = (base_seq.count('G') + base_seq.count('C')) / len(base_seq)

    if rate >=0.4 and rate <= 0.6:
        return 1

    else:
        return 0

def check_reapte_triples(base_seq):
    if 'AAA' in base_seq or 'TTT' in base_seq or 'CCC' in base_seq or 'GGG' in base_seq:
        return 1
    else:
        return 0

def check_GGC(base_seq):
    if 'GGC' in base_seq:
        return 1
    else:
        return 0

def get_3mer(base_seq):
    return set([base_seq[i:i+3] for i in range(0,len(base_seq)-2)])

def check_selfComple(base_seq):
    selfComple_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    selfComple_seq = ''
    for item in base_seq:
        selfComple_seq += selfComple_dict[item]
    
    three_mer_set1 = get_3mer(base_seq)
    three_mer_set2 = get_3mer(selfComple_seq)

    # print(three_mer_set1, three_mer_set2)

    if len(three_mer_set1.intersection(three_mer_set2)) != 0:
        return 1
    else:
        return 0

def filter(base_barcode_list):
    """description: filter out sequences that do not meet barcode requirements"""
    res_list = []
    for item in base_barcode_list:
        GC_rate = cal_GC_content(item)

def add_str(base_str_list, char_list):
    res_list = []
    for char in char_list:
        for base_str in base_str_list:

            temp_str = base_str + char
            res_list.append(temp_str)
            
    return res_list
    
def generateAllseq(DNA_len, base_char_list = ["A", "T", "C", "G"]):
    """description: get all DNA sequences under a certain length."""
    
    if len(base_char_list[0]) == DNA_len:
        
        return base_char_list
    
    else:
        base_char_list = add_str(base_char_list, ["A", "T", "C", "G"])
        return generateAllseq(DNA_len, base_char_list)

def main():
    print('script name: %s'%sys.argv[0])
    print('script start time: %s\n'%datetime.datetime.now())
    args = initialization_parameters()
    start = time.time()
    res = generateAllseq(args.seq_len, base_char_list = ["A", "T", "C", "G"])
    new_res = []
    _append = new_res.append
    for temp_str in res:
        if cal_GC_content(temp_str) == 1 and check_reapte_triples(temp_str) == 0 and check_GGC(temp_str) == 0 and check_selfComple(temp_str) == 0:
            _append(temp_str)

    file = open(args.output, 'w')
    for i in range(0,len(new_res)):
        file.write('>%d\n'%i)
        file.write(new_res[i])
        file.write('\n')
    file.close()
    
    end = time.time()
    print('Note:')
    print('    When the sequence length is %d: '%args.seq_len)
    print('    The number of qualified sequences is: %d, the run time is %s.\n'%(len(new_res), end-start))
    print('script start time: %s'%datetime.datetime.now())
        
if __name__ == "__main__":
    # test = check_selfComple("GGGCCC")
    # print(test)
    main()
        
        
