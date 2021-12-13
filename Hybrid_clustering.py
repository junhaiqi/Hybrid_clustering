# coding:utf-8
import os
import sys
import math
from numpy.lib.function_base import average
# from pygal.util import merge
from scipy.special import comb
import random
import threading
# from datetime import datetime
from time import time

# Get the number of consensus sequences
def Get_consensus_seq_count(DNA_len):
    DNA_len = int(DNA_len)
    singleBase_probability = math.pow(0.9999,1/DNA_len)
    print(singleBase_probability)
    correct_probability = 0.8
    error_probability = 1-correct_probability
    n = 30
    base_line = int(n/4)
    fail_probability = 0
    # Less than n/4
    for i in range(0,base_line):
        single_correct_probability = math.pow(correct_probability,i)*comb(n,n-i)* \
        math.pow(error_probability,n-i)
        fail_probability += single_correct_probability
    # large than n/4 and less than n/2
    for j in range(base_line, 2*base_line+1):
        local_error_probability = comb(n-j,j+1)*3*math.pow(error_probability/3,j+1)* \
        math.pow(correct_probability,j)*comb(n,j)*math.pow(error_probability,n-2*j-1)                                                      
        print(local_error_probability)
        fail_probability+=local_error_probability
    return 1-fail_probability

def Get_DNA_count_Fromfasta(fasta_file_name):
    fasta_file = open(fasta_file_name,'r')
    t = 0
    for line in fasta_file:
        if '>' in line:
            t += 1
    fasta_file.close()
    return t

def trans(i):
    f = open("%d.clstr" %i)
    c = open("trans_%d.txt"%i,"w")
    for line in f:
        if "Cluster" in line:
            c.write("\n")
        else:
            line = line.strip("\n")
            line = line.split(">")
            line = line[1].split(".")
            c.write(line[0])
            c.write("	")

def TransList(i):
    trans(i)
    b = []
    c = open("trans_%d.txt"%i,"r")
    for line in c.readlines()[1:]:
        line = line.strip("\n")
        line = line.strip("\t")
        line = [int(a) for a in line.split('\t')]
        b.append(line)
    return b

def define_good_index(clustering_list):
    clustering_len_list = []
    for item in clustering_list:
        clustering_len_list.append(len(item))
        
    clustering_len_list.sort()
    index = clustering_len_list[int(len(clustering_list)/1.6)]
    return index

def write_info_file(T,threshold,K,goodIndex):
    th_good_count_file = open('th_goodcount_file.txt','w+')
    th_good_count_file.write(str(T))
    th_good_count_file.write('\n')
    th_good_count_file.write(str(threshold))
    th_good_count_file.write('\n')
    th_good_count_file.write(str(K))
    th_good_count_file.write('\n')
    th_good_count_file.write(str(goodIndex-1))
    th_good_count_file.close()

def Get_goodcluster_list():
    goodclusterfile = open('goodclusterfile.txt','r')
    goodclusterlist = []
    for line in goodclusterfile:
        line = line.strip('\n')
        line = line.split(' ')
        temp_list = []
        for item in line:
            if item != '':
                temp_list.append(int(item) + 1)
        goodclusterlist.append(temp_list)
    return goodclusterlist

def Get_refinecluster_len(refinegoodcluster):
    sum = 0
    for item in refinegoodcluster:
        sum += len(item)
    return sum

def write_clusteringresult_toFile(clustering_result,res_name):
    res_file = open(res_name,'w')
    for i in range(0,len(clustering_result)):
        for j in range(0,len(clustering_result[i])):
            res_file.write('%s '%str(clustering_result[i][j]))
        res_file.write('\n')
    res_file.close()
    
def fromFastaGetATCG(fasta_name):
    fasta_file = open(fasta_name,'r')
    fasta_list = []
    for line in fasta_file:
        if '>' not in line:
            line = line.strip('\n')
            fasta_list.append(line)
    fasta_file.close()
    return fasta_list

def Get_less50_list(cluster_list):
    temp_list = []
    for item in cluster_list:
        if len(item) > 50:
            temp_list.append(item[0:50])
        else:
            temp_list.append(item)
    return temp_list

def Get_label_set(cluster_list):
    temp_list = []
    for item in cluster_list:
        temp_list = temp_list + item
    return set(temp_list)

def Get_all_label_set(fasta_name):
    fasta_file = open(fasta_name,'r')
    all_label_list = []
    for line in fasta_file:
        if '>' in line:
            line = line.strip('\n')
            line = line.split('>')
            all_label_list.append(int(line[1]))
    fasta_file.close()
    return set(all_label_list)

def Get_TwoSet_difference(set_1,set_2):
    return list(set_1.difference(set_2))

def write_Refine_info(cluster_list,ReadyToSortlist):
    ReadyToSortfile = open('ReadyToSortfile.txt','w+')
    ReadyToSortfile.write(str(len(cluster_list)))
    ReadyToSortfile.write('\n')
    ReadyToSortfile.write(str(T))
    ReadyToSortfile.write('\n')
    for item in ReadyToSortlist:
        ReadyToSortfile.write(str(item-1))
        ReadyToSortfile.write('\n')
    ReadyToSortfile.close()
    
def write_mkdir_file_info(s_1,s_2):
    argv_file = open('argv_file.txt','w')
    argv_file.write(s_1)
    argv_file.write('\n')
    argv_file.write(s_2)
    argv_file.close()
    
def getminvalue(dislist):
    minvalue = dislist[0]
    for i in range(1,len(dislist)):
        if dislist[i] < minvalue:
            minvalue = dislist[i]
    return minvalue
        
def findminposition(minvalue,dislist):
    for i in range(0,len(dislist)):
        if minvalue == dislist[i]:
            return i
# threshold        
def Get_maxlen_list(a):
    indicator_list = []
    index = 0
    for item in a:
        if len(item) > index:
            indicator_list = item
            index = len(item)
    return indicator_list

def getdata_2(signal_file1,signal_file2):
    A = 0
    result = os.popen("./dtw_zscore3 %s %s"%(signal_file1,signal_file2))
    context = result.read()
    for line in context.splitlines():
        c = line
        result.close()
        A = c
    return (float(A))

def Get_not0_count(list):
    t = 0

    for item in list:
        if item != 0:
            t += 1
    return t

def Get_Signal_list(i,Signal_folder):
    test_file = open('%s/signal_%d.txt'%(Signal_folder,i))
    signal_list = []
    for line in test_file:
        signal_list.append(float(line))
    test_file.close()
    return signal_list

def Get_merge_threshold(OldCluster,Signal_folder,fasta_file):
    
    max_list = Get_maxlen_list(OldCluster)
    fasta_list = fromFastaGetATCG(fasta_file)
    
    fasta_len = 0
    fasta_len = 0
    no_0_fasta_len = 0
    for item in fasta_list:
        if len(item) != 0:
            fasta_len+=len(item)
            no_0_fasta_len+=1
    
    average_1 = fasta_len/no_0_fasta_len
    
    len_max_list = len(max_list)
    signallength = 0
    no_0_sig_len = 0
    for i in range(0,len_max_list):
        sig = Get_Signal_list(max_list[i] - 1,Signal_folder)
        if len(sig) != 0:
            signallength += len(sig)
            no_0_sig_len += 1
            
    average_2 = (signallength/no_0_sig_len)/8
    c = int((average_1 + average_2)/2) - 5
    
    temp_list = []
    if len(max_list) > 10:
       RandomCheck1 = random.sample(max_list,10)
       for i in range(0,len(RandomCheck1)-1):
           for j in range(i+1,len(RandomCheck1)):
               s1 = Signal_folder + '/' + 'signal_' + str(RandomCheck1[i]) + '.txt'
               s2 = Signal_folder + '/' + 'signal_' + str(RandomCheck1[j]) + '.txt'
               temp_list.append(getdata_2(s1,s2))

    else:
        for i in range(0,len(max_list)-1):
            for j in range(i+1,len(max_list)):
                s1 = Signal_folder + '/' + 'signal_' + str(max_list[i]) + '.txt'
                s2 = Signal_folder + '/' + 'signal_' + str(max_list[j]) + '.txt'
                temp_list.append(getdata_2(s1,s2))

    sum = 0
    for item in temp_list:
        sum += item
    sum = sum - max(temp_list)*2
    len_not0 = Get_not0_count(temp_list) - 2

    t = int(sum/len_not0)
    threshold = t + c
    
    return threshold

def mafft(i):
    os.system("mafft --globalpair --maxiterate 16 --inputorder \
              Goodcluster/%d.fasta >aligned/%d.fasta"%(i,i))

def Multithreading_mafft(thread_num,start_file_id):
    threads = []
    for i in range(start_file_id,thread_num + start_file_id):
        threads.append(threading.Thread(target=mafft,args=(i,)))

    for t in threads:
        t.setDaemon(True)
        t.start()

    for t in threads:
        t.join()
        
def consensus(i):
    os.system("./consensus.py -i aligned/%d.fasta -o %d -c 1"%(i,i))
    # print('The consistent sequence of aligned/%d.fasta is generated!'%i)
    
def Multithreading_consensus(thread_num,start_file_id):
    threads = []
    for i in range(start_file_id,thread_num + start_file_id):
        threads.append(threading.Thread(target=consensus,args=(i,)))
        
    for t in threads:
        t.setDaemon(True)
        t.start()

    for t in threads:
        t.join()
        
def get_consensussignal(i):
    os.system("./get_consensus_signal.sh -i Consensus_Outputs/%d_consensus.fasta \
              -o consensus_sim/%d_sim"%(i,i))
    
def Multithreading_get_consensussignal(thread_num,start_file_id):
    threads = []
    for i in range(start_file_id,thread_num + start_file_id):
        threads.append(threading.Thread(target=get_consensussignal,args=(i,)))
    
    for t in threads:
        t.setDaemon(True)
        t.start()

    for t in threads:
        t.join()

def process_additional_output(clustering_count):
    os.system('rm argv_file.txt')
    os.system('rm th_goodcount_file.txt')
    os.system('rm ReadyToSortfile_1.txt')
    os.system('rm ReadyToSortfile.txt')
    os.system('rm trans_*.txt')
    os.system('rm -r consensus_*')
    os.system('rm -r Consensus_*')
    os.system('rm -r Good*')
    os.system('rm *luster*.txt')
    os.system('rm -r aligned')
    os.system('rm OnetoNdisfile.txt')
    os.system('rm %d*'%clustering_count)   

def cp_consensus_signal(consensus_list):
    for i in range(0,len(consensus_list)):
        os.system("cp consensus_sim/%d_sim/signal/signal_0.txt consensus_sig_%s/consensus_sig_%d.txt"%(i,sys.argv[1],i))
        print('The %dth consensus signal has been generated! There are %d consensus signals.'%(i+1,len(consensus_list)))
        print('--------------------------------------------------------------------------------')
        print('################################################################################')

# start = datetime.now()
start = time()
if len(sys.argv) != 4:
    print('usage: python Hybrid_clustering.py *.fasta Signal_folder out.txt')
    print('*.fasta ---FASTA file to be clustered.')
    print('Signal_folder ---Nanopore signal files corresponding to DNA. Naming rules: signal_0.txt,signal_1.txt...')
    print('out.txt ---Result output file.')
    # print('thread_number ---Threads used to form a consensus sequence')
else:
    # initial clustering
    T = Get_DNA_count_Fromfasta(sys.argv[1])
    os.system('./initial_clustering -i %s -o %d -c %f'%(sys.argv[1],T,0.85))
    OldCluster = TransList(T)
    write_clusteringresult_toFile(OldCluster,'initial_%s'%sys.argv[3])
    print('The merge threshold is being determined...')
    threshold = Get_merge_threshold(OldCluster,sys.argv[2],sys.argv[1])
    # threshold = 86
    print('Merge threshold determined!')
    print('Merge threshold is %f'%threshold)
    goodIndex = define_good_index(OldCluster)
    goodIndex = max(3,goodIndex)
    print('Initial clustering completed!')
    end_initial = time()
    print('good index is: %d'%goodIndex)
    K = min(3, goodIndex + 1)
    write_info_file(T,threshold,K,goodIndex)
    print('Good clusters merging...')
    write_mkdir_file_info(sys.argv[2],'consensus_sig_' + sys.argv[1])
    # start mergering
    start_merge = time()
    os.system('./mergering_cluster')
    end_good_cluster_merge = time()
    # os.system('./merge')
    refineGoodCluster = Get_goodcluster_list()
    # end mergering
    #==========================================================
    # get mergering result
    print('Outputting merge results...')
    mergering_file = open('merge_%s'%sys.argv[3],'w')
    refine_copy = []
    refine_list = []
    
    for item in refineGoodCluster:
        refine_list = refine_list + item
        
    for item in refineGoodCluster:
        refine_copy.append(item)
        
    for item in OldCluster:
        if item[0] not in refine_list:
            for item_1 in item:
                temp = []
                temp.append(item_1)
                refine_copy.append(temp)
    
    refine_copy_len = Get_refinecluster_len(refine_copy)
    
    if refine_copy_len == T:
        for item in refine_copy:
            for item_1 in item:
                mergering_file.write(str(item_1))
                mergering_file.write(' ')
            mergering_file.write('\n')
        mergering_file.close()
    
    else:
       print('Error merger result!')
    
    #==========================================================
    lenRefinecluster = Get_refinecluster_len(refineGoodCluster)
    if lenRefinecluster == T:
        print('Clustering complete! Total seq: %d'%T)
        write_clusteringresult_toFile(refineGoodCluster,sys.argv[3])
        
    # start refinement
    else:
        start_consensus = time()
        print('consensus sequences will be generated...')
        fasta_list = fromFastaGetATCG(sys.argv[1])
        os.system('mkdir Goodcluster')
        consensus_list = Get_less50_list(refineGoodCluster)
        for i in range(0,len(consensus_list)):
            goodClusterfile = open('Goodcluster/%d.fasta'%i,'w')
            t = 1
            for j in range(0,len(consensus_list[i])):
                goodClusterfile.write('>%d\n'%t)
                goodClusterfile.write(fasta_list[consensus_list[i][j]-1])
                goodClusterfile.write('\n')
                t+=1
            goodClusterfile.close()
            
        os.system('mkdir aligned')
        os.system('mkdir consensus_sim')
        os.system('mkdir consensus_sig_%s'%sys.argv[1])
        
        #Multithreaded version
        max_thread = 8
        split_stage = int(len(consensus_list)/max_thread)
        
        if split_stage == 0:
            max_thread = len(consensus_list)
            Multithreading_mafft(max_thread,0)
            Multithreading_consensus(max_thread,0)
            Multithreading_get_consensussignal(max_thread,0)
            cp_consensus_signal(consensus_list)
        else:
            addtion_num = len(consensus_list) - split_stage*8
            for i in range(0,split_stage):
                Multithreading_mafft(max_thread,i*8)
                Multithreading_consensus(max_thread,i*8)
                Multithreading_get_consensussignal(max_thread,i*8)
            
            Multithreading_mafft(addtion_num,split_stage*8)
            Multithreading_consensus(addtion_num,split_stage*8)
            Multithreading_get_consensussignal(addtion_num,split_stage*8)
            cp_consensus_signal(consensus_list)
        end_merge = time()
        all_label_set = Get_all_label_set(sys.argv[1])
        now_label_set = Get_label_set(refineGoodCluster)
        ReadyToSort = Get_TwoSet_difference(all_label_set,now_label_set)
        write_Refine_info(refineGoodCluster,ReadyToSort)
        os.system('./cudamnDTW')
        
        disfile = open('OnetoNdisfile.txt','r')
        OnetoNdislist = []
        for line in disfile:
            line = line.strip('\n')
            line = line.split(' ')
            line.remove('')
            OnetoNdislist.append(line)

        for i in range(0,len(OnetoNdislist)):
            for j in range(0,len(OnetoNdislist[i])):
                OnetoNdislist[i][j] = float(OnetoNdislist[i][j])
        Test_Readytosort = []
        for i in range(0,len(OnetoNdislist)):
            minvalue = min(OnetoNdislist[i])
            if minvalue < threshold:
                Position = OnetoNdislist[i].index(minvalue)
                refineGoodCluster[Position].append(ReadyToSort[i])
                Test_Readytosort.append(ReadyToSort[i])
                
        if Get_refinecluster_len( refineGoodCluster) == T:
             print('Clustering complete! Total seq: %d'%T)
             write_clusteringresult_toFile(refineGoodCluster,sys.argv[3])
        else:
            ReadyToSort = list(set(ReadyToSort).difference(set(Test_Readytosort)))
            ReadyToSort.sort()
            
            ReadyToSortfile_1 = open('ReadyToSortfile_1.txt','w+')
            ReadyToSortfile_1.write(str(T))
            ReadyToSortfile_1.write('\n')
            ReadyToSortfile_1.write(str(threshold))
            ReadyToSortfile_1.write('\n')
            for i in range(0,len(ReadyToSort)):
                ReadyToSortfile_1.write(str(ReadyToSort[i]))
                ReadyToSortfile_1.write('\n')
            ReadyToSortfile_1.close()
            
            os.system('./cudannDTW_cluster')
            
            file = open('gpuClusterResult.txt','r')
            refineGoodCluster_2 = []
            for line in file:
                line = line.strip('\n')
                line = line.split(' ')
                line.remove('')
                temp = []
                for item in line:
                    temp.append(int(item))
                refineGoodCluster_2.append(temp)
            file.close()
            
            for item in refineGoodCluster_2:
                refineGoodCluster.append(item)
            clustering_count = Get_refinecluster_len(refineGoodCluster)
            print('Clustering complete! Total seq: %d'%clustering_count)
            write_clusteringresult_toFile(refineGoodCluster,sys.argv[3])
            
# end = datetime.now()
    end = time()
    print('run time: %s'%str(end-start))
    print('initial clustering time: %s'%str(end_initial - start))
    print('cluster merging time: %s'%str(end_merge-start_merge))
    print('refinement time: %s'%str(end - start - (end_initial - start) - (end_merge - start_merge)))
    print('good clusters merging: %s'%str(end_good_cluster_merge - start_merge))
    print('consensus time: %s'%(end_merge - start_consensus))
    
    print('Processing additional output...')
    process_additional_output(clustering_count)
    print('Finished!')
    # os.system('rm -r Consensus_Outputs')
            
            
            
        
