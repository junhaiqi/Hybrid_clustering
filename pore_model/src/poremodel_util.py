import numpy as np
import scipy.stats as st
from random import *
import scipy.signal
import uuid
from shutil import copyfile
import os
import h5py


#--------------- step 1: load input sequence ---------------#
def get_seq_list(file_name):
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    seq_list = filter(lambda x: x!='', lines)
    seq_list = filter(lambda x: '>' not in x, seq_list)
    return seq_list

def get_id_list(file_name):
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    lines = filter(lambda x: '>' in x, lines)
    id_list = map(lambda x: x.split('|')[0][1:], lines)
    return id_list

#------ output functions ------#
def write_output(result, file_name):
    with open(file_name, 'w') as f:
        for i in result:
            temp = str(i)+'\n'
            f.write(temp)

def write_alignment(result, file_name):
    with open(file_name, 'w') as f:
        for i in result:
            temp = str(i[0]+1)+' '+str(i[1]+1)+'\n'
            f.write(temp)

def signal2fasta5(template_file, data_in, fast5_root, fast5_base):
    uid = str(uuid.uuid4())
    fast5_fn = os.path.join(fast5_root,
        fast5_base+'_'+uid+'.fast5')
    copyfile(template_file, fast5_fn)
    ##Open file
    try:
        fast5_data = h5py.File(fast5_fn, 'r+')
    except IOError:
        raise IOError, 'Error opening file. Likely a corrupted file.'

    #Get raw data
    try:
        raw_dat   = fast5_data['/Raw/Reads/'].values()[0]
        raw_attrs = raw_dat.attrs
        del raw_dat['Signal']
        raw_dat.create_dataset('Signal',data=data_in, dtype='i2', compression='gzip', compression_opts=9)  #-> with compression
        raw_attrs['duration'] = data_in.size
        raw_attrs['read_id'] = uid
    except:
        raise RuntimeError, (
            'Raw data is not stored in Raw/Reads/Read_[read#] so ' +
            'new segments cannot be identified.')
    fast5_data.close()


#---------- step 2: repeat length sample -----------#
def rep_rvs(size,a, more, seed=0):
    a = a*5
    array_1 = np.ones(int(size*(0.075-0.015*a))).astype(int)
    samples = st.alpha.rvs(3.3928495261646932+a,
        -7.6451557771999035+(2*a), 50.873948369526737,
        size=(size-int(size*(0.075-0.015*a))), random_state=seed).astype(int)
    samples = np.concatenate((samples, array_1), 0)
    samples[samples<1] = 1
    samples[samples>40] = 40
    if more == 1:
        np.random.seed(seed)
        addi = np.array(abs(np.random.normal(2,1,size))).astype(int)
        samples[samples<8] += addi[samples<8]
        np.random.shuffle(samples)
        samples[samples<8] += addi[samples<8]
    return samples

def repeat_n_time(a, result, more, seed=0):
    rep_times = rep_rvs(len(result), a, more, seed)
    out = list()
    ali = list()
    pos = 0
    for i in range(len(result)):
        k = rep_times[i]
        cur = [result[i]] * k
        out.extend(cur)
        for j in range(k):
            ali.append((pos,i))
            pos = pos + 1
    event_idx = np.repeat(np.arange(len(result)), rep_times)
    return out,ali,event_idx
    
def repeat_k_time(k, result):
    out = list()
    ali = list()
    pos = 0
    for i in range(len(result)):
        cur = [result[i]] * k
        out.extend(cur)
        for j in range(k):
            ali.append((pos,i))
            pos = pos + 1
    return out,ali


#------------- step 3: low pass filter for signal simulation -----#
#-> low pass filter
#   sampling_rate = 4000.0, cut_off_freq = 1750.0, bandwidth_freq = 40.0
def low_pass_filter(sampling_rate, cut_off_freq, bandwidth_freq):
    # Read input parameter
    fS = sampling_rate  # Sampling rate.
    fL = cut_off_freq   # Cutoff frequency.
    fb = bandwidth_freq # Bandwidth frequency

    # Generate frequency bin
    b = fb / fS
    N = int(np.ceil((4 / b)))
    if not N % 2: N += 1  # Make sure that N is odd.
    n = np.arange(N)

    # Compute sinc filter.
    h = np.sinc(2 * fL / fS * (n - (N - 1) / 2.))

    # Compute Blackman window.
    w = 0.42 - 0.5 * np.cos(2 * np.pi * n / (N - 1)) + \
        0.08 * np.cos(4 * np.pi * n / (N - 1))

    # Compute h and h_start
    h = h * w
    h /= np.sum(h)
    impulse = np.repeat(0., len(h))
    impulse[0] = 1.
    h_response = scipy.signal.lfilter(h, 1, impulse)
    h_start = np.argmax(h_response)

    # return
    return h,h_start,N


#---------- step 4: add Gaussian noise ----------#
def add_noise(std, l, seed=0):
    np.random.seed(seed)
    noise = np.random.normal(0, std, l)
    return noise


