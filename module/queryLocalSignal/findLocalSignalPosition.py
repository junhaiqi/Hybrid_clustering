import numpy as np
import sys
if 'module' not in sys.path:
    sys.path.append('module')

from dtw_semi_global import semi_global_dtw
from dtw_semi_global import semi_global_dtw_with_rescaling


def drawFigureForTest(position1, position2):

    signal1 = get_signal_file(filetxt_path = 'testData/signal_0.txt')[position1:position2]

    signal2 = get_signal_file(filetxt_path = 'testData/top_adapter_signal.txt')

    plt.plot(normalise(signal1))
    plt.plot(normalise(signal2))

    plt.savefig('testData/test.png')
    plt.show()

def drawFigureForCheckbarcodeSig(position1):

    signal1 = get_signal_file(filetxt_path = 'testData/signal_0.txt')[position1:position1 + 250]

    signal2 = get_signal_file(filetxt_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/BarcodeSignal15bp831/timeSeries_0.txt')

    plt.plot(normalise(signal1))
    plt.plot(normalise(signal2))

    plt.savefig('testData/test.png')
    plt.show()


def normalise(signal):
    if len(signal) == 0:
        return signal
    mean = np.mean(signal)
    stdev = np.std(signal)
    if stdev > 0.0:
        return (signal - mean) / stdev
    else:
        return signal - mean

def get_signal_file(filetxt_path):
    signal_list = list()
    with open(filetxt_path, 'r') as f:
        signal_file = f.readlines()
        for line in signal_file:
            signal_Value = float(line.rstrip())
            signal_list.append(signal_Value)

    signal = np.array(signal_list)
    return normalise(signal)
    # return signal

def fromLongRefFindShortQuery(ref_signal, query_signal):

    # query_signal = get_signal_file(query_signal_path)
    # ref_signal = get_signal_file(ref_siganl_path)

    normalise_query_signal = normalise(query_signal)
    normalise_ref_signal = normalise(ref_signal[0:2000])

    # topRef_signal = normalise_ref_signal[0:1000]

    position_start, position_end = \
        semi_global_dtw(normalise_ref_signal, normalise_query_signal)

    return (position_start, position_end + 10)

def fromLongRefFindShortQuery_test(ref_signal_path, query_signal_path):

    query_signal = get_signal_file(query_signal_path)
    ref_signal = get_signal_file(ref_signal_path)

    normalise_query_signal = normalise(query_signal)
    normalise_ref_signal = normalise(ref_signal[0:2000])

    # topRef_signal = normalise_ref_signal[0:1000]

    position_start, position_end = \
        semi_global_dtw(normalise_ref_signal, normalise_query_signal)

    return (position_start, position_end + 10)


def RefromLongRefFindShortQuery(ref_signal, query_signal):

    # query_signal = get_signal_file(query_signal_path)
    # ref_signal = get_signal_file(ref_siganl_path)

    normalise_query_signal = normalise(query_signal)
    normalise_ref_signal = normalise(ref_signal)

    # topRef_signal = normalise_ref_signal[0:1000]

    position_start, position_end = \
        semi_global_dtw_with_rescaling(normalise_ref_signal, normalise_query_signal)

    return (position_start, position_end + 10)

def checkFindFlankingSeqPosition(sigId = 0, sigId2 = 0, mode = 'top'):
    test1 = fromLongRefFindShortQuery_test(ref_signal_path = 
            '../../exper_15bp_randomVSOurs/data/SelectedBarcodeLongFlanking100_15bp/barcode_%d/signal/signal_%d.txt'%(sigId, sigId2) \
    , query_signal_path = '../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/%s/timeSeries_0.txt'%mode)
    
    print(test1)
    print('position1: %d, position2: %d'%(test1[0], test1[1]))
    sig1 = get_signal_file('../../exper_15bp_randomVSOurs/data/SelectedBarcodeLongFlanking100_15bp/barcode_%d/signal/signal_%d.txt'%(sigId, sigId2))[test1[0]:test1[1]]
    sig2 = get_signal_file('../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/%s/timeSeries_0.txt'%mode)
    
    plt.plot(sig1)
    plt.plot(sig2)
    plt.show()

def RefromLongRefFindShortQuery_twoEnd(ref_signal_path, bottom_signal_path, topAndAdapter_signal_path):
    query_signal1 = get_signal_file(topAndAdapter_signal_path)
    # print(len(query_signal1))
    query_signal2 = get_signal_file(bottom_signal_path)

    ref_signal = get_signal_file(ref_signal_path)

    normalise_query_signal1 = normalise(query_signal1)
    normalise_query_signal2 = normalise(query_signal2)
    normalise_ref_signal = normalise(ref_signal[0:2000])

    # topRef_signal = normalise_ref_signal[0:1000]

    position_start1, position_end1 = \
        semi_global_dtw(normalise_ref_signal, normalise_query_signal1)

    position_start2, position_end2 = \
        semi_global_dtw(normalise_ref_signal, normalise_query_signal2)

    # print(position_start1, position_end1, position_start2, position_end2)
    return (position_end1, position_start2, ref_signal[position_end1 + 60:position_start2-40])

if __name__ == "__main__":
    from tqdm import tqdm
    import os
    # t = 0
    # for i in tqdm(range(0, 100)):
    #     for j in range(0, 100):
    #         os.system('cp ../../exper_15bp_randomVSOurs/data/RandomBarcodeLongFlanking100_15bp/barcode_%d/signal/signal_%d.txt \
    #         ../../exper_15bp_randomVSOurs/data/RandomBarcodeLongFlanking100_15bp_allSig/timeSeries_%d.txt'%(i, j, t))
    #         t += 1

    # os.system('mkdir ../../exper_15bp_randomVSOurs/data/RandomBarcodeLongFlanking100_15bp_allSig_barcode')
    # for i in tqdm(range(0, 10000)):
    #     file = open('../../exper_15bp_randomVSOurs/data/RandomBarcodeLongFlanking100_15bp_allSig_barcode/timeSeries_%d.txt'%i, 'w')
    #     test = RefromLongRefFindShortQuery_twoEnd(ref_signal_path = '../../exper_15bp_randomVSOurs/data/RandomBarcodeLongFlanking100_15bp_allSig/timeSeries_%d.txt'%i, 
    #                                     bottom_signal_path = '../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/bottom/timeSeries_0.txt', 
    #                                     topAndAdapter_signal_path = '../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/adapter/timeSeries_0.txt')

    #     # print(i, len(test[2]), test[0], test[1])
    #     for item in test[2]:
    #         file.write(str(item))
    #         file.write('\n')

    #     file.close()

    # os.system('mkdir ../../exper_15bp_randomVSOurs/data/SelectedBarcodeLongFlanking100_15bp_allSig_barcode')
    # for i in tqdm(range(0, 10000)):
    #     file = open('../../exper_15bp_randomVSOurs/data/SelectedBarcodeLongFlanking100_15bp_allSig_barcode/timeSeries_%d.txt'%i, 'w')
    #     test = RefromLongRefFindShortQuery_twoEnd(ref_signal_path = '../../exper_15bp_randomVSOurs/data/SelectedBarcodeLongFlanking100_15bp_allSig/timeSeries_%d.txt'%i, 
    #                                     bottom_signal_path = '../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/bottom/timeSeries_0.txt', 
    #                                     topAndAdapter_signal_path = '../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/adapter/timeSeries_0.txt')

    #     # print(i, len(test[2]), test[0], test[1])
    #     for item in test[2]:
    #         file.write(str(item))
    #         file.write('\n')

    #     file.close()

    # sig1 = get_signal_file('../../exper_15bp_randomVSOurs/data/SelectedBarcodeLongFlanking100_15bp_allSig_barcode/timeSeries_123.txt')
    # sig1 = sig1[0:len(sig1) - 50]
    # sig2 = get_signal_file('../../exper_15bp_randomVSOurs/data/selectedLongflankingBarcode15bpWithTopFlanking/timeSeries_2.txt')
    # sig3 = get_signal_file('../../OursbarcodeSignalLongFlanking/timeSeries_0.txt')
    # # sig3 = get_signal_file('../../exper_15bp_randomVSOurs/data/SelectedBarcodeLongFlanking100_15bp_allSig_barcode/timeSeries_1.txt')
    
    # plt.plot(sig1)
    # plt.plot(sig2)
    # plt.plot(sig3)
    # plt.show()
    # if list(sig1) == list(sig2):
    #     print('good')
    # query_sig1 = get_signal_file('../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/top_and_adapter/timeSeries_0.txt')
    # test1 = fromLongRefFindShortQuery(ref_signal = sig1, query_signal = query_sig1)
    # test2 = fromLongRefFindShortQuery(ref_signal = sig2, query_signal = query_sig1)
    # print(test1, test2)

    # test3 = fromLongRefFindShortQuery_test(ref_signal_path = '../../exper_15bp_randomVSOurs/data/SelectedBarcodeLongFlanking100_15bp/barcode_0/signal/signal_6.txt', 
    # query_signal_path = '../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/top_and_adapter/timeSeries_0.txt')
    # print(test3)

    # test = RefromLongRefFindShortQuery_twoEnd(ref_signal_path = '../../exper_15bp_randomVSOurs/data/SelectedBarcodeLongFlanking100_15bp_allSig/timeSeries_6.txt', 
    #                                     bottom_signal_path = '../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/bottom/timeSeries_0.txt', 
    #                                     topAndAdapter_signal_path = '../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/top_and_adapter/timeSeries_0.txt')
    # print(test, len(test[2]))

    plt.figure(figsize=(8, 4))
    plt.subplot(2,2,1)
    trueONT12Barcode1Sig = [0.057, 0.057, 0.057, 0.057, 0.057, 0.057, 0.057, -0.346, -0.346, -0.346, -0.346, -0.346, -0.346, -0.346, -0.346, -0.346, -1.260, -1.260, -1.260, -1.260, -1.260, -1.260, -1.260, -1.260, -1.260, -1.971, -1.971, -1.971, -1.971, -1.971, -1.971, -1.971, -1.971, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.748, 0.748, 0.748, 0.748, 0.748, 0.748, -0.506, -0.506, -0.506, -0.506, -0.506, -0.506, -0.047, -0.047, -0.047, -0.047, -0.047, -0.047, -0.047, -0.047, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.264, 0.264, 0.264, 0.264, 0.264, 0.264, 1.148, 1.148, 1.148, 1.148, 1.148, 1.148, 1.148, 1.148, 0.209, 0.209, 0.209, 0.209, 0.209, 0.209, -0.299, -0.299, -0.299, -0.299, -0.299, -0.299, -0.299, -0.299, -0.693, -0.693, -0.693, -0.693, -0.693, -0.693, -0.693, -0.693, -1.315, -1.315, -1.315, -1.315, -1.315, -1.315, -1.315, -1.315, -0.440, -0.440, -0.440, -0.440, -0.440, -0.440, -0.440, 0.551, 0.551, 0.551, 0.551, 0.551, 0.551, 0.551, 0.078, 0.078, 0.078, 0.078, 0.078, 0.078, 0.078, 1.079, 1.079, 1.079, 1.079, 1.079, 1.079, 0.321, 0.321, 0.321, 0.321, 0.321, 0.321, 0.321, -1.001, -1.001, -1.001, -1.001, -1.001, -1.001, -1.001, -0.222, -0.222, -0.222, -0.222, -0.222, -0.222, -0.222, 1.003, 1.003, 1.003, 1.003, 1.003, 1.003, 1.003, 1.003, -0.558, -0.558, -0.558, -0.558, -0.558, -0.558, -0.558, -0.558, -0.558, -0.463, -0.463, -0.463, -0.463, -0.463, -0.463, -0.212, -0.212, -0.212, -0.212, -0.212, -0.212, -0.212, 0.079, 0.079, 0.079, 0.079, 0.079, 0.079, 0.079, 0.344, 0.344, 0.344, 0.344, 0.344, 0.344, 1.271, 1.271, 1.271, 1.271, 1.271, 1.271, 1.271, 1.271, -0.498, -0.498, -0.498, -0.498, -0.498, -0.498, -0.498, 0.487, 0.487, 0.487, 0.487, 0.487, 0.487, 0.487, 0.487, 1.174, 1.174, 1.174, 1.174, 1.174, 1.174, 1.174, 1.174, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, -0.518, -0.518, -0.518, -0.518, -0.518, -0.518, -0.518, -0.390, -0.390, -0.390, -0.390, -0.390, -0.390, -0.390, -0.390, 0.392, 0.392, 0.392, 0.392, 0.392, 0.392, 0.392, 0.392, 0.392, 0.392, 0.392, -0.048, -0.048, -0.048, -0.048, -0.048, -0.048, -0.048, -0.048, -0.048, -0.048, -0.048, -0.048, 0.722, 0.722, 0.722, 0.722, 0.722, 0.722, 0.722, 0.722, 0.722, 0.722, 0.722, 0.722, 0.722, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.129, 0.079, 0.079, 0.079, 0.079, 0.079, 0.079, 0.079, 0.079, 0.079, 0.079, 0.079, 0.079, 0.079]
    refSignal = get_signal_file('../../data/ONTBarcode12/signal_0.txt')
    querySignalAdapter = get_signal_file('../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/top_and_adapter/timeSeries_0.txt')
    querySignalAdapter_kit = [-0.588, -0.588, -0.588, -0.588, -0.588, -0.588, -0.588, -0.588, -0.787, -0.787, -0.787, -0.787, -0.787, -0.787, -0.787, 1.136, 1.136, 1.136, 1.136, 1.136, 1.136, 1.136, -0.578, -0.578, -0.578, -0.578, -0.578, -0.578, -0.578, -0.578, 0.709, 0.709, 0.709, 0.709, 0.709, 0.709, 0.709, 0.709, 0.709, -0.197, -0.197, -0.197, -0.197, -0.197, -0.197, -0.197, -0.197, -0.482, -0.482, -0.482, -0.482, -0.482, -0.482, -0.482, 0.407, 0.407, 0.407, 0.407, 0.407, 0.407, 0.407, 1.156, 1.156, 1.156, 1.156, 1.156, 1.156, 1.156, 1.156, 1.156, 1.156, 0.181, 0.181, 0.181, 0.181, 0.181, 0.181, 0.181, -1.904, -1.904, -1.904, -1.904, -1.904, -1.904, -1.904, -1.904, -1.904, -1.904, 0.345, 0.345, 0.345, 0.345, 0.345, 0.345, 1.197, 1.197, 1.197, 1.197, 1.197, 1.197, 1.197, 1.197, 1.197, 0.466, 0.466, 0.466, 0.466, 0.466, 0.466, 0.466, 0.466, -0.960, -0.960, -0.960, -0.960, -0.960, -0.960, -0.960, -0.960, -1.879, -1.879, -1.879, -1.879, -1.879, -1.879, -1.879, -1.879, 0.641, 0.641, 0.641, 0.641, 0.641, 0.641, 0.641, 0.641, 0.562, 0.562, 0.562, 0.562, 0.562, 0.562, 0.562, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.348, 0.348, 0.348, 0.348, 0.348, 0.348, 0.348, 0.348, 0.348, -0.860, -0.860, -0.860, -0.860, -0.860, -0.860, -0.860, -0.860, -0.860, -0.860, 0.881, 0.881, 0.881, 0.881, 0.881, 0.881, 0.881, 0.881, 0.881, 0.881, -0.972, -0.972, -0.972, -0.972, -0.972, -0.972, -0.972, -0.972, -0.972, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783, 0.783, 0.686, 0.686, 0.686, 0.686, 0.686, 0.686, 0.686, 0.686, 0.686, 0.686, -0.467, -0.467, -0.467, -0.467, -0.467, -0.467, -0.467, -0.467, -0.467, -0.467, 0.169, 0.169, 0.169, 0.169, 0.169, 0.169, 0.169, 0.169, 0.169, 0.169, 0.169, 0.169, 0.140, 0.140, 0.140, 0.140, 0.140, 0.140, 0.140, 0.140, 0.140, 0.140, 0.140, 0.140]
    position_start = fromLongRefFindShortQuery(refSignal, querySignalAdapter_kit)[1]
    barcodeSig = refSignal[position_start + 40: position_start + 40 + 310]
    trueBarcodeSigPositon = fromLongRefFindShortQuery(refSignal, trueONT12Barcode1Sig)
    plt.plot(refSignal[0:2000])
    plt.plot([i for i in range(position_start + 40, position_start + 40 + 310)], refSignal[position_start + 40: position_start + 40 + 310])
    plt.axvline(0, color= 'r')
    plt.axvline(position_start + 40, color= 'r')
    plt.axvline(position_start + 40 + 310, color= 'r')
    plt.text(position_start + 40 + 20, -2.7, position_start + 40,color='r')
    plt.text(position_start + 40 + 310 + 20, -2.7, position_start + 40 + 310,color='r')
    # plt.text(20, -2.7, position_start + 40 + 310,color='r')
    plt.subplot(2,2,2)
    plt.plot(refSignal[0:2000])
    plt.axvline(trueBarcodeSigPositon[0], color= 'r')
    plt.axvline(trueBarcodeSigPositon[1], color= 'r')
    plt.text(trueBarcodeSigPositon[0] + 20, -2.7, trueBarcodeSigPositon[0],color='r')
    plt.text(trueBarcodeSigPositon[1] + 20, -2.7, trueBarcodeSigPositon[1],color='r')
    plt.plot([i for i in range(trueBarcodeSigPositon[0], trueBarcodeSigPositon[1])], refSignal[trueBarcodeSigPositon[0]:trueBarcodeSigPositon[1]])
    plt.subplot(2,2,3)
    plt.plot(barcodeSig)
    plt.subplot(2,2,4)
    plt.plot(trueONT12Barcode1Sig)
    plt.savefig('cutBarcodeSig.png', dpi=350)
    plt.show()
    # test = RefromLongRefFindShortQuery_twoEnd(ref_signal_path = '../../exper_15bp_randomVSOurs/data/SelectedBarcodeLongFlanking100_15bp/barcode_0/signal/signal_6.txt', 
    #                                     bottom_signal_path = '../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/bottom/timeSeries_0.txt', 
    #                                     topAndAdapter_signal_path = '../../exper_15bp_randomVSOurs/data/longFlankingSig/noise/top_and_adapter/timeSeries_0.txt')
    # print(test, len(test[2]))

    # id1, id2 = 0, 0
    # # # checkFindFlankingSeqPosition(sigId = id, mode = 'top')
    # checkFindFlankingSeqPosition(sigId = id1, sigId2 = id2, mode = 'adapter')
    # checkFindFlankingSeqPosition(sigId = id1, sigId2 = id2, mode = 'bottom')
    

    
    

    # test1 = fromLongRefFindShortQuery_test(ref_signal_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignal/signal_1.txt'
    # , query_signal_path = 'testData/top_adapter_signal.txt')
    
    # print(test1)
    
    # # test2 = fromLongRefFindShortQuery_test(ref_signal_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignal/signal_1.txt'
    # # , query_signal_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/BarcodeSignal15bp831/timeSeries_0.txt')
    
    # test2 = fromLongRefFindShortQuery_test(ref_signal_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignal/signal_1.txt'
    # , query_signal_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/OursBarcodeSignal15bp831_deepsim/timeSeries_0.txt')

    # print(test2)
    

    # # test2 = fromLongRefFindShortQuery_test(ref_signal_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignal/signal_2.txt'
    # # , query_signal_path = 'testData/top_adapter_signal.txt')
    # sig1 = get_signal_file(filetxt_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignal/signal_55.txt')[test1[1] + 40:test1[1] + 40 + 100 + 150]
    # sig = get_signal_file(filetxt_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/OursBarcodeSignal15bp831_deepsim/timeSeries_0.txt')
    # sig2 = get_signal_file(filetxt_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignal/signal_1.txt')[test2[0]:test2[1]]

    # # # print(test1, test2)

    # # # sig = get_signal_file(filetxt_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignalParsedBarcodeSig/timeSeries_0.txt')
    # # # sig = get_signal_file(filetxt_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignal/signal_0.txt')[500:800]
    # plt.plot(sig)
    # plt.plot(sig1)
    # # # sig2 = get_signal_file(filetxt_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignal/signal_2.txt')[test2[0]:test2[1]]
    # # # sig2 = get_signal_file(filetxt_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignal/signal_10.txt')[500:800]
    # # # sig2 = get_signal_file(filetxt_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/UMI831AllSignalParsedBarcodeSig/timeSeries_6.txt')
    # plt.plot(sig2)
    # plt.show()
    # drawFigureForTest(position1 = test[0], position2 = test[1])

    # drawFigureForCheckbarcodeSig(position1 = test[1])


