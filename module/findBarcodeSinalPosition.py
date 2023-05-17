import sys
if 'queryLocalSignal/module' not in sys.path:
    sys.path.append('queryLocalSignal/module')
if 'queryLocalSignal/' not in sys.path:
    sys.path.append('queryLocalSignal/')

import numpy as np

from findLocalSignalPosition import fromLongRefFindShortQuery

def get_signal_file(filetxt_path):

    signal_list = list()
    with open(filetxt_path, 'r') as f:
        signal_file = f.readlines()
        for line in signal_file:
            signal_Value = float(line.rstrip())
            signal_list.append(signal_Value)

    signal = np.array(signal_list)
    return signal

def findBarcodeSinalPosition(signalFile='signal.txt', outAdapterSignalDir='tempoutput/adpterSig', \
                                 outBarcodeSigFile='testBarcode.txt', BarcodeLength = 40):
    """get the barcode signal from the read nanopore signal."""
    estimateBarcodeLength = int(BarcodeLength * 7.75) # a super parameter.
    refSignal = get_signal_file(filetxt_path = signalFile)[0:1800]
    querySignalAdapterPath = outAdapterSignalDir + '/' + 'timeSeries_0.txt'
    
    querySignalAdapter = get_signal_file(filetxt_path = querySignalAdapterPath)
    position_start = fromLongRefFindShortQuery(refSignal, querySignalAdapter)[1]
    barcodeSig = refSignal[position_start + 41: position_start + 41 + estimateBarcodeLength]
    with open(outBarcodeSigFile, 'w') as f:
        for item in barcodeSig:
            f.write('%s\n'%str(item))
    return barcodeSig