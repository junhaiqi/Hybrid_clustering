import os
import sys
from time import time
import argparse
if 'module' not in sys.path:
    sys.path.append('module')

from mainHybridClusteringClassRefine import DNASeqAndSig

def initializationParameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('--barSigDir', type=str, required=True,
                        help='Indicates a path to a folder containing only extracted barcode signal files (format: txt).')
    parser.add_argument('--barSeqFile', type=str, required=True,
                        help='Indicates a fasta file that containing extracted barcode sequences.')
    parser.add_argument('--sigRootName', type=str, required=True,
                        help='Indicates the prefix name of the signal file.')
    parser.add_argument('--oclusterFile', type=str, required=True,
                        help='Indicates a file for storing final clustering result.')
    parser.add_argument('--precise', type=int, required=False, default=1,
                        help='When the estimated number of clusters is greater than 100, it is recommended to set it to 1, otherwise, set it to 0.')
    args = parser.parse_args()
    return args

def main():
    args = initializationParameters()
    clusteringStructure =  DNASeqAndSig(DNAFilePath = args.barSeqFile, \
                            SigDirPath = args.barSigDir, \
                            finalClusteringFile = args.oclusterFile, \
                            sigRootName = args.sigRootName)
    
    if args.precise == 1:
        # thresFactorBG = 4.3
        # thresFactorRefine = 4
        # thresFactorGM = 4
        thresFactorBG = 4.7
        thresFactorRefine = 4
        thresFactorGM = 4
        
    else:
        thresFactorBG = 3
        thresFactorRefine = 3
        thresFactorGM = 3

    clusteringRes = clusteringStructure.mainHybridClustering(thresFactorBG = thresFactorBG, \
                                                             thresFactorRefine = thresFactorRefine, \
                                                             thresFactorGM = thresFactorGM)  # defalut: 4.3, 3.3, 3.3
                                                             # When the number of barcode is small(<100), 
                                                             # 'thresFactorGM', 'thresFactorRefine' and 'thresFactorBG' can be set 3.
                                                    

if __name__ == "__main__":
    main()
    
    
