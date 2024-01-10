
import os
from subprocess import call, PIPE, STDOUT
import argparse

def clusterRes2DemultiplexRes(clusterFile, barcodeSigDirPath, trueBarcodeSigDirPath, outFile, sigRootName, mode='voting'):
    """get the demultiplexing results based on clustering result."""
    barcodeNum = len(os.listdir(trueBarcodeSigDirPath))
    if mode == "voting":  # randomly select signals, voting to get results.
        command = "./bin/mainCluster2DemRes %s %s %s %s %s %d" \
                %(clusterFile, barcodeSigDirPath, sigRootName, trueBarcodeSigDirPath, outFile, barcodeNum)

        runcode = call(command,
        stdout=PIPE,
        stderr=STDOUT,
        shell=True
        )
        if runcode != 0:
            raise Exception('Refining clusters failed! Please check the input file!')
        demultiplexResList = []
        with open(outFile) as f:
            lines = f.readlines()
            for line in lines[1:]:
                info = line.strip('\n').split(' ')
                demIdx = int(info[-1])
                demultiplexResList.append(demIdx)
        return demultiplexResList

def initializationParameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('--barSigDir', type=str, required=True,
                        help='Indicates the path to a folder containing only extracted barcode signal files (format: txt).')
    parser.add_argument('--sigRootName', type=str, required=True,
                        help='Indicates the prefix name of the signal file.')
    parser.add_argument('--sbarSigDir', type=str, required=True,
                        help='Indicates the folder where standard signals are stored.')
    parser.add_argument('--clusterFile', type=str, required=True,
                        help='Indicates the file that storing the clustering result.')
    parser.add_argument('--oDemFile', type=str, required=True,
                        help='Indicates the file for storing final demultiplexing result.')
    
    args = parser.parse_args()
    return args

def main():
    args = initializationParameters()
    clusterRes2DemultiplexRes(clusterFile = args.clusterFile, \
                              barcodeSigDirPath = args.barSigDir, \
                              trueBarcodeSigDirPath = args.sbarSigDir, \
                              outFile = args.oDemFile, sigRootName = args.sigRootName, \
                              mode='voting')

if __name__ == "__main__":
    main()
