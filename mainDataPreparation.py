import sys
import os
if 'module' not in sys.path:
    sys.path.append('module')
if 'module/generate_nanoTruesig/module/' not in sys.path:
    sys.path.append('module/generate_nanoTruesig/module/')
from mainHybridClusteringClass import DNASeqAndSig
from generatNoiselessSignal import sequence_to_true_signal
from multiprocessing import Pool
import argparse

class ClusteringDataPipeline(DNASeqAndSig):
    """A class for preparing data for hybrid clustering."""
    def __init__(self, DNAFilePath, SigDirPath, finalClusteringFile = '', sequencedFilePath = '',
                 adpterSeq = 'GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT',
                 sigRootName = 'timeSeries', 
                 barcodeLen = 40, 
                 threadNum = 8,
                 barcodeNoflankingLen = 24,
                 fastqDirPath = '',
                 fastaFile = '',
                 barcodeFile = '',
                 seqPrefixLen = 200):

        super(ClusteringDataPipeline, self).__init__(DNAFilePath, SigDirPath, finalClusteringFile, sequencedFilePath, 
                                                adpterSeq, sigRootName, barcodeLen, threadNum, barcodeNoflankingLen)

        self.fastqDirPath = fastqDirPath
        self.seqPrefixLen = seqPrefixLen
        self.fastaFile = fastaFile
        self.barcodeFile = barcodeFile

    def newProcessSequencedFile(self, fastqFilePath):
        """get the prefix of sequenced DNA sequences in fastq file. Barcode in this prefix sequence."""
        file = fastqFilePath  # A fastq file that include some sequences.
        DNAList = ['A', 'C', 'T', 'G']  # DNA character.
        seqList = []  # A list to load DNA sequences.
        _append = seqList.append
        with open(file) as f:
            lines = f.readlines()
            for line in lines:
                if line[0] in DNAList:
                    _append(line.strip('\n')[0:self.seqPrefixLen])
        return seqList

    def extractBarcodeSeqFromFastq(self, fastqDirPath = '', fastqList = [], mode = 'path'):
        """Input a floder path that only include fastq files, output a fasat file that include pseudo-barcode sequences. """
        if mode == "path":
            fastqFiles = [fastqDirPath + '/' + item for item in os.listdir(path=fastqDirPath)]
        else:
            fastqFiles = fastqList

        allSeqList = []
        for fastq in fastqFiles:
            seqList = self.newProcessSequencedFile(fastq)
            allSeqList += seqList

        return allSeqList

    def getBarcodePosInAllSeqList(self, outFile=None, fastqList = [], mode = 'fasta'):
        """Get barcodes from a list that include DNA sequences and write them into 'outFile'."""
        if mode == 'fastqPath':
            seqList = self.extractBarcodeSeqFromFastq(self.fastqDirPath)  # Sequence list.
        elif mode == 'list':
            if fastqList == []:
                print("Empty file list! Please check it!")
                return 0
            else:
                seqList = self.extractBarcodeSeqFromFastq(mode = 'list', fastqList = fastqList)
        elif mode == 'fasta':
            seqList = []
            file = open(self.fastaFile)
            for line in file:
                if line[0] != '>':
                    line = line.strip('\n')
                    seqList.append(line)
            file.close()
        else:
            print("There are only three modes A, B and C here! Please enter the correct pattern to extract the barcode sequence!")
            exit(-1)

        _func = self.findBarPosByAdpterSeq
        pool = Pool(self.threadNum)
        barcodeList = list(pool.imap(_func, seqList))
        pool.close()
        pool.join()

        if outFile == None:
            outFile = self.fastqDirPath.split('/')[-1] + '.fasta'

        writeFile = open(outFile, 'w')
        t = 0
        for seq in barcodeList:
            writeFile.write('>%d\n'%t)
            writeFile.write('%s\n'%seq)
            t += 1
        writeFile.close()
        return barcodeList

    def getBarSeqFromBarFasta(self):
        seqAndIdxList = []
        with open(self.barcodeFile) as fr:
            lines = fr.readlines()
            t = 0
            for line in lines:
                if line[0] != '>':
                    line = line.strip('\n')
                    seqAndIdxList.append((line, t))
                    t += 1
        return seqAndIdxList

    def getStandradBarSigs(self, outSigDir):
        if not os.path.exists(outSigDir):
            os.makedirs(outSigDir)
        seqAndIdxList = self.getBarSeqFromBarFasta()
        for i in range(len(seqAndIdxList)):
            sequence_to_true_signal(seqAndIdxList[i], output_folder = outSigDir, sigroot=self.sigRootName)

    def getBarcodeSigFiles(self, adapterSigsDir = 'tempoutput/barcodeONT12AdapterSig',
                            outBarcodeDir = '12BarSigs'):

        self.getAdpterSig(output_folder = adapterSigsDir)
        self.mutiFindBarcodeSig(outAdapterSignalDir=adapterSigsDir, barcodeSigDir = outBarcodeDir, past=0)


def initializationParameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sigDir', type=str, required=True,
                        help='Indicates a path to a folder containing only nanopore signal files (format: txt).')
    parser.add_argument('--seqFile', type=str, required=True,
                        help='Indicates a fasta file that containing sequenced DNA sequences.')
    parser.add_argument('--barSeqFile', type=str, required=True,
                        help='Indicates a fasta file that containing barcode sequences.')
    parser.add_argument('--sigRootName', type=str, required=True,
                        help='Indicates the prefix name of the signal file.')
    parser.add_argument('--adapSeq', type=str, required=True,
                        help='Indicates a adapter sequence.')
    parser.add_argument('--oADir', type=str, required=True,
                        help='Indicates a folder for storing standard nanopore signal corresponding to the adapter sequence.')
    parser.add_argument('--oBDir', type=str, required=True,
                        help='Indicates a folder for storing extracted barcode signals from nanopore signals.')
    parser.add_argument('--oTBDir', type=str, required=True,
                        help='Indicates a folder for storing strandard barcode signals.')
    parser.add_argument('--oBF', type=str, required=True,
                        help='Indicates a fasta file to record the extracted barcode sequences.')
    parser.add_argument('--bl', type=int, required=False, default=40,
                        help='Indicates the sequence length of the barcode (including flanking sequences).')
    parser.add_argument('--spl', type=int, required=False, default=200,
                        help='Indicates the length of the prefix sequence. The barcode sequence is in this prefix sequence.')

    args = parser.parse_args()
    return args

                        
def main():
    args = initializationParameters()
    mainPipeLine = ClusteringDataPipeline(DNAFilePath = '', \
                                SigDirPath = args.sigDir, \
                                sigRootName = args.sigRootName, \
                                fastaFile = args.seqFile, \
                                barcodeLen = args.bl, \
                                barcodeFile = args.barSeqFile, \
                                seqPrefixLen = args.spl, \
                                adpterSeq = args.adapSeq)

    mainPipeLine.getBarcodePosInAllSeqList(outFile=args.oBF)  # Get all extracted barcode sequences.
    mainPipeLine.getBarcodeSigFiles(adapterSigsDir = args.oADir, outBarcodeDir = args.oBDir)  # Get all extracted barcode signals.
    mainPipeLine.getStandradBarSigs(outSigDir = args.oTBDir)  # Get all strandard barcode signals from barcode fasta file.

if __name__ == "__main__":
    main()
