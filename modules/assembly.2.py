#!/usr/bin/python
import pysam,gzip, argparse,datetime, collections, csv, os, sys, filetype, logging
import numpy as np
from Bio import SeqIO, Seq

# ========================================================================================
PROGRAM_NAME="update_assembly_with_breaks"
VERSION=1
MIN_VERSION=0
FIX_VERSION=0
# ========================================================================================
APPLOGGER = logging.getLogger(PROGRAM_NAME)
APPLOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
loggerFormatter = logging.Formatter(\
  fmt = "[%(asctime)s]  %(levelname)s (%(funcName)s|%(lineno)d): %(message)s",\
  datefmt = "%d/%m/%Y %I:%M:%S %p")
ch.setFormatter(loggerFormatter)
APPLOGGER.addHandler(ch)
#=========================================================================================
def readGenomeAssembly(parameters):
    # read genome assembly file
    scaffolds=dict()
    APPLOGGER.info("Reading genome file...")
    handle=None
    filekind=filetype.guess(parameters['genome_filename'])
    if filekind and filekind.extension in ['gz','GZ','gZ','Gz']:
        APPLOGGER.info("Running gzip...")
        handle=gzip.open(parameters['genome_filename'], "rt")
    else:
        handle=open(parameters['genome_filename'], "rt")
    scaffolds=SeqIO.to_dict(SeqIO.parse(handle, "fasta")) 
    handle.close()
    return scaffolds

def readBedFile(parameters):
    #read bedfile
    APPLOGGER.info("Reading BED file...")
    handle=None
    filekind=filetype.guess(parameters['bed_filename'])
    if filekind and filekind.extension in ['gz','GZ','gZ','Gz']:
        APPLOGGER.info("Running gzip...")
        handle=gzip.open(parameters['bed_filename'], "rt")
    else:
        handle=open(parameters['bed_filename'], "rt")
    bed=dict()
    bedCounter=1
    for line in handle:
        bedline=line.strip().split()
        try:
            gapId=bed[bedline[0]][-1][1]+1
            bed[bedline[0]]+=[(bedCounter,gapId, int(bedline[1]), int(bedline[2]), int(bedline[2])-int(bedline[1]))]     
        except:
            bed[bedline[0]]=[(bedCounter,0, int(bedline[1]), int(bedline[2]), int(bedline[2])-int(bedline[1]))]     
        bedCounter+=1
    handle.close()
    return bed

def readEditsFile(parameters, scaffolds):
    bed=readBedFile(parameters)
    edits=dict()
    with open(parameters['edits_filename'], "rt") as handle:
        lines=handle.readlines()
        for lineIndex in range(0,len(lines)):
            line=lines[lineIndex]
            s=line.strip().split()
            APPLOGGER.info("[{:>4}/{:>4}]\t {}".format(lineIndex,len(lines),s))
            if s:
                scaffold=s[0]
                gapid=int(s[1])
                start=bed[scaffold][gapid][2]
                end=bed[scaffold][gapid][3]
                sequence=scaffolds[scaffold][start:end]
                try:
                    edits[scaffold]+=[{'gapid':gapid,'start':int(start),'end':int(end),'sequence':sequence}]
                except:
                    edits[scaffold]=[{'gapid':gapid,'start':int(start),'end':int(end),'sequence':sequence}]
    # adding scaffolds with no edits
    return edits

def run(parameters):
    APPLOGGER.info("Start")
    try:
        scaffolds=readGenomeAssembly(parameters)
        edits=readEditsFile(parameters,scaffolds)
        newscaffolds=[]
        seqkeys=dict(); counter=0
        APPLOGGER.info("Setting scaffold keys...")
        for sc in sorted(scaffolds.keys()):
            seqkeys[sc]=counter
            counter+=1
        APPLOGGER.info("Getting scaffolds without breaks")        
        for sc in scaffolds:
            if sc not in edits.keys():
                # k=int(sc.split("=")[1])
                newscaffolds+=[(scaffolds[sc], seqkeys[sc])]
        APPLOGGER.info("Generating new assembly...")                
        # ====================================================================================
        # generating new assembly
        for sc in edits:
            sequence=str(scaffolds[sc].seq)
            allcoordinates=[dd for dd in edits[sc]]
            allcoordinates.sort(key=lambda k: k['start'], reverse=True)
            for i in allcoordinates:
                APPLOGGER.debug(i)
            for coordinates in allcoordinates:
                APPLOGGER.info("{}\t{}\t{}\t{}\t{}".format(\
                    sc,\
                    coordinates['start'],\
                    coordinates['end'],\
                    len(coordinates),\
                    coordinates['sequence']))
                start=int(coordinates['start'])
                end=int(coordinates['end'])
                sequence=sequence[0:start]+"X"+sequence[end+1:]
            # k=int(sc.split("=")[1])
            k=seqkeys[sc]
            seqs=[s for s in sequence.split("X") if s != ""]
            APPLOGGER.info(len(seqs))
            for cIndex in range(0,len(seqs)):
                APPLOGGER.debug("{:5}/{:5}".format(cIndex,len(seqs)))
                newSeq2=SeqIO.SeqRecord(\
                    seq=Seq.Seq(seqs[cIndex]),\
                    id="{};Part={}".format(sc,cIndex),\
                )
                newscaffolds+=[(newSeq2, k, cIndex)]
        # ====================================================================================
        # Update order of the scaffolds based on name info
        APPLOGGER.info("Updating order of the scaffolds based on name info...")                
        newscaffolds.sort(key=lambda k:k[1])
        newscaffolds=[sc[0] for sc in newscaffolds]
        # ====================================================================================
        # Writing output
        APPLOGGER.info("Writing output...({})".format(parameters['output_filename']))                
        with open(parameters['output_filename'], "w") as handle: 
            _=SeqIO.write(newscaffolds,handle,"fasta")
    except Exception as ex:
        return -1, ex.message
    return 0,"Done."

def main():
    parser = argparse.ArgumentParser(\
    prog=PROGRAM_NAME,\
        description='Updates (generates new version) of a genome assembly file given the' +\
            'gaps positions in a BED file, the gap identifiers per scaffold (starting at '+\
            '0, per scaffold) and the current genome assembly.',\
        add_help=False\
    )
    requiredArgs=parser.add_argument_group("{0}Required arguments{1}".format("\033[1m","\033[0m"))
    requiredArgs.add_argument(\
        '-i','--input',\
        metavar = '<sequences.fasta>',\
        type = str,\
        required = True,\
        help = 'Filepath of the reference genome file to be used (accepts .gz.)')
    requiredArgs.add_argument(\
        '-o','--output',\
        metavar = '<file_path>',\
        type = str,\
        required = True,\
        help = 'Output file path.')        
    requiredArgs.add_argument(\
        '-b','--bed',\
        metavar = '<file_path>',\
        type = str,\
        required = True,\
        help = 'BED file with the position of the gaps (accepts .gz).')        
    requiredArgs.add_argument(\
        '-e','--edits',\
        metavar = '<file_path>',\
        type = str,\
        required = True,\
        help = 'Tab separated file with the following format:\n<scaffold_name>   <gap_id>\nGap identifier 0...num. gaps per scaffold.')     
    optionalArts=parser.add_argument_group("{0}Optional arguments{1}".format("\033[1m","\033[0m"))
    optionalArts.add_argument(\
        '-l','--log',\
        metavar = '<log_level>',\
        type = str,\
        choices = ['INFO','DEBUG'],\
        default= 'INFO',\
        help = 'Verbosity levels.')
    informationGroup= parser.add_argument_group("{0}Information arguments{1}".format("\033[1m","\033[0m"))
    informationGroup.add_argument('-v', '--version',\
        action='version',\
        version='{0}: Version {1}.{2}.{3}'.format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
        help="Show program's version number and exit")
    informationGroup.add_argument('-h', '--help',\
        action='store_true',\
        help="Show this help message and exit")         
    try:
        args=parser.parse_args()
        APPLOGGER.debug(args)
        parameters=dict()
        if args.help: 
            return(0,"Finished")
        # Checking files
        if args.log == 'DEBUG':
            APPLOGGER.setLevel(logging.DEBUG)
        else:
			APPLOGGER.setLevel(logging.INFO)
        # Checking input file
        if (os.path.exists(os.path.abspath(args.input))):
			REFFILE=os.path.abspath(args.input)
        else:
			message="Reference file does not exist. Please verify. Exiting."
			return (-1, message)
        # Checking output file
        if not (os.path.exists(os.path.abspath(args.output))):
			OUTPUT=os.path.abspath(args.output)
        else:
			message="Output file exist. File will NOT be overwritten. Please verify. Exiting."
			return (-1, message)
        # Checking bed file
        if (os.path.exists(os.path.abspath(args.bed))):
			BEDFILE=os.path.abspath(args.bed)
        else:
            message="BED file, with the description of the gaps of the given genome, does not exist. Please verify. Exiting."
            return (-1, message)
        # Checking edits file
        if (os.path.exists(os.path.abspath(args.edits))):
			EDITS=os.path.abspath(args.edits)
        else:
			message="Edits file does not exist. Please verify. Exiting."
			return (-1, message)
        
        parameters['genome_filename']=REFFILE
        parameters['edits_filename']=EDITS
        parameters['bed_filename']=BEDFILE
        parameters['output_filename']=OUTPUT
        return run(parameters)
    except Exception as e:
        message="Parsing error. Exiting"
        parser.print_help()
        return (-1, e)


if __name__ == "__main__":
    SIGNAL_EXIT,MESSAGE=main()
    if SIGNAL_EXIT == -1:
        APPLOGGER.error(MESSAGE)
    else:
        APPLOGGER.info(MESSAGE)
    sys.exit(SIGNAL_EXIT)