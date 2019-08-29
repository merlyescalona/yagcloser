#!/usr/bin/python
import pysam,gzip, argparse,datetime, collections, csv, os, sys, filetype, logging, copy
import numpy as np
from Bio import SeqIO, Seq, AlignIO
# ========================================================================================
PROGRAM_NAME="update_assembly_edits_contigs"
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
    # read genome assembly file    d=d[which(d$coord1<start),]
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

#=========================================================================================

def readEditsFile(parameters, scaffolds):
    APPLOGGER.info("Reading edits file...")
    edits=dict()
    with open(parameters['edits_filename'], "rt") as handle:
        lines=handle.readlines()
        for line in lines:
            s=line.strip().split()
            if s:
                scaffold=s[0]
                start=float(s[1])
                end=float(s[2])
                sequence=""
                print(s)
                # if len(s)==6:
                if len(s)==4:
                    sequence=s[3].replace("N", "X")
                try:
                    edits[scaffold]+=[{'start':int(start),'end':int(end),'sequence':sequence}]
                except:
                    edits[scaffold]=[{'start':int(start),'end':int(end),'sequence':sequence}]
    return edits

def run(parameters):
    try:
        APPLOGGER.info("Starting...")
        scaffolds=readGenomeAssembly(parameters)
        edits=readEditsFile(parameters,scaffolds)
        generateContigs=parameters['generate_contigs']
        newscaffolds=[]
        contigs=[]
        newoutput=parameters['output_filename_full']
        newoutput2=parameters['output_filename_contigs']
        # Generating scaffolds/contigs  with no modifications
        APPLOGGER.info("Getting scaffolds/contigs with no modifications...")
        for sc in scaffolds:
            if sc not in edits.keys():
                k=int(sc.split("=")[1])
                # k=int(sc.split("_")[1])
                newscaffolds+=[(scaffolds[sc], k)]
                if generateContigs:
                    seqs=[s for s in str(scaffolds[sc].seq).split("N") if s != ""]
                    for cIndex in range(0,len(seqs)):
                        newSeq2=SeqIO.SeqRecord(\
                            seq=Seq.Seq(seqs[cIndex]),\
                            id="{};Part={}".format(sc,cIndex),\
                        )
                        contigs+=[(newSeq2, k, cIndex)]
        # SOI recording
        soicontigs=[]
        soifull=[]
        # Generating new assembly - applying changes
        APPLOGGER.info("Applying changes...")
        for sc in edits:
            sequence=str(scaffolds[sc].seq)
            newSeq=None
            counterClosing=0
            counterFilling=0
            allcoordinates=[dd for dd in edits[sc]]
            allcoordinates.sort(key=lambda k: k['start'], reverse=True)
            for i in allcoordinates:
                APPLOGGER.debug(i)
            for coordinates in allcoordinates:
                APPLOGGER.debug("{}\t{}\t{}\t{}\t{}".format(sc,coordinates['start'],coordinates['end'],len(coordinates), coordinates['sequence']))
                start=int(coordinates['start'])
                end=int(coordinates['end'])
                APPLOGGER.debug(sequence[start-1:end+1])
                if not coordinates['sequence'] == "":
                    # check >ScNM3vo_3;HRSCAF=5:34934647-34934648
                    # if sc=="ScNM3vo_2470;HRSCAF=2499" and start==29753528:
                        # start+=5;end+=5
                    APPLOGGER.info("Filling gap...")
                    seq=sequence[0:start]+"Z"+coordinates['sequence']+"Z"+sequence[end:]
                    sequence=seq
                    counterFilling+=1
                else:
                    APPLOGGER.info("Closing gap...")
                    seq=sequence[0:start]+"Z"+"Z"+sequence[end:]
                    sequence=seq
                    counterClosing+=1
            newSeq=SeqIO.SeqRecord(\
                seq=Seq.Seq(sequence.replace("X", "N").replace("Z","")),\
                id="{}-MOD".format(sc),\
                description="Closings={}, Fillings={}".format(counterClosing,counterFilling)\
            )
            # k=int(sc.split("_")[1])
            k=int(sc.split("=")[1])
            newscaffolds+=[(newSeq, k)]
            sequence2=copy.copy(sequence)
            APPLOGGER.info("{}\t{}".format(sequence[sequence.find("Z")],sequence.find("Z")))
            while not (sequence.find("Z") == -1):
                startz=sequence.find("Z")
                sequence=sequence[0:startz]+sequence[(startz+1):]
                APPLOGGER.info("{}\t{}".format(sequence[sequence.find("Z")],sequence.find("Z")))               
                endz=sequence.find("Z")
                sequence=sequence[0:endz]+sequence[(endz+1):]
                soifull+=[[sc,startz,endz]]
            if generateContigs:
                APPLOGGER.info("Generating contig-only assembly")
                seqs=[s for s in sequence2.split("N") if s != ""]
                for cIndex in range(0,len(seqs)):
                    APPLOGGER.info("{}\t{}".format(seqs[cIndex][seqs[cIndex].find("Z")],seqs[cIndex].find("Z")))
                    while not (seqs[cIndex].find("Z") == -1):
                        startz=seqs[cIndex].find("Z")
                        seqs[cIndex]=seqs[cIndex][0:startz]+seqs[cIndex][(startz+1):]
                        APPLOGGER.info("{}\t{}".format(seqs[cIndex][startz],startz))
                        endz=seqs[cIndex].find("Z")
                        seqs[cIndex]=seqs[cIndex][0:endz]+seqs[cIndex][(endz+1):]
                        soicontigs+=[[sc,cIndex,startz,endz]]
                    newSeq2=SeqIO.SeqRecord(\
                        seq=Seq.Seq(seqs[cIndex].replace("X", "N")),\
                        id="{};Part={}-MOD".format(sc,cIndex),\
                    )
                    contigs+=[(newSeq2, k, cIndex)]
        # Sorting
        APPLOGGER.info("Sorting scaffolds...")
        newscaffolds.sort(key=lambda k:k[1])
        newscaffolds=[sc[0] for sc in newscaffolds]
        # Output
        APPLOGGER.info("Writing new genome assembly...")
        with open(newoutput, "w") as handle: 
            _=SeqIO.write(newscaffolds,handle,"fasta")
        # If contigs....
        if generateContigs:
            APPLOGGER.info("Writing new genome assembly (contigs-only)...")
            contigs.sort(key=lambda k:(k[1], k[2]))
            contigs=[sc[0] for sc in contigs]
            with open(newoutput2, "w") as handle: 
                _=SeqIO.write(contigs,handle,"fasta")
        APPLOGGER.info("Writing coordinates of interests...")
        with open(parameters['soi_contig_filename'], 'a') as handle:
            for item in soicontigs:
                handle.write("{soi[0]}\t{soi[1]}\t{soi[2]}\t{soi[3]}\n".format(soi=item))
        with open(parameters['soi_full_filename'], 'a') as handle:
            for item in soifull:
                handle.write("{soi[0]}\t{soi[1]}\t{soi[2]}\n".format(soi=item))
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
        '-c','--contigs',\
        required = True,\
        action='store_true',\
        default=True,\
        help = 'Generate extra output file of a contig-only assembly. ')        
    requiredArgs.add_argument(\
        '-e','--edits',\
        metavar = '<file_path>',\
        type = str,\
        required = True,\
        help = 'Tab separated file with the following format:\n<scaffold_name>   <gap_id>')     
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
        parameters=dict()
        if args.help: 
            return(0,"Finished")
        # Checking files
        if args.log == 'DEBUG':
            APPLOGGER.setLevel(logging.DEBUG)
        else:
			APPLOGGER.setLevel(logging.INFO)
        APPLOGGER.debug(args)
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
        # Checking edits file
        if (os.path.exists(os.path.abspath(args.edits))):
			EDITS=os.path.abspath(args.edits)
        else:
			message="Edits file does not exist. Please verify. Exiting."
			return (-1, message)
        parameters['genome_filename']=REFFILE
        parameters['edits_filename']=EDITS
        parameters['output_filename_full']=OUTPUT
        parameters['generate_contigs']=args.contigs
        filename, file_extension = os.path.splitext(OUTPUT)
        parameters['output_filename_contigs']="{}.contigs{}".format(filename, file_extension)
        parameters['soi_full_filename']="{}.soi.full.txt".format(filename)
        parameters['soi_contig_filename']="{}.soi.contig.txt".format(filename)
        return run(parameters)
    except:
        message="Parsing error. Exiting"
        # parser.print_usage()
        return (-1, message)


if __name__ == "__main__":
    SIGNAL_EXIT,MESSAGE=main()
    if SIGNAL_EXIT == -1:
        APPLOGGER.error(MESSAGE)
    else:
        APPLOGGER.info(MESSAGE)
    sys.exit(SIGNAL_EXIT)

# APPLOGGER.info("Starting...")
# scaffolds=readGenomeAssembly(parameters)
# edits=readEditsFile(parameters,scaffolds)
# generateContigs=parameters['generate_contigs']
# newscaffolds=[]
# contigs=[]
# newoutput=parameters['output_full']
# newoutput2=parameters['output_contigs']


# assemblyfile="/scratch1/merly/Mende/HG03486.assembly.3.fasta"
# editsfile="/scratch1/merly/Mende/gapfilling/mende.3.full.edits.txt"
# newoutput="/scratch1/merly/Mende/HG03486.assembly.4.fasta"
# newoutput2="/scratch1/merly/Mende/HG03486.assembly.5.fasta"



