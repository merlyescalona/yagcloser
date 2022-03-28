#!/usr/bin/python
import pysam,gzip, argparse,datetime, collections, csv, os, sys, filetype, logging, copy, re
import numpy as np
from Bio import SeqIO, Seq, AlignIO
# ========================================================================================
PROGRAM_NAME="update_assembly_edits_and_breaks"
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
                APPLOGGER.debug(s)
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
        regexY=re.compile(r"Y")
        regexZ=re.compile(r"Z")
        APPLOGGER.info("Starting...")
        scaffolds=readGenomeAssembly(parameters)
        edits=readEditsFile(parameters,scaffolds)
        newscaffolds=[]
        newoutput=parameters['output_filename_full']
        # Generating scaffolds/contigs  with no modifications
        APPLOGGER.info("Getting scaffolds/contigs with no modifications...")
        k=0
        for sc in sorted(scaffolds):
            if sc not in edits.keys():
                k+=1
                # k=int(sc.split("_")[1].split(";")[0])
                newscaffolds+=[(scaffolds[sc], k)]
        # SOI recording
        soifull=[]
        # Generating new assembly - applying changes
        APPLOGGER.info("Applying changes...")
        for sc in edits:
            sequence=str(scaffolds[sc].seq)
            newSeq=None
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
                    APPLOGGER.info("Filling gap...")
                    tofill="Y"+coordinates['sequence']+"Y"
                    if(coordinates['sequence']=="-"):
                        tofill=coordinates['sequence']
                    seq=sequence[0:start]+tofill+sequence[end:]
                    sequence=seq
                else:
                    APPLOGGER.info("Closing gap...")
                    seq=sequence[0:start]+"Z"+"Z"+sequence[end:]
                    sequence=seq
            # Breaking 
            seqs=[s for s in sequence.split("-") if s != ""]
            if (len(seqs)>1):
                APPLOGGER.info("Breaking...[{}]{:5}".format(sc,len(seqs)))
                for cIndex in range(0,len(seqs)):
                    k+=1 
                    APPLOGGER.info("Breaking ...{:5}/{:5}".format(cIndex,len(seqs)))
                    counterClosing=len(regexY.findall(seq[cIndex]))/2
                    counterFilling=len(regexZ.findall(seq[cIndex]))/2
                    newSeq2=SeqIO.SeqRecord(\
                        seq=Seq.Seq(seqs[cIndex].replace("X", "N").replace("Z","").replace("Y", "")),\
                        id="{};YAGPart={} Closings={}, Fillings={}".format(sc,cIndex,counterClosing,counterFilling),\
                    )
                    newscaffolds+=[(newSeq2, k, cIndex)]
                    currentSequence=seqs[cIndex].replace("Y","Z")
                    while not (currentSequence.find("Z") == -1):
                        startz=currentSequence.find("Z")
                        currentSequence=currentSequence[0:startz]+currentSequence[(startz+1):]
                        APPLOGGER.info("{}\t{}".format(currentSequence[currentSequence.find("Z")],currentSequence.find("Z")))               
                        endz=currentSequence.find("Z")
                        currentSequence=currentSequence[0:endz]+currentSequence[(endz+1):]
                        edittype="FILL"
                        if (startz+1 ==endz):
                            edittype="CLOSE"
                        soifull+=[[sc,startz,endz, edittype]]
            else:
                APPLOGGER.info("No breaks...{}".format(sc))
                k+=1 
                APPLOGGER.debug("No breaks...{}: Length of sequence={}".format(sc,len(sequence)))
                counterClosing=len(regexY.findall(sequence))/2
                counterFilling=len(regexZ.findall(sequence))/2
                APPLOGGER.debug("No breaks...{}: Closings={}, Fillings={}".format(sc,counterClosing,counterFilling))
                newSeq=SeqIO.SeqRecord(\
                    seq=Seq.Seq(sequence.replace("X", "N").replace("Z","").replace("Y","")),\
                    id="{}_{}".format(sc, parameters['output_suffix']),\
                    description="Closings={}, Fillings={}".format(counterClosing,counterFilling)\
                )   
                newscaffolds+=[(newSeq, k)]
                APPLOGGER.debug("{}\t{}".format(sequence[sequence.find("Z")],sequence.find("Z")))
                sequence=sequence.replace("Y","Z")
                APPLOGGER.debug("No breaks...{}: Replaced Ys and Zs".format(sc))
                while not (sequence.find("Z") == -1):
                    APPLOGGER.debug("No breaks...{}: Finding Zs".format(sc))
                    startz=sequence.find("Z")
                    sequence=sequence[0:startz]+sequence[(startz+1):]
                    APPLOGGER.debug("No breaks...{}: Sequence chunk".format(sc))
                    APPLOGGER.info("{}\t{}".format(sequence[sequence.find("Z")],sequence.find("Z")))               
                    endz=sequence.find("Z")
                    sequence=sequence[0:endz]+sequence[(endz+1):]
                    edittype="FILL"
                    if (startz+1 ==endz):
                        edittype="CLOSE"
                    soifull+=[[sc,startz,endz, edittype]]
        # Sorting
        APPLOGGER.info("Sorting scaffolds...")
        newscaffolds.sort(key=lambda k:k[1])
        newscaffolds=[sc[0] for sc in newscaffolds]
        # Output
        APPLOGGER.info("Writing new genome assembly...")
        with open(newoutput, "w") as handle: 
            _=SeqIO.write(newscaffolds,handle,"fasta")
        # If contigs....
        APPLOGGER.info("Writing coordinates of interests...")
        with open(parameters['soi_full_filename'], 'a') as handle:
            for item in soifull:
                handle.write("{soi[0]}\t{soi[1]}\t{soi[2]}\n".format(soi=item))
    except Exception as ex:
        return -1, ex.message
    return 0,"Done."


def main():
    parser = argparse.ArgumentParser(\
    prog=PROGRAM_NAME,\
        description='Updates (generates new version) of a genome assembly' +\
            ' file given the scaffold, start change position, end change position'+\
            ' and sequence (if any) in a  BED-like (TSV) file, and the current genome'+\
            ' assembly. Output will include the updated version of the given genome'+\
            ' and a file with the edits relative to the new version of the genome.'+\
            ' Composition table output is under development',\
        add_help=True\
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
        '-e','--edits',\
        metavar = '<file_path>',\
        type = str,\
        required = True,\
        help = 'BED-like (Tab separated) file with the following format:\n<scaffold_name>   <start_change_position>    <end_change_position>    <sequence_if_any>\n'+\
            'If the edits that will be made is a break, <sequence_if_any> should be the symbol: "-".')     
    optionalArts=parser.add_argument_group("{0}Optional arguments{1}".format("\033[1m","\033[0m"))
    optionalArts.add_argument(\
        '-s','--suffix',\
        metavar = '<modified_sequence_suffix>',\
        type = str,\
        default= 'MOD',\
        help = 'Suffix used in the sequence description to identify modified scaffolds. Default: MOD')
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
    try:
        args=parser.parse_args()
        parameters=dict()
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
        # Checkiong output file
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
        filename, file_extension = os.path.splitext(OUTPUT)
        parameters['output_suffix']=args.suffix
        parameters['soi_full_filename']="{}.soi.full.txt".format(filename)
        parameters['soi_contig_filename']="{}.soi.contig.txt".format(filename)
        return run(parameters)
    except:
        message="Parsing error. Exiting"
        return (-1, message)


if __name__ == "__main__":
    SIGNAL_EXIT,MESSAGE=main()
    if SIGNAL_EXIT == -1:
        APPLOGGER.error(MESSAGE)
    else:
        APPLOGGER.info(MESSAGE)
    sys.exit(SIGNAL_EXIT)
