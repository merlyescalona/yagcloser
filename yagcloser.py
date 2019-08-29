#! /usr/bin/python
"""
Generates a list of gaps that are potentially to be filled or closed.
(c) 2018 Merly Escalona <mmescalo@ucsc.edu>
"""

import argparse,collections,copy,csv,datetime,filetype
import glob,gzip,logging,os,pysam,random,string,sys
import numpy as np
from Bio import SeqIO, Seq, AlignIO,Alphabet
from StringIO import StringIO
from Bio.Align.Applications import MuscleCommandline, MafftCommandline
import Bio.Align.AlignInfo

###############################################################################
PROGRAM_NAME="yagcloser"
VERSION=1
MIN_VERSION=0
FIX_VERSION=0
###############################################################################

APPLOGGER = logging.getLogger(PROGRAM_NAME)
APPLOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
loggerFormatter = logging.Formatter(\
  fmt = "[%(asctime)s]  %(levelname)s (%(funcName)s|%(lineno)d): %(message)s",\
  datefmt = "%d/%m/%Y %I:%M:%S %p")
ch.setFormatter(loggerFormatter)
APPLOGGER.addHandler(ch)

###############################################################################

SCAFFOLD_ID=0
CIGARVALUES={
    "0":"M","1":"I","2":"D","3":"N",\
    "4":"S","5":"H","6":"P","7":"=",\
    "8":"X","9":"B" \
}

###############################################################################

def scaffold_name_generator(scaffold_name_length=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase+string.ascii_uppercase
    scaffold_name=''.join([random.choice(letters) for i in range(0,scaffold_name_length)])
    return "S{}".format(scaffold_name)

###############################################################################

def parse_bed_file(bedfile):
    """ Reading bedfile to extract gap info"""
    APPLOGGER.info("{}".format("="*80))
    APPLOGGER.info("Reading BED file ({})".format(bedfile))
    APPLOGGER.info("{}".format("="*80))
    # print("Reading BED file ({})".format(bedfile))
    start=datetime.datetime.now()
    bed=dict()
    # BED Dictionary:
    # BED[SCAFFOLD]: GAPID,START,END,SIZE
    bedCounter=1
    filekind=filetype.guess(bedfile)
    handle=None
    if filekind and filekind.extension in ['gz','GZ', 'gZ','Gz']:
        handle=gzip.open(bedfile, 'rb') 
    else:
        handle=open(bedfile, "rt")
    for line in handle:
        bedline=line.strip().split()
        try:
            bed[bedline[0]]+=[(bedCounter,int(bedline[1]), int(bedline[2]), int(bedline[2])-int(bedline[1]))]     
        except:
            bed[bedline[0]]=[(bedCounter,int(bedline[1]), int(bedline[2]), int(bedline[2])-int(bedline[1]))]     
        bedCounter+=1
    end=datetime.datetime.now()
    handle.close()
    APPLOGGER.info("Done reading BED file > {}".format( end-start))
    return bed

###############################################################################


def parse_reference_file(bed, parameters):
    reffile=parameters['REFFILE']
    """ Get reference file info"""
    APPLOGGER.info("{}".format("="*80))  
    APPLOGGER.info("Reading scaffolds from reference file...")
    APPLOGGER.info("{}".format("="*80))
    # print("Reading scaffolds from reference file...")
    new_scaffold_name=scaffold_name_generator(10)
    start=datetime.datetime.now()
    scaffolds=dict()
    counter=0
    filekind=filetype.guess(reffile)
    handle=None
    if filekind and filekind.extension in ['gz','GZ']:
        handle=gzip.open(reffile, 'rb') 
    else:
        handle=open(reffile, "rt")
    for record in SeqIO.parse(handle, "fasta"):
        if (counter%100)==0:APPLOGGER.debug("Scaffold:\t{}".format(counter))
        newid=record.id.split("|")[0]
        try:
            scaffolds[newid]={'id':counter,\
                'newname':"{}_{}".format(new_scaffold_name,counter),\
                'seq':str(record.seq),\
                'size':len(record.seq),\
                'flanks':[],\
                'numgaps':len(bed[newid])}  
        except:
            scaffolds[newid]={'id':counter,\
                'newname':"{}_{}".format(new_scaffold_name,counter),\
                'seq':str(record.seq),\
                'size':len(record.seq),\
                'flanks':[],\
                'numgaps':0}
        counter+=1
    handle.close()
    flankTable=generate_flank_table(scaffolds,bed,parameters['FLANKSIZE'])
    # Updating scaffolds data structure
    for f in flankTable:
        scaffolds[f[2]]['flanks']+=[f[1]]
    for f in flankTable:
        scaffolds[f[2]]['numgaps']=len(scaffolds[f[2]]['flanks'])/2
    end=datetime.datetime.now()
    scaffoldsOrder=sorted(scaffolds.keys())
    for index in range(0,len(scaffoldsOrder)):
        s=scaffoldsOrder[index]
        scaffolds[s]['id']=index
    APPLOGGER.info("Done reading reference file ({}) > {}".format(reffile, end-start))
    return scaffolds

###############################################################################


def generate_flank_table(scaffolds, bed, READFLANK):
    APPLOGGER.info("{}".format("="*80))
    APPLOGGER.info("Extracting flank regions...")
    APPLOGGER.info("{}".format("="*80))
    start=datetime.datetime.now()
    flankTable=[]
    flankCounter=1
    for g in bed:
        for gap in bed[g]:
            gapID=gap[0]
            scaffold=g
            gapStart=gap[1]
            gapEnd=gap[2]
            FSGS=gapStart-READFLANK
            FEGS=gapStart
            FSGE=gapEnd
            FEGE=gapEnd+READFLANK
            if FSGS<0: FSGS=1
            if FEGE>scaffolds[scaffold]['size']: FEGE=scaffolds[scaffold]['size']-1
            flankTable+=[\
            (\
                gapID,\
                flankCounter,\
                scaffold,\
                FSGS,\
                FEGS\
            ),(\
                gapID,\
                flankCounter+1,\
                scaffold,\
                FSGE,\
                FEGE\
            )]
            flankCounter+=2
    end=datetime.datetime.now()
    APPLOGGER.info("Done with flanks extraction: > {}".format(end-start))
    return flankTable

###############################################################################


def compute_intergap_distance(scaffolds, bed,parameters):
    OUTPUT=parameters['OUTPUT']
    SAMPLENAME=parameters['SAMPLENAME']
    start=datetime.datetime.now()
    intergapdistancefile='{}/{}gap.interdistance.txt'.format(OUTPUT, SAMPLENAME)
    distances=[]
    with open(intergapdistancefile, 'a') as f: 
        for sindex in range(0,len(scaffolds)):
            s=[s for s in scaffolds if scaffolds[s]['id']==sindex][0]
            if scaffolds[s]['numgaps'] >0:
                gaps=bed[s]
                for g in range(1,len(gaps)):
                    APPLOGGER.debug(gaps[g])
                    value=gaps[g][1]-gaps[g-1][2]
                    f.write("{:25}\t{:5}\t{:5}\t{}\n".format(s, g, g-1,value))
                    distances+=[value]
    end=datetime.datetime.now()
    APPLOGGER.debug("Intergap distances... ({}) > {}".format(intergapdistancefile, end-start))

###############################################################################


def identify_potential_gaps(bed,scaffolds,parameters):
    BAMFILE=parameters['BAMFILE']
    READFLANK=parameters['FLANKSIZE']
    OUTPUT=parameters['OUTPUT']
    SAMPLENAME=parameters['SAMPLENAME']
    MAPQ_THRESHOLD=parameters['MAPQ_THRESHOLD']
    MIN_SUPPORT_THRESHOLD=parameters['MIN_SUPPORT_THRESHOLD']
    APPLOGGER.info("{}".format("="*80))
    APPLOGGER.info("{}".format("Identifying potential gaps"))
    APPLOGGER.info("{}".format("="*80))
    APPLOGGER.debug("Data structure initialization... ")
    APPLOGGER.debug("{}".format("-"*80))
    # Data structure initialization
    gapCounter=0
    start=datetime.datetime.now()
    numGaps=len([j for i in bed.keys() for j in range(0,len(bed[i]))])
    data=dict()
    for index in range(0,len(bed.keys())):
        scaffold=sorted(bed.keys())[index]
        swgtable=bed[scaffold]
        for gapIndex in range(0,len(swgtable)):
            gapCounter+=1
            gapStart=bed[scaffold][gapIndex][1]
            gapEnd=bed[scaffold][gapIndex][2]
            gapSize=bed[scaffold][gapIndex][3]
            region="{}:{}-{}".format(\
                scaffold,\
                gapStart,\
                gapEnd)
            if gapCounter % 100 == 0:
                APPLOGGER.debug("Gap {:4}/{:4}:  || {:>10,} | {:>10,} ||".format(\
                    gapCounter,numGaps,gapStart,gapEnd))
            data[region]={\
                'scaffold':"",\
                'gapstart':0,\
                'gapend':0,\
                'gapsize':0,\
                'alncov':0,\
                'mapq':dict(),\
                'AS':dict(),\
                'supp':dict(),\
                'primary':dict(),\
                'positions':dict(),\
                'refpositions':dict(),\
                'querylen':dict(),\
                'side':dict()\
            }
    end=datetime.datetime.now()
    APPLOGGER.debug("Done with initialization... > {}".format( end-start))
    ######################################################################################
    APPLOGGER.debug("{}".format("-"*80))
    start=datetime.datetime.now()
    samfile = pysam.AlignmentFile(BAMFILE, "rb")
    numGaps=len([j for i in bed.keys() for j in range(0,len(bed[i]))])
    gapCounter=0
    for scaffold in bed.keys():
        swgtable=bed[scaffold]
        for gapIndex in range(0,len(swgtable)):
            gapCounter+=1
            gapStart=bed[scaffold][gapIndex][1]
            gapEnd=bed[scaffold][gapIndex][2]
            gapSize=bed[scaffold][gapIndex][3]
            region="{}:{}-{}".format(\
                scaffold,\
                max(1,gapStart-READFLANK),\
                min(gapEnd+READFLANK, scaffolds[scaffold]['size']))
            regionKey="{}:{}-{}".format(scaffold,gapStart,gapEnd)
            if gapCounter % 100 == 0:
                APPLOGGER.debug("Gap {:4}/{} ({readFlank:}): [{:<25}] ||{:>12,}|{:>12,}||".format(\
                    gapCounter,\
                    numGaps,\
                    scaffold,\
                    gapStart,\
                    gapEnd,\
                    readFlank=READFLANK\
                ))
            # Checking only PRIMARY alignments, with MAPQ > 20, no DUPLICATES and no SUPPLEMENTARY
            for read in samfile.fetch(region=region):
                side="N"
                data[regionKey]['scaffold']=scaffold
                data[regionKey]['alncov']+=1
                data[regionKey]['gapstart']=int(gapStart)
                data[regionKey]['gapend']=int(gapEnd)
                data[regionKey]['gapsize']=int(gapSize)
                if  (not read.is_secondary) and \
                    (not read.is_supplementary) and \
                    (read.mapping_quality > MAPQ_THRESHOLD) and \
                    (not read.is_duplicate):            
                    if read.reference_start < gapStart and  gapEnd < read.reference_end: 
                        side="B"
                    if read.reference_start <= gapEnd and read.reference_end <= gapEnd: 
                        side="L"
                    if gapEnd <= read.reference_start: 
                        side="R"
                    try:
                        data[regionKey]['mapq'][read.query_name]+=[read.mapping_quality]
                        data[regionKey]['side'][read.query_name]+=[side]
                        data[regionKey]['supp'][read.query_name]+=[read.is_supplementary]
                        data[regionKey]['primary'][read.query_name]+=[False]
                        data[regionKey]['positions'][read.query_name]+=[(read.query_alignment_start,read.query_alignment_end)]
                        data[regionKey]['refpositions'][read.query_name]+=[(read.reference_start,read.reference_end)]
                        data[regionKey]['querylen'][read.query_name]+=[read.query_length]
                        data[regionKey]['AS'][read.query_name]+=[read.get_tag("AS")]
                    except:
                        data[regionKey]['mapq'][read.query_name]=[read.mapping_quality]
                        data[regionKey]['side'][read.query_name]=[side]
                        data[regionKey]['supp'][read.query_name]=[read.is_supplementary]
                        data[regionKey]['primary'][read.query_name]=[False]
                        data[regionKey]['positions'][read.query_name]=[(read.query_alignment_start,read.query_alignment_end)]
                        data[regionKey]['refpositions'][read.query_name]=[(read.reference_start,read.reference_end)]
                        data[regionKey]['querylen'][read.query_name]=[read.query_length]
                        data[regionKey]['AS'][read.query_name]=[read.get_tag("AS")]
    samfile.close()
    end=datetime.datetime.now()
    APPLOGGER.info("Done with the info extraction... > {}".format(end-start))
    ###############################################################################
    # Removing 0-support data
    totalGaps=len(data.keys())
    APPLOGGER.info("{}".format("="*80))
    APPLOGGER.info("{}".format("Removing gaps with no support from analysis... "))
    APPLOGGER.info("{}".format("="*80))
    start=datetime.datetime.now()
    removedGaps=[]
    removedCounter=0
    gapCounter=0
    for gapKey in data.keys():
        if gapCounter % 100 == 0:   
            APPLOGGER.debug("Gap {:4}/{}: {:35}".format(gapCounter,totalGaps,gapKey)) 
        if len(data[gapKey]['mapq']) ==  0: 
            removedGaps+=[(gapKey,data[gapKey]['gapsize'] )]
            removedCounter+=1    
            del data[gapKey]
        gapCounter+=1
    end=datetime.datetime.now()
    APPLOGGER.info("Removed {} gaps.. > {}".format(removedCounter, end-start))
    APPLOGGER.info("-"*80)
    APPLOGGER.info("Reporting removed gaps to file: {}/{}.no.support.gaps.txt".format(OUTPUT, SAMPLENAME))
    with open('{}/{}.no.support.gaps.txt'.format(OUTPUT, SAMPLENAME), "w") as handle:
        for pair in removedGaps:
            gap=pair[0]; size=pair[1]
            handle.write("{}\t{}\n".format(gap, size))
    ###############################################################################
    gapCounter=0
    totalGaps=len(data.keys())
    gapAlnCoverage=dict()
    gapReadCoverage=dict()
    gapAlnCoverageB=dict()
    counterB=dict()
    counterN=dict()
    counterL=dict()
    counterR=dict()
    for k in data.keys():
        gapAlnCoverage[k]=0; gapReadCoverage[k]=0;gapAlnCoverageB[k]=0
        counterB[k]=0; counterR[k]=0; counterL[k]=0; counterN[k]=0
    for gapKey in data.keys():
        if gapCounter % 100 == 0:   
            APPLOGGER.debug("Gap {:4}/{}: {:35}".format(gapCounter,totalGaps,gapKey)) 
            # print("Gap {:4}/{}: {:35}".format(gapCounter,totalGaps,gapKey)) 
        gapCounter+=1
        for readKey in data[gapKey]['side']: 
            indexL=list(np.where(np.array([item for item in data[gapKey]['side'][readKey]])=="L")[0])
            indexR=list(np.where(np.array([item for item in data[gapKey]['side'][readKey]])=="R")[0])
            suppL = [ data[gapKey]['supp'][readKey][i] for i in range(0, len(data[gapKey]['supp'][readKey])) if i in indexL ]
            suppR = [ data[gapKey]['supp'][readKey][i] for i in range(0, len(data[gapKey]['supp'][readKey])) if i in indexR ]
            suppL.sort();         suppR.sort()
            tmp=suppL and suppR
            countCov=len(tmp)-sum(np.array(tmp))
            if countCov> 0: gapReadCoverage[gapKey]+=1
            gapAlnCoverage[gapKey]+=countCov
            valB=sum(np.array([item for item in data[gapKey]['side'][readKey]])=="B")
            valL=sum(np.array([item for item in data[gapKey]['side'][readKey]])=="L")
            valR=sum(np.array([item for item in data[gapKey]['side'][readKey]])=="R")
            valN=sum(np.array([item for item in data[gapKey]['side'][readKey]])=="N")
            counterB[gapKey]+=valB
            counterL[gapKey]+=valL
            counterR[gapKey]+=valR
            counterN[gapKey]+=valN
            gapAlnCoverageB[gapKey]+=min(1, min(valR,valL))
    #=========================================================================================
    toBeKeptBothSides={i:counterB[i] for i in counterB if counterB[i] > MIN_SUPPORT_THRESHOLD }
    #=========================================================================================
    APPLOGGER.info("Reporting potential fillable gaps ({}/{}.potential.fillable.gaps.txt)...".format(OUTPUT, SAMPLENAME))
    with open('{}/{}.potential.fillable.gaps.txt'.format(OUTPUT, SAMPLENAME), "w") as handle:
        handle.write("{:25}\t{:15}\t{:15}\t{}\t{}\t{}\t{}\t{}\n".format(\
                "scaffold","start", "end","gapsize",\
                "supportB","supportL","supportR","supportN"))
        for gapKey in toBeKeptBothSides:
            handle.write("{:25}\t{:12}\t{:12}\t{}\t{}\t{}\t{}\t{}\n".format(\
                data[gapKey]['scaffold'],\
                data[gapKey]['gapstart'],\
                data[gapKey]['gapend'],\
                data[gapKey]['gapsize'],\
                counterB[gapKey],\
                counterL[gapKey],\
                counterR[gapKey],\
                counterN[gapKey]))
    data={gapKey:data[gapKey] for gapKey in toBeKeptBothSides if not gapKey == None}
    return data

###############################################################################

def extract_support_data(data, parameters, reference):
    BAMFILE=parameters['BAMFILE']
    OUTPUT=parameters['OUTPUT']
    SAMPLENAME=parameters['SAMPLENAME']
    READFLANK=parameters['FLANKSIZE']
    starttime=datetime.datetime.now()
    try:
        allreads=dict()
        APPLOGGER.info("{}".format("="*80))
        APPLOGGER.info("{}".format("Extracting support data"))
        APPLOGGER.info("{}".format("="*80))
        samfile = pysam.AlignmentFile(BAMFILE, "rb")
        gapCounter=0
        for gapKey in data.keys():
            gapCounter+=1
            scaffold=gapKey.split(":")[0]
            start=int(gapKey.split(":")[1].split("-")[0])-READFLANK
            end=int(gapKey.split(":")[1].split("-")[1])+READFLANK
            allreads[gapKey]={}
            readlist=[]
            for x in samfile.fetch(contig=scaffold,start=start,end=end): 
                readlist+=[x.query_name]
            readlistfiltered=[r for r in  list(set(readlist)) if r in data[gapKey]['side']]
            allreads[gapKey]={rn:{\
                    "positions":[],\
                    "positionslflank":[],\
                    "positionsrflank":[],\
                    "sequence":"",\
                    "sequence_and_flanks":"",\
                    "reference":"",\
                    "read":"",\
                    "lflank":"",\
                    "rflank":""} for rn in readlistfiltered}
            # region="{}:{}-{}".format(scaffold,start,end)
            APPLOGGER.debug("[{:>10,}/{:>10,}]\t{:50} Spanning/all\t({:>6,}/{:6,})".format(gapCounter, len(data),gapKey,len(readlistfiltered),len(readlist)))
            for read in samfile.fetch(contig=scaffold,start=start,end=end): # Iterator
                if  (not read.is_secondary) and \
                    (not read.is_supplementary) and \
                    (read.mapping_quality > parameters['MAPQ_THRESHOLD']) and \
                    (not read.is_duplicate) and \
                    (read.query_name in readlistfiltered) and \
                    (len(data[gapKey]['side'][read.query_name])>0):
                    # print(data[gapKey]['side'][read.query_name])
                    side=data[gapKey]['side'][read.query_name].pop(0)
                    # refpositions=data[gapKey]['refpositions'][read.query_name].pop(0)
                    # APPLOGGER.debug("[{:>10,}/{:>10,}]\t{}\t{}/{}({})".format(gapCounter, len(data),gapKey,read.query_name, len(data[gapKey]['side'][read.query_name]), side))
                    # print(refpositions,(read.reference_start,read.reference_end))
                    if (side == "B"): #and  (refpositions == (read.reference_start,read.reference_end)):
                        # APPLOGGER.debug("[{:>10,}/{:>10,}]\t{}\t{} BSide".format(gapCounter, len(data),gapKey,read.query_name))
                        cigartupleseqs=[ CIGARVALUES[str(t[0])]*t[1] for t in read.cigartuples]
                        cigarsequence="".join(cigartupleseqs)
                        positionRead=0
                        # left most mapping position # column 4 (1-based in bam|0-based in pysam)
                        positionReference=read.reference_start 
                        readlist=[]
                        # APPLOGGER.debug("[{:>10,}/{:>10,}]\t{}\t{} CIGAR len:{}".format(gapCounter,len(data),gapKey,read.query_name,len(cigarsequence)))
                        for i in range(0, len(cigarsequence)):
                            CIGAR=cigarsequence[i]
                            if positionReference >= start and positionReference <= end:
                                allreads[gapKey][read.query_name]['positions']+=[positionRead]
                            if positionReference >= start and  positionReference < start+READFLANK:
                                allreads[gapKey][read.query_name]['positionslflank']+=[positionRead]
                            if positionReference > end-READFLANK and positionReference <= end:
                                allreads[gapKey][read.query_name]['positionsrflank']+=[positionRead]                        
                            if CIGAR == "M":
                                positionRead+=1 # consumes query
                                positionReference+=1 # consumes reference
                            if CIGAR == "I":
                                positionRead+=1 # consumes query
                            if CIGAR == "D":
                                positionReference+=1 # consumes reference
                            if CIGAR == "N":
                                positionReference+=1 # consumes reference
                            if CIGAR == "S":
                                positionRead+=1 # consumes query
                            if CIGAR == "H":
                                pass
                            if CIGAR == "P":
                                pass
                            if CIGAR == "=":
                                positionRead+=1 # consumes query
                                positionReference+=1 # consumes reference
                            if CIGAR == "X":
                                positionRead+=1 # consumes query
                                positionReference+=1 # consumes reference
                        allreads[gapKey][read.query_name]['positions']=list(set(allreads[gapKey][read.query_name]['positions']))
                        allreads[gapKey][read.query_name]['positionslflank']=list(set(allreads[gapKey][read.query_name]['positionslflank']))
                        allreads[gapKey][read.query_name]['positionsrflank']=list(set(allreads[gapKey][read.query_name]['positionsrflank']))                  
                        minseq=np.min(allreads[gapKey][read.query_name]['positions'])
                        # APPLOGGER.debug("[{:>10,}/{:>10,}]\t{}\t{} MINSEQ".format(gapCounter, len(data),gapKey,read.query_name))
                        maxlflank=np.max(allreads[gapKey][read.query_name]['positionslflank'])
                        # APPLOGGER.debug("[{:>10,}/{:>10,}]\t{}\t{} MAXLFLANK".format(gapCounter, len(data),gapKey,read.query_name))
                        minrflank=np.min(allreads[gapKey][read.query_name]['positionsrflank'])
                        # APPLOGGER.debug("[{:>10,}/{:>10,}]\t{}\t{} MINRFLANK".format(gapCounter, len(data),gapKey,read.query_name))
                        maxseq=np.max(allreads[gapKey][read.query_name]['positions'])
                        # APPLOGGER.debug("[{:>10,}/{:>10,}]\t{}\t{} MAXSEQ".format(gapCounter, len(data),gapKey,read.query_name))
                        allreads[gapKey][read.query_name]['read']=read.query_sequence
                        # APPLOGGER.debug("[{:>10,}/{:>10,}]\t{}\t{} {}".format(gapCounter, len(data)548,gapKey,read.query_name, read.query_sequence))
                        allreads[gapKey][read.query_name]['sequence']=allreads[gapKey][read.query_name]['read'][maxlflank:minrflank+1]
                        allreads[gapKey][read.query_name]['sequence_and_flanks']=allreads[gapKey][read.query_name]['read'][minseq:maxseq]
                        allreads[gapKey][read.query_name]['reference']=reference[scaffold]['seq'][read.reference_start:read.reference_end+1]
                        allreads[gapKey][read.query_name]['lflank']=allreads[gapKey][read.query_name]['read'][minseq:maxlflank]
                        allreads[gapKey][read.query_name]['rflank']=allreads[gapKey][read.query_name]['read'][minrflank:maxseq]
                    else:
                        APPLOGGER.debug("[{:>10,}/{:>10,}]\t{}\t{} Skipping aligned segment...".format(gapCounter, len(data),gapKey,read.query_name))
        samfile.close()
        endtime=datetime.datetime.now()
        APPLOGGER.info("End of extracting support data > {}".format(endtime-starttime))
        APPLOGGER.info("{}".format("="*80))
        APPLOGGER.info("{}".format("Checking for ambiguos decisions..."))
        ambiguosdecisions=[]
        APPLOGGER.debug("="*80)
        APPLOGGER.debug("[{}]\t| {} - {} |\t{}".format("GapKey","MostFrequentSequenceSize","EmptyFlankRatio","Support"))
        for gapKey in allreads.keys():
            flanksizes=[]
            seqsizes=[]
            for readKey in allreads[gapKey]:
                flanksizes+=[len(allreads[gapKey][readKey]['positionslflank']),len(allreads[gapKey][readKey]['positionsrflank'])]
                seqsizes+=[len(allreads[gapKey][readKey]['sequence'])]
            summarySizes=collections.Counter(flanksizes)
            summaryseqsizes=collections.Counter(seqsizes)
            mostFrequenceSequenceSize=[d for d in summaryseqsizes if summaryseqsizes[d] == max(summaryseqsizes.values())][0]
            # Removig - if the emtpy_flanks_ration > emtpy_flanks_ratio_allowed
            if int(sorted(summarySizes.keys())[0]) == 0:
                empty_flanks_ratio=(summarySizes[sorted(summarySizes.keys())[0]]*1.0)/sum([summarySizes[d] for d in summarySizes])
            else:
                empty_flanks_ratio=0
            if empty_flanks_ratio > parameters['EMPTY_FLANKS_THRESHOLD']:
                gap=allreads.pop(gapKey)
                _=data.pop(gapKey)
                ambiguosdecisions+=[(gapKey, len(gap.keys()))]
            APPLOGGER.debug("[{:>42}]\t|{:4,} - {:1.4f}|\t{:4}".format(\
                    gapKey,mostFrequenceSequenceSize,\
                    empty_flanks_ratio,\
                    len(flanksizes)/2))
                    # { d:summarySizes[d] for d in sorted(summarySizes.keys())}))
        APPLOGGER.debug("-"*80)      
        with open("{}/{}.ambiguous.txt".format(OUTPUT,SAMPLENAME), "wb" ) as handle:
            handle.write("gapkey\tsupport\n")
            for index in range(0,len(ambiguosdecisions)):
                APPLOGGER.info("Removing gaps with ambiguous decisions... [{:>40} (Support = {})] ".format(ambiguosdecisions[index][0],ambiguosdecisions[index][1] ))          
                handle.write("{}\t{}\n".format(ambiguosdecisions[index][0],ambiguosdecisions[index][1] ))
        APPLOGGER.info("{}".format("="*80))
        APPLOGGER.info("{}".format("Writing support files"))
        APPLOGGER.info("{}".format("="*80))
        gapCounter=0 
        for gapKey in allreads.keys():
            gapCounter+=1
            scaffold=data[gapKey]['scaffold']
            basefilename="{}-{}-{}".format(reference[scaffold]['newname'],data[gapKey]['gapstart'],data[gapKey]['gapend'])
            APPLOGGER.debug("[{:>10,}/{:>10,}]\tWriting {}...".format(gapCounter,len(data),basefilename))
            gapFilename0="{}/{}.support/{}.fasta".format(OUTPUT,SAMPLENAME, basefilename)
            gapFilename1="{}/{}.fullsupport/{}.fasta".format(OUTPUT,SAMPLENAME, basefilename)
            gapFilename2="{}/{}.flanks/{}.fasta".format(OUTPUT,SAMPLENAME, basefilename)
            gapFilename3="{}/{}.reads/{}.fasta".format(OUTPUT,SAMPLENAME, basefilename)
            MSA=[]
            MSA2=[]
            FLANKS=[]
            PILE=[]
            for readKey in allreads[gapKey]:
                F1=allreads[gapKey][readKey]['lflank']
                F2=allreads[gapKey][readKey]['rflank']
                MSA+=[
                    SeqIO.SeqRecord(
                        seq=Seq.Seq(allreads[gapKey][readKey]['sequence']), \
                        id=readKey\
                    )
                ]
                MSA2+=[
                    SeqIO.SeqRecord(
                        seq=Seq.Seq(allreads[gapKey][readKey]['sequence_and_flanks']), \
                        id="{}_{}_{}".format(readKey, len(F1),len(F2))\
                    )
                ]
                FLANKS+=[
                    SeqIO.SeqRecord(
                        seq=Seq.Seq(F1), \
                        id="{}_F1".format(readKey)\
                    ),\
                    SeqIO.SeqRecord(
                        seq=Seq.Seq(F2), \
                        id="{}_F2".format(readKey)\
                    )
                ]
                PILE+=[
                    SeqIO.SeqRecord(
                        seq=Seq.Seq(allreads[gapKey][readKey]['read']), \
                        id="{}_FULLREAD".format(readKey)\
                    )
                ]
            with open(gapFilename0, "wb") as handle:
                _=SeqIO.write(MSA,handle,"fasta")
            with open(gapFilename1, "wb") as handle:
                _=SeqIO.write(MSA2,handle,"fasta")
            with open(gapFilename2, "wb") as handle:
                _=SeqIO.write(FLANKS,handle,"fasta")
            with open(gapFilename3, "wb") as handle:
                _=SeqIO.write(PILE,handle,"fasta")
    except Exception as e:
        return None, "{}\t{}".format(e.__class__,e.message)
    return data, "Done with file writing..."

###############################################################################

def length_agreement(data,reference, parameters):
    starttime=datetime.datetime.now()
    APPLOGGER.info("Length agreement")
    APPLOGGER.info("="*80)
    message="Done"
    OUTPUT=parameters['OUTPUT']
    SAMPLENAME=parameters['SAMPLENAME']
    ratioThreshold=parameters['RATIO_THRESHOLD']
    fastafiles=glob.glob("{}/{}.fullsupport/*".format(OUTPUT, SAMPLENAME))
    logfile="{}/{}.log".format(OUTPUT,SAMPLENAME)
    errfile="{}/{}.alignment.err".format(OUTPUT, SAMPLENAME)
    with open(errfile,"w") as handle: handle.write("Alignment error output.")
    # Re-populate list with filename, size tuples
    for i in range(0,len(fastafiles)):
        fastafiles[i] = (os.path.basename(fastafiles[i]).split(".")[0], os.path.getsize(fastafiles[i]))
    fastafiles.sort(key=lambda filename: filename[1])
    gapsToBeFilled={}
    for gapKey in data.keys():
        newscaffold=reference[data[gapKey]['scaffold']]['newname']
        gapStart=data[gapKey]['gapstart']
        gapEnd=data[gapKey]['gapend']
        newGapKey="{}:{}-{}".format(newscaffold,gapStart,gapEnd)
        gapsToBeFilled[newGapKey]=data[gapKey]
    seqlengths=dict()
    for gapIndex in range(0,len(fastafiles)):
        gapKey=fastafiles[gapIndex][0]
        illuminafile="{}/{}.fullsupport/{}.fasta".format(OUTPUT,SAMPLENAME, gapKey.replace(":","-"))
        with open(illuminafile, "rt") as handle:
            reads=SeqIO.to_dict(SeqIO.parse(handle, "fasta")) 
        MSA=[reads[record] for record in reads]
        aligner="NONE"
        timing=0
        # Sorting of the sequences based on length (longest first) 
        MSA.sort(key = lambda s: len(s.seq), reverse=True)
        gapalnfilename="{}/{}.pre/{}.fasta".format(OUTPUT,SAMPLENAME, gapKey.replace(":","-"))
        seqlengths[gapKey]=[len(i)for i in MSA]
        #=====================================================================================
        # Sequence selection - based on size agreement
        #=====================================================================================
        if len(MSA) > 1:
            interseqdists=np.zeros((len(MSA),len(MSA)))
            # 1. calculate min distance between sequences to find starting avg. sequence length
            seqlengths[gapKey].sort(reverse=True)
            testSet=copy.deepcopy(seqlengths[gapKey])
            for item in range(0,interseqdists.shape[0]):
                for jtem in range(item,interseqdists.shape[0]):
                    interseqdists[item,jtem]=np.abs(testSet[item]-testSet[jtem])
                    interseqdists[jtem,item]=np.nan
                    if item==jtem: interseqdists[item,jtem]=np.nan
            # 2. select pair of lengths that are closest together
            minDistance=np.nanmin(interseqdists)
            minDistPosition=np.unique([j for i in np.where(interseqdists==minDistance) for j in i])
            # 3. calculate mean length between those closest points
            newlength=np.mean(np.array(testSet)[minDistPosition])
            # 4. get a threshold
            threshold=newlength*0.1
            thresholdRange=[newlength-threshold, newlength+threshold]
            # 5. start calculating if there are sequences within the thresholdRange
            selectedPositions=[j for i in np.where(np.logical_and(np.array(testSet) >= thresholdRange[0],np.array(testSet) <= thresholdRange[1])) for j in i]
            # 6. iteration to get the threshold range
            while(len(selectedPositions) >0):
                selected=[testSet[sp] for sp in selectedPositions]
                testSet=[testSet[t] for t in range(0, len(testSet)) if not t in selectedPositions]
                newlength=np.mean(np.array(selected))
                threshold=newlength*0.1
                thresholdRange=[newlength-threshold, newlength+threshold]
                selectedPositions=[j for i in np.where(np.logical_and(np.array(testSet) >= thresholdRange[0],np.array(testSet) <= thresholdRange[1])) for j in i]
            thresholdRange=[newlength-threshold, newlength+threshold]
            selectedPositions=[j for i in np.where(np.logical_and(np.array(seqlengths[gapKey]) >= thresholdRange[0],np.array(seqlengths[gapKey]) <= thresholdRange[1])) for j in i]
            ratio=len(MSA)*1.0/len(seqlengths[gapKey])
            MSA=[MSA[i] for i in selectedPositions]
            with open(logfile,"a") as handle: 
                handle.write("{}\t{}\t{}\n".format(gapKey, len(MSA)+len(selectedPositions), len(MSA)))
            if len(MSA) > 1 and ratio >= ratioThreshold:
                with open(gapalnfilename, "w") as handle: 
                    _=SeqIO.write(MSA,handle,"fasta")
                #=====================================================================================        
                # Also I can do the alignment
                msafile="{}/{}.msa/{}.fasta".format(OUTPUT,SAMPLENAME, gapKey.replace(":","-"))
                startMUSCLE=datetime.datetime.now()
                aligner_cline = MuscleCommandline(\
                    input=gapalnfilename)
                try:
                    stdout, stderr = aligner_cline()
                    aligner="MUSCLE"
                except:
                    aligner_cline=MafftCommandline(\
                        input=gapalnfilename)
                    stdout, stderr = aligner_cline()
                    aligner="MAFFT"
                with open(msafile,"w") as handle: handle.write(stdout)
                endMUSCLE=datetime.datetime.now()
                timing=endMUSCLE-startMUSCLE
                with open(errfile,"a") as handle: 
                    handle.write("{separator}\n{}\t{}\nInput: {}\nOutput:{}\n{}\n".format(\
                        aligner,timing, stderr, gapalnfilename,msafile,separator="="*80))
                APPLOGGER.debug("[{:40}]\t{}\t({})\t{}".format(gapKey,aligner,len(MSA),timing))
    APPLOGGER.info("Done length agreement process > {}".format(datetime.datetime.now()-starttime))
    APPLOGGER.info("="*80)
    return gapsToBeFilled, message

###############################################################################

def consensus_generation(data,parameters):
    # this data has as gapkeys the new scaffold name
    try:
        starttime=datetime.datetime.now()
        APPLOGGER.info("Consensus generation...")
        APPLOGGER.info("="*80)
        message="Done"
        OUTPUT=parameters['OUTPUT']
        SAMPLENAME=parameters['SAMPLENAME']
        coverageConsensus=parameters['MIN_CONVERAGE_CONSENSUS']
        readFlank=parameters['FLANKSIZE']
        editfile="{}/{}.edits.txt".format(OUTPUT,SAMPLENAME)
        logfile="{}/{}.consensus.log.txt".format(OUTPUT,SAMPLENAME)
        logfile2="{}/{}.gaps.closed.original_coordinates.txt".format(OUTPUT,SAMPLENAME)
        logfile3="{}/{}.no.length.agreement.txt".format(OUTPUT,SAMPLENAME)
        with open(editfile,'w') as log: log.write("")
        with open(logfile2,'w') as log: log.write("")
        with open(logfile,'w') as log: log.write("")
        matches=dict()
        alignedpos=dict()
        index=0
        APPLOGGER.info("Reporting gaps with no length agreement data...")
        with open(logfile3,'w') as log:
            for gapKey in data.keys():
                f="{}/{}.msa/{}.fasta".format(OUTPUT,SAMPLENAME, gapKey.replace(":","-"))
                log.write("{}\n".format(gapKey))
                if not os.path.exists(f): 
                    del data[gapKey]
        totalGaps=len(data.keys())
        for gapKey in data.keys():
            index+=1
            f="{}/{}.msa/{}.fasta".format(OUTPUT,SAMPLENAME, gapKey.replace(":","-"))
            scaffold=data[gapKey]['scaffold']
            start=data[gapKey]['gapstart']
            end=data[gapKey]['gapend']
            align=AlignIO.read(f, "fasta")
            summary=Bio.Align.AlignInfo.SummaryInfo(align)
            cnss=consensus(summary,threshold=0.1,ambiguous='N', require_multiple=coverageConsensus)
            ncount=np.sum(np.array([i for i in cnss])=='N')
            gapcount=np.sum(np.array([i for i in cnss])=='-')
            abase=len(cnss)-(ncount+gapcount)
            matches[gapKey]=dict()
            alignedpos[gapKey]=dict()
            APPLOGGER.info("[{:5,}/{:5,}]\t{}\t[{:>4}|{:>4}|{:>4}]".format(index,totalGaps, f,abase,ncount,gapcount))
            dalign=SeqIO.to_dict(align)
            for iseq in dalign.keys(): 
                for ic in range(0,len(dalign[iseq])):
                    if dalign[iseq][ic]==cnss[ic] and dalign[iseq][ic]!="-":
                        try:
                            matches[gapKey][iseq]+=1
                        except:
                            matches[gapKey][iseq]=1
                    else:
                        matches[gapKey][iseq]=0
                    if dalign[iseq][ic]!="-":
                        try:
                            alignedpos[gapKey][iseq]+=1
                        except:
                            alignedpos[gapKey][iseq]=1
            cnssfile="{}/{}.consensus/{}.fasta".format(OUTPUT,SAMPLENAME,gapKey.replace(":","-"))
            dumbconsensus=SeqIO.SeqRecord(\
                seq=Seq.Seq(str(cnss).replace("-","")),\
                id="{})PLAEC_FillingConsensus".format(gapKey.replace(":","-"))
            )
            sequence=str(dumbconsensus.seq)
            with open(cnssfile,"w") as handle:
                _=SeqIO.write(dumbconsensus,handle,"fasta")
            with open(editfile,'a') as log:
                if len(sequence) <= 2*readFlank:
                    log.write("{:30}\t{:10}\t{:10}\t{}\n".format(scaffold, start,end,""))
                else:
                    log.write("{:30}\t{:10}\t{:10}\t{}\n".format(scaffold, start-readFlank,end+readFlank,sequence))
            with open(logfile2,'a') as log:
                log.write("{:30}\t{:10}\t{:10}\n".format(scaffold, start,end))
            with open(logfile,'a') as log:
                if len(sequence) > 0:
                    for iline in dalign.keys():
                        log.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            gapKey  ,\
                            matches[gapKey][iline],\
                            alignedpos[gapKey][iline],\
                            1.0*matches[gapKey][iline]/alignedpos[gapKey][iline],\
                            len(sequence),\
                            ncount,\
                            gapcount,\
                            abase))
        APPLOGGER.info("Done consensus process > {}".format(datetime.datetime.now()-starttime))
        APPLOGGER.info("="*80)
    except Exception as e:
        APPLOGGER.error("{}\t{}".format(e.__class__, e))
        data=None
    return data, message

###############################################################################
# Adapted from the BioPython function 

def consensus(summary, threshold=.1, ambiguous="X",require_multiple=2):
    consensus = '' 
    # find the length of the consensus we are creating 
    con_len = summary.alignment.get_alignment_length() 
    # go through each seq item 
    for n in range(con_len): 
        # keep track of the counts of the different atoms we get 
        atom_dict = {} 
        num_atoms = 0 
        for record in summary.alignment: 
            # make sure we haven't run past the end of any sequences 
            # if they are of different lengths 
            if n < len(record.seq): 
                if record.seq[n] not in atom_dict: 
                    atom_dict[record.seq[n]] = 1 
                else: 
                    atom_dict[record.seq[n]] += 1 
                num_atoms += 1 
        # APPLOGGER.info(atom_dict)
        max_atoms = [] 
        max_size = 0 
        for atom in atom_dict: 
            if atom_dict[atom] > max_size: 
                max_atoms = [atom] 
                max_size = atom_dict[atom] 
            elif atom_dict[atom] == max_size: 
                max_atoms.append(atom) 
        if require_multiple and num_atoms == 1:
            consensus += ambiguous 
        elif (len(max_atoms) == 1) and ((float(max_size) / float(num_atoms)) >= threshold): 
            consensus += max_atoms[0] 
        else: 
            ambiguous=np.random.choice([it for it in atom_dict.keys() if it !="-"])
            consensus += ambiguous 
    consensus_alpha = summary._guess_consensus_alphabet(ambiguous) 
    return Seq.Seq(consensus, consensus_alpha) 

###############################################################################

def create_output_folders(parameters):
    output=parameters['OUTPUT']
    samplename=parameters['SAMPLENAME']
    outputsupport="{}/{}.support".format(output, samplename)
    APPLOGGER.info("Creating output directory ({}).".format(outputsupport))
    APPLOGGER.warning(\
        "Take into account that if these folders already exist files will be overwritten."
    )
    try:
        os.makedirs(outputsupport)
    except:
        APPLOGGER.info("Output directory ({}) exists.".format(outputsupport))
    try:
        os.makedirs(outputsupport.replace("support","flanks"))
    except:
        APPLOGGER.info("Output directory ({}) exists.".format(outputsupport.replace("support","flanks")))
    try:
        os.makedirs(outputsupport.replace("support","fullsupport"))
    except:
        APPLOGGER.info("Output directory ({}) exists.".format(outputsupport.replace("support","fullsupport")))
    try:
        os.makedirs(outputsupport.replace("support","reads"))
    except:
        APPLOGGER.info("Output directory ({}) exists.".format(outputsupport.replace("support","reads")))
    try:
        os.makedirs(outputsupport.replace("support","pre"))
    except:
        APPLOGGER.info("Output directory ({}) exists.".format(outputsupport.replace("support","pre")))
    try:
        os.makedirs(outputsupport.replace("support","msa"))
    except:
        APPLOGGER.info("Output directory ({}) exists.".format(outputsupport.replace("support","msa")))
    try:
        os.makedirs(outputsupport.replace("support","consensus"))
    except:
        APPLOGGER.info("Output directory ({}) exists.".format(outputsupport.replace("support","consensus")))

###############################################################################

def run(parameters):
    status=0; message="Done"
    create_output_folders(parameters)
    # Get the gap info
    bed=parse_bed_file(parameters['BEDFILE'])
    scaffolds=parse_reference_file(bed, parameters)
    # Getting stats
    # computeIntergapDistance(scaffolds, bed,OUTPUT,SAMPLENAME)
    data=identify_potential_gaps(bed,scaffolds, parameters)
    data, message=extract_support_data(data, parameters, scaffolds)
    if (data):
        data, message=length_agreement(data,scaffolds,parameters)
        if (data):
            data, message=consensus_generation(data,parameters)
        else:
            status = -1
    else:
        status = -1
    return (status, message)

def print_args(args):
    APPLOGGER.info("{}".format("="*80))
    APPLOGGER.info('{:^80}'.format("Settings"))
    APPLOGGER.info("{}".format("="*80))
    APPLOGGER.info("{:>40} = {:<}".format("aln", args.aln))
    APPLOGGER.info("{:>40} = {:<}".format("bed", args.bed))
    APPLOGGER.info("{:>40} = {:<}".format("genome", args.genome))
    APPLOGGER.info("{:>40} = {:<}".format("flanksize", args.flanksize))
    APPLOGGER.info("{:>40} = {:<}".format("output", args.output))
    APPLOGGER.info("{:>40} = {:<}".format("samplename", args.samplename))
    APPLOGGER.info("{:>40} = {:<}".format("mapping_guality_threshold", args.mapping_guality_threshold))
    APPLOGGER.info("{:>40} = {:<}".format("min_coverage_consensus", args.min_coverage_consensus))
    APPLOGGER.info("{:>40} = {:<}".format("min_support", args.min_support))
    APPLOGGER.info("{:>40} = {:<}".format("percent_reads_threshold", args.percent_reads_threshold))
    APPLOGGER.info("{:>40} = {:<}".format("empty_flank_threshold", args.empty_flanks_threshold))
    APPLOGGER.info("{:>40} = {:<}".format("log", args.log))
    APPLOGGER.info("{}".format("="*80))

def main():
    parser = argparse.ArgumentParser(\
        prog=PROGRAM_NAME,\
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
        description='Generates list of potential gaps to be closed or filled.',\
        add_help=True\
    )
    requiredArgs=parser.add_argument_group("{0}Required arguments{1}".format("\033[1m","\033[0m"))
    requiredArgs.add_argument(\
        '-g','--genome',\
        metavar = 'FASTA FILE PATH',\
        type = str,\
        required = True,\
        help = 'Filepath of the reference genome file to be used. Accepts compressed files (GZIP)')
    requiredArgs.add_argument(\
        '-b','--bed',\
        metavar = 'BED FILE PATH',\
        type = str,\
        required = True,\
        help = 'Filepath of the bed file describing the gaps of the genome.'
               ' Accepts compressed files (GZIP)')
    requiredArgs.add_argument(\
        '-a','--aln',\
        metavar  =  'BAM FILE PATH',\
        type  =  str,\
        required =   True,\
        help =  'Filepath of the alignment of reads to the reference genome '
                'in BAM format. This file needs to be indexed before running.')
    requiredArgs.add_argument(\
        '-o','--output',\
        metavar = 'FOLDER_PATH',\
        type = str,\
        required = True,\
        help = 'Output path folder.')
    requiredArgs.add_argument(\
        '-s','--samplename',\
        metavar = 'STR',\
        type = str,\
        required = True,\
        help = 'Short sample name that will be used for naming OUTPUT files.')
    optionalArgs = parser.add_argument_group("{0}Optional arguments{1}".format("\033[1m","\033[0m"))
    optionalArgs.add_argument(\
        '-mins','--min-support',\
        metavar = 'INT',\
        type = int,\
        help =  'Minimum number of reads needed spanning a gap to be considered for closing or filling.',\
        default=5)
    optionalArgs.add_argument(\
        '-f','--flanksize',\
        metavar = 'INT',\
        type = int,\
        help =  'Flank size to be used to select the reads that are in the '
                'surroundings of the gap and determine whether there are '
                'reads that span the gap or not.',\
        default=20)
    optionalArgs.add_argument(\
        '-mapq','--mapping-guality-threshold',\
        metavar = '<MAPQ_threshold>',\
        type = int,\
        help =  'MAPQ value used to filter alignments for posterior processing.'
                'Discarding alingments where: alignment_mapq < MAPQ_threshold.',\
        default=20)
    optionalArgs.add_argument(\
        '-mcc','--min-coverage-consensus',\
        metavar = 'INT',\
        type = int,\
        help =  'Require that more than INT sequences to be part of an alignment to put in the consensus.',\
        default=2)
    optionalArgs.add_argument(\
        '-prt','--percent-reads-threshold',\
        metavar = 'FLOAT',\
        type = float,\
        help =  'Require that more than INT sequences to remain after the length agreement to be considered for consensus. ',\
        default=0.5)
    optionalArgs.add_argument(\
        '-eft','--empty-flanks-threshold',\
        metavar = 'FLOAT',\
        type = float,\
        help =  'Percentage of empty flanks required to skip an ambiguous decision on a gap.',\
        default=0.2)
    optionalArgs.add_argument(\
        '-l','--log',\
        metavar = '<log_level>',\
        type = str,\
        default= "INFO",\
        help = 'Verbosity levels.')
    informationGroup= parser.add_argument_group("{0}Information arguments{1}".format("\033[1m","\033[0m"))
    informationGroup.add_argument('-v', '--version',\
        action='version',\
        version='{0}: Version {1}.{2}.{3}'.format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
        help="Show program's version number and exit")
    try: 
        parameters=dict()
        args = parser.parse_args()      
        print_args(args)
        # Checking files
        if args.log.upper() in "DEBUG":
            APPLOGGER.setLevel(logging.DEBUG)
        else:
            APPLOGGER.setLevel(logging.INFO)
        if (os.path.exists(os.path.abspath(args.genome))):
            parameters['REFFILE']=os.path.abspath(args.genome)
        else:
            message="Reference file does not exist. Please verify. Exiting."
            return (-1, message)
        if (os.path.exists(os.path.abspath(args.bed))):
            parameters['BEDFILE']=os.path.abspath(args.bed)
        else:
            message="BED file does not exist. Please verify. Exiting."
            return (-1, message)
        if (os.path.exists(os.path.abspath(args.aln))):
            parameters['BAMFILE']=os.path.abspath(args.aln)
        else:
            message="BAM file does not exist. Please verify. Exiting."
            return (-1, message)
        parameters['OUTPUT']=os.path.abspath(args.output)
        try:
            APPLOGGER.info("Creating main output directory ({}).".format(parameters['OUTPUT']))
            os.makedirs(parameters['OUTPUT'])
        except:
            APPLOGGER.info("Output directory ({}) exists.".format(parameters['OUTPUT']))
        parameters['FLANKSIZE']=int(args.flanksize)
        parameters['MAPQ_THRESHOLD']=int(args.mapping_guality_threshold)
        parameters['SAMPLENAME']=args.samplename
        parameters['MIN_SUPPORT_THRESHOLD']=int(args.min_support)
        parameters['MIN_CONVERAGE_CONSENSUS']=int(args.min_coverage_consensus)
        parameters['RATIO_THRESHOLD']=int(args.percent_reads_threshold)
        parameters['EMPTY_FLANKS_THRESHOLD']=float(args.empty_flanks_threshold)
        return run(parameters)
    except argparse.ArgumentTypeError as d:
        parser.print_help()
        return -1, "{}\t{}>".format(d.__class__, d.message)
    except Exception as e:
        return -1, "{}\t{}".format(e.__class__, e.message)

if __name__ == "__main__":
    SIGNAL_EXIT,MESSAGE=main()
    if SIGNAL_EXIT == -1:
        APPLOGGER.error(MESSAGE)
    if SIGNAL_EXIT == 2:
        APPLOGGER.warning(MESSAGE)
    else:
        APPLOGGER.info(MESSAGE)
    sys.exit(SIGNAL_EXIT)
