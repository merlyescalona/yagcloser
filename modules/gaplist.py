#! /usr/bin/python
"""
Generates a list of gaps that are potentially to be filled or closed.
2018 - 2019 (c) Merly Escalona <mmescalo@ucsc.edu>  
"""

import datetime,os,csv,glob, argparse, logging,pysam,gzip, filetype, sys
import numpy as np
from Bio import SeqIO, Seq, AlignIO
from StringIO import StringIO
from Bio.Align.Applications import MuscleCommandline, MafftCommandline
import Bio.Align.AlignInfo

# =============================================================================

PROGRAM_NAME="generate_gap_list"
VERSION=1
MIN_VERSION=0
FIX_VERSION=0
# =============================================================================
APPLOGGER = logging.getLogger(PROGRAM_NAME)
APPLOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
loggerFormatter = logging.Formatter(\
  fmt = "[%(asctime)s]  %(levelname)s (%(funcName)s|%(lineno)d): %(message)s",\
  datefmt = "%d/%m/%Y %I:%M:%S %p")
ch.setFormatter(loggerFormatter)
APPLOGGER.addHandler(ch)
# =============================================================================
# # IO Methods

def writeToFile(csvfile, res):
    """ Write table from CSV FILE """
    #Assuming res is a flat list
    with open(csvfile, "w") as OUTPUT:
        writer = csv.writer(OUTPUT, lineterminator='\n')
        writer.writerows(res)

def readFromFile(csvfile):
    """ Read table from CSV FILE """
    #Assuming res is a flat list
    res=[]
    with open(csvfile) as csv_file:
        csvreader = csv.reader(csv_file, delimiter=',')
        for row in csvreader:
            res+=[[float(i) for i in row]]
    return np.array(res)

# =============================================================================

def getBed(bedfile):
    """ Reading bedfile to extract gap info"""
    APPLOGGER.info("Reading BED file ({})".format(bedfile))
    # print("Reading BED file ({})".format(bedfile))
    start=datetime.datetime.now()
    bed=dict()
    # BED Dictionary:
    # BED[SCAFFOLD]: GAPID,START,END,SIZE
    bedCounter=1
    with open(bedfile) as f:
        for line in f:
            bedline=line.strip().split()
            try:
                bed[bedline[0]]+=[(bedCounter,int(bedline[1]), int(bedline[2]), int(bedline[2])-int(bedline[1]))]     
            except:
                bed[bedline[0]]=[(bedCounter,int(bedline[1]), int(bedline[2]), int(bedline[2])-int(bedline[1]))]     
            bedCounter+=1
        end=datetime.datetime.now()
    APPLOGGER.info("Done reading BED file > {}".format( end-start))
    # print("Done reading BED file > {}".format( end-start))
    return bed

# =============================================================================

def getRefereceScaffolds(reffile,bed):
    """ Get reference file info"""
    APPLOGGER.info("Reading scaffolds from reference file...")
    # print("Reading scaffolds from reference file...")
    start=datetime.datetime.now()
    scaffolds=dict()
    counter=0
    # ScpqPdg
    filekind=filetype.guess(reffile)
    handle=None
    if filekind and filekind.extension in ['gz','GZ']:
        handle=gzip.open(reffile, 'rb') 
    else:
        handle=open(reffile, "rt")
    for record in SeqIO.parse(handle, "fasta"):
        if (counter%10)==0:APPLOGGER.debug("Scaffold:\t{}".format(counter))
        # if (counter%10)==0: print("Scaffold:\t{}".format(counter))
        newid=record.id.split("|")[0]
        try:
            scaffolds[newid]={'id':counter,'seq':str(record.seq),'size':len(record.seq), 'flanks':[], 'numgaps':len(bed[newid])}  
        except:
            scaffolds[newid]={'id':counter,'seq':str(record.seq),'size':len(record.seq), 'flanks':[], 'numgaps':0}
        counter+=1
    handle.close()
    end=datetime.datetime.now()
    # Sorting the scaffold by LENGTH
    scaffoldnames=scaffolds.keys()
    scaffoldnames=[s.split(";") for s in scaffoldnames]
    scaffoldnames=[s[0].split("_")+s[1].split("=") for s in scaffoldnames]
    scaffoldnames=[ [s[0], int(s[1]), s[2], int(s[3])] for s in scaffoldnames]
    scaffoldnames=sorted(scaffoldnames, key=lambda k: k[1])
    scaffoldsOrder=[ "{}_{};{}={}".format(i[0], i[1], i[2], i[3]) for i in scaffoldnames]
    for index in range(0,len(scaffoldsOrder)):
        s=scaffoldsOrder[index]
        scaffolds[s]['id']=index
    APPLOGGER.info("Done reading reference file ({}) > {}".format(reffile, end-start))
    # print("Done reading reference file ({}) > {}".format(reffile, end-start))
    return scaffolds

# =============================================================================

def getFlankTable(scaffolds, bed, READFLANK):
    APPLOGGER.info("Extracting flank regions...")
    # print("Extracting flank regions...")
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
    # print("Done with flanks extraction: > {}".format(end-start))
    return flankTable

# =============================================================================

def computeIntergapDistance(scaffolds, bed,OUTPUT,SAMPLENAME):
    start=datetime.datetime.now()
    intergapdistancefile='{}/{}.gap.interdistance.txt'.format(OUTPUT,SAMPLENAME)
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
    print("Intergap distances... ({}) > {}".format(intergapdistancefile, end-start))

# =============================================================================

def readLengthDistribution(pathpattern):
    filelist=glob.glob(pathpattern)
    readlen=[]
    for item in filelist:
        with open(item,'r') as handle:
            readlen+=[ [n.strip().split()[0]]*int(n.strip().split()[1]) for n in handle.readlines()]

# =============================================================================

def identify_potential_fillable_gaps(\
    bed, scaffolds,BAMFILE, READFLANK,OUTPUT,SAMPLENAME,\
    MAPQ_THRESHOLD,MIN_SUPPORT_THRESHOLD):
    APPLOGGER.info("Data structure initialization... ")
    # print("Data structure initialization... ")
    gapCounter=0
    start=datetime.datetime.now()
    numGaps=len([j for i in bed.keys() for j in range(0,len(bed[i]))])
    data=dict()
    for scaffold in bed.keys():
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
                APPLOGGER.debug("Gap {:4}/{}:  || {:>10} | {:<10} ||".format(\
                # print("Gap {:4}/{}:  || {:>10} | {:<10} ||".format(\
                    gapCounter,\
                    numGaps,\
                    gapStart,\
                    gapEnd\
                ))
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
    APPLOGGER.info("Done with initialization... > {}".format( end-start))
    # print("Done with initialization... > {}".format( end-start))
    ######################################################################################
    APPLOGGER.info("Extracting SAM/BAM information... ")
    # print("Extracting SAM/BAM information... ")
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
            side=""
            if gapCounter % 100 == 0:
                APPLOGGER.debug("Gap {:4}/{} ({readFlank:}): [{:<25}] ||{:>10}|{:<10}||".format(\
                # print("Gap {:4}/{} ({readFlank:}): [{:<25}] ||{:>10}|{:<10}||".format(\
                    gapCounter,\
                    numGaps,\
                    scaffold,\
                    gapStart,\
                    gapEnd,\
                    readFlank=READFLANK\
                ))
            # Checking only PRIMARY alignments, with MAPQ > 20 and no DUPLICATES
            for read in samfile.fetch(region=region):
                data[regionKey]['scaffold']=scaffold
                data[regionKey]['alncov']+=1
                data[regionKey]['gapstart']=int(gapStart)
                data[regionKey]['gapend']=int(gapEnd)
                data[regionKey]['gapsize']=int(gapSize)
                if (not read.is_secondary) and (read.mapping_quality > MAPQ_THRESHOLD) and (not read.is_duplicate):            
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
    APPLOGGER.info("Done with the info extraction... > {}".format( end-start))
    # print("Done with the info extraction... > {}".format( end-start))
    # =============================================================================
    # Removing 0-support data
    totalGaps=len(data.keys())
    APPLOGGER.info("Removing gaps with no support from analysis... ")
    # print("Removing gaps with no support from analysis... ")
    start=datetime.datetime.now()
    removedGaps=[]
    removedCounter=0
    gapCounter=0
    for gapKey in data.keys():
        if gapCounter % 100 == 0:   
            APPLOGGER.debug("Gap {:4}/{}: {:35}".format(gapCounter,totalGaps,gapKey)) 
            # print("Gap {:4}/{}: {:35}".format(gapCounter,totalGaps,gapKey)) 
        if len(data[gapKey]['mapq']) ==  0: 
            removedGaps+=[(gapKey,data[gapKey]['gapsize'] )]
            removedCounter+=1    
            del data[gapKey]
        gapCounter+=1
    end=datetime.datetime.now()
    APPLOGGER.info("Removed {} gaps.. > {}".format(removedCounter, end-start))
    # print("Removed {} gaps.. > {}".format(removedCounter, end-start))
    APPLOGGER.info("Writing removed gaps to file: {}/{}.no.support.gaps.txt".format(OUTPUT,SAMPLENAME))
    # print("Writing removed gaps to file: {}/{}.no.support.gaps.txt".format(OUTPUT,SAMPLENAME))
    with open('{}/{}.no.support.gaps.txt'.format(OUTPUT,SAMPLENAME), "w") as handle:
        for pair in removedGaps:
            gap=pair[0]; size=pair[1]
            handle.write("{}\t{}\n".format(gap, size))
    # =============================================================================
    gapCounter=0
    totalGaps=len(data.keys())
    gapAlnCoverage=dict()
    gapReadCoverage=dict()
    gapAlnCoverageB=dict()
    counterB=dict()
    counterL=dict()
    counterR=dict()
    for k in data.keys():
        gapAlnCoverage[k]=0
        gapReadCoverage[k]=0
        gapAlnCoverageB[k]=0
        counterB[k]=0
        counterR[k]=0
        counterL[k]=0
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
            counterB[gapKey]+=valB
            counterL[gapKey]+=valL
            counterR[gapKey]+=valR
            gapAlnCoverageB[gapKey]+=min(1, min(valR,valL))
    #=========================================================================================
    toBeKeptBothSides={i:counterB[i] for i in counterB if counterB[i] > MIN_SUPPORT_THRESHOLD }
    #=========================================================================================
    gapCounter=0
    for gapKey in toBeKeptBothSides:
        if gapCounter % 100 == 0:   
            APPLOGGER.debug("Gap {:4}/{}: {:35}".format(gapCounter,totalGaps,gapKey)) 
        gapCounter+=1
        keys=data[gapKey]['side'].keys()
        for readKey in keys:
            indexL=list(np.where(np.array([item for item in data[gapKey]['side'][readKey]])=="L")[0])
            indexR=list(np.where(np.array([item for item in data[gapKey]['side'][readKey]])=="R")[0])
            alnToRemove=sorted(indexL+indexR, reverse=True)
            # print(data[gapKey]['side'][readKey],alnToRemove)
            for ind in alnToRemove:
                del data[gapKey]['side'][readKey][ind]
                del data[gapKey]['mapq'][readKey][ind]
                del data[gapKey]['supp'][readKey][ind]
                del data[gapKey]['primary'][readKey][ind]
                del data[gapKey]['positions'][readKey][ind]
                del data[gapKey]['refpositions'][readKey][ind]
                del data[gapKey]['querylen'][readKey][ind]
                del data[gapKey]['AS'][readKey][ind]
            if len(data[gapKey]['side'][readKey]) == 0:
                del data[gapKey]['side'][readKey]
    APPLOGGER.info("Writing into file potential fillable gaps ({}/{}.potential.fillable.gaps.txt)...".format(OUTPUT,SAMPLENAME))
    # print("Writing into file potential fillable gaps ({}/{}.potential.fillable.gaps.txt)...".format(OUTPUT,SAMPLENAME))
    with open('{}/{}.potential.fillable.gaps.txt'.format(OUTPUT,SAMPLENAME), "w") as handle:
        handle.write("{:25}\t{:15}\t{:15}\t{}\t{}\n".format("scaffold","start", "end","gapsize","support"))
        for gapKey in toBeKeptBothSides:
            handle.write("{:25}\t{:12}\t{:12}\t{}\t{}\n".format(data[gapKey]['scaffold'],data[gapKey]['gapstart'],data[gapKey]['gapend'],data[gapKey]['gapsize'], len(data[gapKey]['side'])))
    data={gapKey:data[gapKey] for gapKey in toBeKeptBothSides}
    return data

def generateSupportDataFolder(data, BAMFILE, OUTPUT, SAMPLENAME,READFLANK, scaffolds):
    try:
        allreads=dict()
        fullreads=dict()
        APPLOGGER.info("Generating support files")
        samfile = pysam.AlignmentFile(BAMFILE, "rb")
        gapCounter=1
        APPLOGGER.info("Total before filtering {}".format(len(data)))
        APPLOGGER.info("Filtering oversampling (support > 1000)")
        data2={ sc:data[sc] for sc in data.keys() if len(data[sc]['side']) < 1000 }
        APPLOGGER.info("Total after filtering {}".format(len(data2)))
        data=data2
        index=0
        for gapKey in data:
            # APPLOGGER.debug("[{:>10,}/{:>10,}]\t{}".format(gapCounter, len(data),gapKey))
            scaffold=gapKey.split(":")[0]
            start=np.max([int(gapKey.split(":")[1].split("-")[0])-READFLANK,1])
            end=np.min([int(gapKey.split(":")[1].split("-")[1])+READFLANK,len( scaffolds[scaffold]['seq'])])
            if (index%100 ==0):
                APPLOGGER.debug("[{:>10,}/{:>10,}]\t({}){}-{}".format(gapCounter, len(data),gapKey, start,end))
            index+=1
            allreads[gapKey]={}
            fullreads[gapKey]={}
            readlist=[]
            regionKey="{}:{}-{}".format(scaffold,start,end)
            iter=samfile.pileup(contig=scaffold, start=start, end=end)
            for x in iter: 
                if x.reference_pos > start and x.reference_pos < end:
                    readlist+=x.get_query_names()
            readlist=[r for r in  list(np.unique(readlist)) if r in data[gapKey]['side']]
            allreads[gapKey]={rn:{"positions":[], "sequence":"", "reference":"", "read":"", "lflank":"","rflank":""} for rn in readlist}
            fullreads[gapKey]={rn:{"positions":[], "sequence":"", "read":"", "lflank":"","rflank":""} for rn in readlist}
            for x in samfile.pileup(contig=scaffold, start=start, end=end):
                if x.reference_pos > start and x.reference_pos < end:
                    reads=x.get_query_names()
                    positions=x.get_query_positions()
                    sequences=x.get_query_sequences()
                    for rindex in range(0, len(reads)):
                        if reads[rindex] in readlist:
                            allreads[gapKey][reads[rindex]]['positions']+=[int(positions[rindex])]
                            # APPLOGGER.info("{}\t{}\t{}".format(gapKey, x.reference_pos, reads[rindex]))
                            if x.reference_pos < start+READFLANK or x.reference_pos >= end-READFLANK:
                                allreads[gapKey][reads[rindex]]['sequence']+=sequences[rindex].upper()
                                allreads[gapKey][reads[rindex]]['reference']+=scaffolds[scaffold]['seq'][x.reference_pos].upper()
                                if x.reference_pos <start+READFLANK:
                                    allreads[gapKey][reads[rindex]]['lflank']+=scaffolds[scaffold]['seq'][x.reference_pos].upper()                      
                                    fullreads[gapKey][reads[rindex]]['lflank']+=sequences[rindex]
                                if x.reference_pos >= end-READFLANK:
                                    allreads[gapKey][reads[rindex]]['rflank']+=scaffolds[scaffold]['seq'][x.reference_pos].upper()
                                    fullreads[gapKey][reads[rindex]]['rflank']+=sequences[rindex]                              
                            else:
                                allreads[gapKey][reads[rindex]]['sequence']+=sequences[rindex].lower()
                                allreads[gapKey][reads[rindex]]['reference']+=scaffolds[scaffold]['seq'][x.reference_pos].lower()
                        # APPLOGGER.info("{}\t{}\t{}\tPILEUPS".format(gapKey, x.reference_pos, reads[rindex]))
                    for pile in x.pileups:
                        readname=pile.alignment.query_name
                        if readname in readlist:
                            minpos=np.min(allreads[gapKey][readname]['positions'])
                            maxpos=np.max(allreads[gapKey][readname]['positions'])
                            allreads[gapKey][readname]['read']=pile.alignment.query_sequence[minpos-1:maxpos+1]
                            # second part
                            fullreads[gapKey][readname]['read']=pile.alignment.query_sequence
                            if pile.query_position:
                                fullreads[gapKey][readname]['positions']+=[pile.query_position]
            # APPLOGGER.debug("[{:>10,}/{:>10,}]\tGetting positions".format(gapCounter, len(data)))
            for readKey in fullreads[gapKey]:
                if len(fullreads[gapKey][readKey]['positions'])>0:
                    minpos=min(fullreads[gapKey][readKey]['positions'])
                    maxpos=max(fullreads[gapKey][readKey]['positions'])
                    fullreads[gapKey][readKey]['sequence']=fullreads[gapKey][readKey]['read'][minpos-1:maxpos+1]
            gapCounter+=1
        samfile.close()
        APPLOGGER.info("Writing support files")
        gapCounter=0
        for gapKey in data:
            APPLOGGER.debug("[{:>10,}/{:>10,}]\tWriting {}...".format(gapCounter,len(data),gapKey))
            gapCounter+=1
            gapFilename0="{}/{}.support/{}.fasta".format(OUTPUT, SAMPLENAME,gapKey.replace(";", "-").replace(":", "-").replace("=", "-"))
            gapFilename1="{}/{}.fullsupport/{}.fasta".format(OUTPUT, SAMPLENAME,gapKey.replace(";", "-").replace(":", "-").replace("=", "-"))
            gapFilename2="{}/{}.flanks/{}.fasta".format(OUTPUT, SAMPLENAME,gapKey.replace(";", "-").replace(":", "-").replace("=", "-"))
            gapFilename3="{}/{}.pileup/{}.fasta".format(OUTPUT, SAMPLENAME,gapKey.replace(";", "-").replace(":", "-").replace("=", "-"))
            MSA=[]
            MSA2=[]
            FLANKS=[]
            PILE=[]
            for readKey in allreads[gapKey]:
                F1=fullreads[gapKey][readKey]['lflank']
                F2=fullreads[gapKey][readKey]['rflank']
                if len(allreads[gapKey][readKey]['read'][len(F1):-len(F2)]) >0:
                    APPLOGGER.debug("{}\t{}\t{}".format(\
                        gapKey,\
                        readKey,\
                        len(allreads[gapKey][readKey]['sequence'])))
                    MSA+=[
                        SeqIO.SeqRecord(
                            seq=Seq.Seq(allreads[gapKey][readKey]['sequence'][len(F1):-len(F2)]), \
                            id=readKey\
                        )
                    ]
                    MSA2+=[
                        SeqIO.SeqRecord(
                            seq=Seq.Seq(fullreads[gapKey][readKey]['sequence']), \
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
                            seq=Seq.Seq(fullreads[gapKey][readKey]['sequence'][len(F1):-len(F2)]), \
                            id="{}_PILEUP_NO_FLANKS".format(readKey)\
                        ),\
                        SeqIO.SeqRecord(
                            seq=Seq.Seq(fullreads[gapKey][readKey]['sequence']), \
                            id="{}_PILEUP_FLANKS".format(readKey)\
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
    except AssertionError as ae:
        return -1, ae.message
    except Exception as e:
        return -1, e.message

def generateOutputFolders(output, samplename):
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
        os.makedirs(outputsupport.replace("support","pileup"))
    except:
        APPLOGGER.info("Output directory ({}) exists.".format(outputsupport.replace("support","pileup")))



def run(REFFILE,BEDFILE,OUTPUT,BAMFILE,SAMPLENAME, READFLANK,MAPQ_THRESHOLD, MIN_SUPPORT_THRESHOLD):
    # Get the gap info
    bed=getBed(BEDFILE)
    numGaps=len([j for i in bed.keys() for j in range(0,len(bed[i]))])
    scaffolds=getRefereceScaffolds(REFFILE,bed)
    flankTable=getFlankTable(scaffolds,bed,READFLANK)
    # Updating scaffolds data structure
    for f in flankTable:
        scaffolds[f[2]]['flanks']+=[f[1]]
    for f in flankTable:
        scaffolds[f[2]]['numgaps']=len(scaffolds[f[2]]['flanks'])/2
    # Getting stats
    computeIntergapDistance(scaffolds, bed,OUTPUT,SAMPLENAME)
    generateOutputFolders(OUTPUT,SAMPLENAME)
    data=identify_potential_fillable_gaps(\
        bed,scaffolds, BAMFILE,READFLANK,OUTPUT,SAMPLENAME,\
        MAPQ_THRESHOLD, MIN_SUPPORT_THRESHOLD)
    status, message=generateSupportDataFolder(data, BAMFILE, OUTPUT, SAMPLENAME,READFLANK, scaffolds)
    if status==0:
        APPLOGGER.info("Done")
        return (0, "")
    else:
        return (status, message)

def main():
    parser = argparse.ArgumentParser(\
        prog="potential_fillable_gaps",\
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
        description='Generates list of potential gaps to be closed or filled.',\
        add_help=False\
    )
    requiredArgs=parser.add_argument_group("{0}Required arguments{1}".format("\033[1m","\033[0m"))
    requiredArgs.add_argument(\
        '-g','--genome',\
        metavar = '<assembly.fasta>',\
        type = str,\
        required = True,\
        help = 'Filepath of the reference genome file to be used.')
    requiredArgs.add_argument(\
        '-b','--bed',\
        metavar = '<gaps.bed>',\
        type = str,\
        required = True,\
        help = 'Filepath of the bed file describing the gaps of the genome')
    requiredArgs.add_argument(\
        '-a','--aln',\
        metavar  =  '<alignment.bam>',\
        type  =  str,\
        required =   True,\
        help =  'Filepath of the alignment of reads to the reference genome'
                'in BAM format. This file needs to be indexed before running this.')
    requiredArgs.add_argument(\
        '-o','--OUTPUT',\
        metavar = '<path_folder>',\
        type = str,\
        required = True,\
        help = 'Output path folder.')
    requiredArgs.add_argument(\
        '-s','--samplename',\
        metavar = '<sample_name>',\
        type = str,\
        required = True,\
        help = 'Short sample name that will be used for naming OUTPUT files.')
    optionalArgs = parser.add_argument_group("{0}Optional arguments{1}".format("\033[1m","\033[0m"))
    optionalArgs.add_argument(\
        '-mins','--min-support',\
        metavar = '<MIN_SUPPORT_THRESHOLD>',\
        type = int,\
        help =  'Minimum number of reads spanning a specific gap to be considered',\
        default=5)
    optionalArgs.add_argument(\
        '-f','--readflank',\
        metavar = 'FLANK_SIZE',\
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
                'Discarding alingments where: alignment.mapq < MAPQ_threshold.',\
        default=20)
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
    informationGroup.add_argument('-h', '--help',\
        action='store_true',\
        help="Show this help message and exit")
    try: 
        args = parser.parse_args()      
        APPLOGGER.info(args)
        if args.help:
            parser.print_help()
        # Checking files
        if args.log.upper() in "DEBUG":
            APPLOGGER.setLevel(logging.DEBUG)
        else:
            APPLOGGER.setLevel(logging.INFO)
        if (os.path.exists(os.path.abspath(args.genome))):
            REFFILE=os.path.abspath(args.genome)
        else:
            message="Reference file does not exist. Please verify. Exiting."
            return (-1, message)
        if (os.path.exists(os.path.abspath(args.bed))):
            BEDFILE=os.path.abspath(args.bed)
        else:
            message="BED file does not exist. Please verify. Exiting."
            return (-1, message)
        if (os.path.exists(os.path.abspath(args.aln))):
            BAMFILE=os.path.abspath(args.aln)
        else:
            message="BAM file does not exist. Please verify. Exiting."
            return (-1, message)
        OUTPUT=os.path.abspath(args.OUTPUT)
        try:
            APPLOGGER.info("Creating OUTPUT directory ({}).".format(OUTPUT))
            # print("Creating OUTPUT directory ({}).".format(OUTPUT))
            os.makedirs(OUTPUT)
        except:
            APPLOGGER.info("Output directory ({}) exists.".format(OUTPUT))
            # print("Output directory ({}) exists.".format(OUTPUT))
        READFLANK=int(args.readflank)
        MAPQ_THRESHOLD=int(args.mapping_guality_threshold)
        SAMPLENAME=args.samplename
        MIN_SUPPORT_THRESHOLD=int(args.min_support)
        return run(REFFILE,BEDFILE,OUTPUT,BAMFILE,SAMPLENAME, READFLANK,MAPQ_THRESHOLD, MIN_SUPPORT_THRESHOLD)
    except argparse.ArgumentTypeError as d:
        APPLOGGER.info(args)
        parser.print_help()
        return (-1,d.message)
    except Exception as e:
        return (-1,e.message)

if __name__ == "__main__":
    SIGNAL_EXIT,MESSAGE=main()
    if SIGNAL_EXIT == -1:
        APPLOGGER.error(MESSAGE)
    else:
        APPLOGGER.info(MESSAGE)
    sys.exit(SIGNAL_EXIT)



# args=argparse.Namespace(\
#       OUTPUT='/scratch1/merly/Gibbon/old-gapclosing/',
#       aln='/scratch5/gibbon/nanopore/1D_Gibbon_Genomic.s.bam',\
#       bed='/scratch1/merly/Gibbon/gaps/gibbon.2.gaps.bed',\
#       genome='/scratch1/merly/Gibbon/gibbon.assembly.2.fasta',\
#       help=False,\
#       log='DEBUG',\
#       mapping_guality_threshold=20,\
#       min_support=2,\
#       readflank=100,\
#       samplename='gibbon.2.ont')