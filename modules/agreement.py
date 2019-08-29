import pysam,gzip,datetime,csv,os, glob,collections, copy,argparse
import numpy as np
from Bio import SeqIO, Seq, AlignIO, Alphabet
from StringIO import StringIO
from Bio.Align.Applications import MuscleCommandline, MafftCommandline
import Bio.Align.AlignInfo

identityThreshold=0.8
ratioThreshold=0.5      # number of reads selected from the supporting set of reads
coverageConsensus=2


samplename="luhya.2.pb"
output="/scratch4/NA19434/gapclosing"
logfile="/scratch4/NA19434/gapclosing/{}.confirm.log.txt".format(samplename)
editfile="/scratch4/NA19434/gapclosing/{}.edit.file.support.txt".format(samplename)


fastafiles=glob.glob("{}/{}.support/*".format(output, samplename))
# Re-populate list with filename, size tuples


for i in xrange(len(fastafiles)):
    fastafiles[i] = (fastafiles[i], os.path.getsize(fastafiles[i]))


fastafiles.sort(key=lambda filename: filename[1])

tmp=[ os.path.basename(fastafiles[ff][0]).split(".")[0] for ff in range(0,len(fastafiles))]

gapsToBeFilled=[]
for i in tmp:
    tt=i.split("-")
    gapsToBeFilled+=["{};{}={}:{}-{}".format(tt[0],tt[1],tt[2], tt[3], tt[4])]
    # gapsToBeFilled+=["{}:{}-{}".format(tt[0],tt[1],tt[2])]




outputsupport="{}/{}.pre".format(output, samplename)
try:
    os.makedirs(outputsupport)
except:
    print("Output directory ({}) exists.".format(outputsupport))


try:
    os.makedirs(outputsupport.replace("pre","msa"))
except:
    print("Output directory ({}) exists.".format(outputsupport.replace("pre","msa")))

try:
    os.makedirs(outputsupport.replace("pre","consensus"))
except:
    print("Output directory ({}) exists.".format(outputsupport.replace("pre","consensus")))




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
        # print(atom_dict)
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

seqlengths=dict()
# for gapKey in gapsToBeFilled:
for gapIndex in range(0,len(gapsToBeFilled)):
    gapKey=gapsToBeFilled[gapIndex]
    print(gapKey)
    illuminafile="{}/{}.support/{}.fasta".format(output,samplename,gapKey.replace(";","-").replace(":","-").replace("=","-") )
    with open(illuminafile, "rt") as handle:
        reads=SeqIO.to_dict(SeqIO.parse(handle, "fasta")) 
    MSA=[reads[record] for record in reads]
    print("{}\tMSA".format(gapKey))    
    #=====================================================================================    
    # print(gapKey,len(MSA))
    #=====================================================================================
    aligner="NONE"
    timing=0
    # Sorting of the sequences based on length (longest first) 
    MSA.sort(key = lambda s: len(s.seq), reverse=True)
    gapalnfilename="{}/{}.pre/{}.fasta".format(output,samplename,gapKey.replace(";","-").replace(":","-").replace("=","-"))
    # with open(gapalnfilename,"w") as handle:
    #     _=SeqIO.write(MSA,handle,"fasta")
    seqlengths[gapKey]=[len(i)for i in MSA]
    print("{}\tMSA Writing".format(gapKey))
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
        print("{}\tBefore while loop".format(gapKey))     
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
        print("{}\tDone length agreement".format(gapKey))    
        with open(logfile,"a") as handle: 
            handle.write("{}\t{}\t{}\n".format(gapKey, len(MSA)+len(selectedPositions), len(MSA)))
        if len(MSA) > 1 and ratio >= ratioThreshold:
            with open(gapalnfilename, "w") as handle: 
                _=SeqIO.write(MSA,handle,"fasta")
            #=====================================================================================        
            # Also I can do the alignment
            msafile="{}/{}.msa/{}.fasta".format(output,samplename,gapKey.replace(";","-").replace("=","-").replace(":","-"))
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
            # align = AlignIO.read(StringIO(stdout), "fasta")
            endMUSCLE=datetime.datetime.now()
            lenMSA=len(MSA)
            timing=endMUSCLE-startMUSCLE
            print("{:50}\t{}".format(gapKey,timing))



readFlank=20
msafiles=glob.glob('{}/{}.msa/*.fasta'.format(output,samplename))
msafiles.sort()
logfile="{}/{}.consensus.log.txt".format(output,samplename)
matches=dict()
alignedpos=dict()
for index in range(0,len(msafiles)):
    f=msafiles[index]
    filebasename=os.path.basename(f)[0:-6]
    s=filebasename.split("-")
    s=s[0:2]+s[2:6]
    scaffold="{scaffold[0]};{scaffold[1]}={scaffold[2]}".format(scaffold=s[0:3])
    start=int(s[3])
    end=int(s[4])
    # scaffold=s[0]
    # start=s[1]
    # end=s[2]
    align=AlignIO.read(f, "fasta")
    summary=Bio.Align.AlignInfo.SummaryInfo(align)
    cnss=consensus(summary,threshold=0.1,ambiguous='N', require_multiple=coverageConsensus)
    ncount=np.sum(np.array([i for i in cnss])=='N')
    gapcount=np.sum(np.array([i for i in cnss])=='-')
    abase=len(cnss)-(ncount+gapcount)
    matches[filebasename]=dict()
    alignedpos[filebasename]=dict()
    print("[{:5}]\t{}\t[{:>4}|{:>4}|{:>4}]".format(index, f,abase,ncount,gapcount))
    dalign=SeqIO.to_dict(align)
    for iseq in dalign.keys(): 
        for ic in range(0,len(dalign[iseq])):
            if dalign[iseq][ic]==cnss[ic] and dalign[iseq][ic]!="-":
                try:
                    matches[filebasename][iseq]+=1
                except:
                    matches[filebasename][iseq]=1
            else:
                matches[filebasename][iseq]=0
            if dalign[iseq][ic]!="-":
                try:
                    alignedpos[filebasename][iseq]+=1
                except:
                    alignedpos[filebasename][iseq]=1
    cnssfile="{}/{}.consensus/{}".format(output,samplename,os.path.basename(f))
    dumbconsensus=SeqIO.SeqRecord(\
        seq=Seq.Seq(str(cnss).replace("-","")),\
        id="{})PLAEC_FillingConsensus".format(os.path.basename(f))
    )
    sequence=str(dumbconsensus.seq)
    with open(cnssfile,"w") as handle:
        _=SeqIO.write(dumbconsensus,handle,"fasta")
    with open(editfile,'a') as log:
        # log.write("{:30}\t{:10}\t{:10}\t{}\n".format(scaffold, start,end,sequence))
        if len(sequence) <= 2*readFlank:
            log.write("{:30}\t{:10}\t{:10}\t{}\n".format(scaffold, start,end,""))
        else:
            log.write("{:30}\t{:10}\t{:10}\t{}\n".format(scaffold, start-readFlank,end+readFlank,sequence))
    with open(logfile,'a') as log:
        if len(sequence) > 0:
            for iline in dalign.keys():
                log.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    filebasename,\
                    matches[filebasename][iline],\
                    alignedpos[filebasename][iline],\
                    1.0*matches[filebasename][iline]/alignedpos[filebasename][iline],\
                    len(sequence),\
                    ncount,\
                    gapcount,\
                    abase))







