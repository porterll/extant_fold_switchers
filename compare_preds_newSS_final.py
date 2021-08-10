#! /usr/bin/python

#Note: this code uses python2.7

#This code compares predicted and experimentally observed secondary structures
#Predictions were made using JPred4.
#It was used to generate Figure 1, Figure S1, and Table S1 in
#Mishra, et. al., Biopolymers 2021:
#A sequence-based method for predicting extant fold switchers that
#undergo alpha-helix<-> beta-strand transitions
#Inputs:
#        Log3 files, which compare 2 experimentally determined secondary
#        structures (taken from Porter and Looger (2018) PNAS)
#        JPred Secondary Structure predictions of:
#        (1) Full-length PDBs taken from Misra, et. al. (2019) Protein Science
#        (2) FSRs generated for this publication
#        List of fasta files used to generate predictions (fastas.txt)
#        List of PDB IDs corresponding to fold-switch pairs(updated_pdb_ids.csv)
#Contact: Lauren Porter (porterll@nih.gov)

import sys, string
from Bio import pairwise2
from matplotlib import pyplot as plt
from difflib import SequenceMatcher
import numpy as np
from scipy import stats
import cPickle

#Input directories
#Predicted secondary structures of FSRs
CDIR = 'seqs_together/'
#Predicted secodary structures of whole PDBs
RDIR = 'JPRED_wholes/'
#Comparisons of experimentally determined secondary structures
LDIR = 'important_LOG3s/'

cdict = {'H':'H','E':'E','-':'C'}
cdict2 = {'H':'H','E':'E',' ':'C','T':'C','S':'C','G':'C','B':'C','I':'C'}

class PDBINFO:

    def __init__(self,info,indices,sequence,ss):

        self.ID          = None
        self.chain       = None
        self.length      = None
        self.description = None
        self.indices     = None
        self.sequence    = None
        self.ss          = None

        self.get_pdb_info(info,indices,sequence,ss)

        self.resis2include = None

    def get_pdb_info(self,info,indices,sequence,ss):

        self.ID = info[0][1:5]
        self.chain = info[0][-1]

        if len(info) == 1:
            
            self.length = None
            self.description = None

        else:
            self.length = int(info[2][7:])
            self.description = string.join(info[3:],' ')
            
        self.indices = (int(indices.split()[1]),int(indices.split()[3]))
        self.sequence = sequence
        self.ss = ss

    def _print(self):

        print 'ID= ',self.ID
        print 'Chain= ',self.chain
        print 'length= ',self.length
        print 'description= ',self.description
        print 'indices= ', self.indices
        print 'sequence= ', self.sequence
        print 'ss= ', self.ss

def max_str(diffstring,maxDiff,sequence):

    l = len(maxDiff)

    if diffstring==maxDiff:
        return sequence, 0, len(sequence)

    for i in xrange(len(diffstring)-l+1):
        if maxDiff == diffstring[i:i+l]:
            return sequence[i:i+l], i, i+l


    return -999,-999,-999

def min_str(str1,str2):

    s1 = string.join([x for x in str1 if x!='-'],'');
    s2 = string.join([x for x in str2 if x!='-'],'');

    if len(s1)<=len(s2):
        return s1

    return s2

#Read in Secondary Structure comparison files
#Generated from Porter and Looger (2018) PNAS

class LOG3:

    def __init__(self,f):

        self.template   = None
        self.pdbInfo    = []
        self.identities = []
        self.hammDists  = []
        self.hammDiffs  = []
        self.maxSeqs    = []
        self.maxStrings = []
        self.diffStrings= []
        self.error1 = 0

        self.parse(f)

        self.idxs = {}

        for i in xrange(len(self.pdbInfo)):
            self.idxs[self.pdbInfo[i].ID+self.pdbInfo[i].chain] = i

        self.bseqs = []
        
        self.fss = []
        self.rfs = []
        self.fps = []

    def parse(self,fname):
        
        i = 0

        f = open(fname).read().splitlines()

        while i < len(f):

            info = f[i].split()

            if info[0] == 'Template:':
                self.template = PDBINFO(f[i+2].split()[:-1],f[i+5],f[i+3],
                                        f[i+4])
                i += 8
                continue
            
            if info[0][0] == '>':
                
                pdb1 = PDBINFO(info,f[i+3],f[i+1],f[i+2])

                sequence = ''
                
                self.pdbInfo.append(pdb1)
                self.identities.append(SequenceMatcher(None,
                                                       self.template.sequence,
                                                       pdb1.sequence).ratio())
                
                self.hammDists.append(float(f[i+4].split()[2][:-1]))
                self.hammDiffs.append(float(f[i+4].split()[5][:-1]))
                self.maxStrings.append(f[i+4].split()[-4][:-1])
                self.diffStrings.append(f[i+5].split()[0])

                sequence,i1,i2 = max_str(self.diffStrings[-1],
                                         self.maxStrings[-1],pdb1.sequence)

                self.pdbInfo[-1].indices = (i1,i2)

                self.maxSeqs.append((f[i+4].split()[-1],sequence))
                
                i+=8

            else:
                print 'File parsing error %s' %(fname)
                self.error1 = 1
                i += 1

def convert(f):

    s= f[1][9:]
    return string.join([cdict[x] for x in s if x !=','],'')

def SS(l):

    return float(len([x for x in l if x in ['H','E']]))

def get_ldict(f):

    ldict = {}
    IDdict = {}

    for i in f:
        info = string.split(i,',')
        fname = string.split(info[0],'/')[1].strip()
        ldict[info[1].strip()] = fname
        ldict[info[2].strip()] = fname
        IDdict[info[1].strip()] = info[2].strip()
        IDdict[info[2].strip()] = info[1].strip()

    return ldict, IDdict

def get_maxString(log3,frag,k):

    seq = string.join([x for x in log3.pdbInfo[Log3Info.idxs[k]].sequence \
                       if x not in ['-',' ']],'')
    seq2 = string.join([x for x in log3.template.sequence \
                        if x not in ['-',' ']],'')
    diffs = string.join([x for x in log3.diffStrings[Log3Info.idxs[k]] \
                         if x not in ['-',' ']],'')

    a = pairwise2.align.localxs(frag,seq,-1,-0.5)[0]
    a2 = pairwise2.align.localxs(frag,seq2,-1,-0.5)[0]

    ss = string.join([cdict2[x] for x in log3.pdbInfo[Log3Info.idxs[k]].ss \
                      if x not in ['-']],'')
    ss2 = string.join([cdict2[x] for x in log3.template.ss \
                       if x not in ['-']],'')


    ssf1 = ss[a[-2]:a[-1]]
    ssf2 = ss2[a2[-2]:a2[-1]]

    minSS = min(len(ssf1),len(ssf2))

    if minSS == 0:
        return 0

    else:
        return float(len([x for x in xrange(min(len(ssf1),len(ssf2))) \
                          if (ssf1[x] == 'H' and ssf2[x] == 'E') or\
                          ssf1[x] == 'E' and ssf2[x] == 'H']))/minSS

if __name__ == '__main__':

    f = open('fastas.txt').read().splitlines()

    of = open('TableS1.csv','w')

    of.write('PDB1, PDB2, Frag_seq\n')

    log3Dict,IDdict = get_ldict(open('updated_pdb_ids.csv').read().splitlines())

    j = 0

    x = []
    y = []
    l = []

    x2 = []
    y2 = []
    l2 = []
    lengths = []

    j = 0

    for i in f:

        ID = string.split(i,'.')[0]

        if not log3Dict.has_key(ID):
            continue

        log3 = log3Dict[ID]

        whole_seq = open(RDIR+ID+'.seq').read().splitlines()[1]
        whole_pred = open(RDIR+ID+'.seq.jnetpred').read().splitlines()[1]

        frag_seq = open(CDIR+i).read().splitlines()[1]
        frag_pred = convert(open(CDIR+ID+'_download.log').read().splitlines())

        #Output information about fold switch pairs for Table S1
        of.write('%s, %s, %s\n' %(ID, IDdict[ID], frag_seq))
        a = pairwise2.align.localxs(frag_seq,whole_seq,-1,-0.5)[0]

        #This eliminates 2gedB, which is not a true fold switcher
        if len(whole_pred) != len(a[1]):
            continue

        wp = whole_pred[a[-2]:a[-1]]

        minSS = min(SS(frag_pred),SS(whole_pred[a[-2]:a[-1]]))

        if minSS == 0:
            nAB = 0

        if minSS == 0:
            nAB = 0
        else:
            nAB = len([i for i in xrange(len(frag_pred)) if \
                       (frag_pred[i] == 'H' and wp[i] == 'E') or \
                       (frag_pred[i] == 'E' and wp[i] == 'H')])/float(len(wp))

        Log3Info = LOG3(LDIR+log3)

        if Log3Info.idxs.has_key(ID):
            k = ID

        else:
            k = IDdict[ID]

        Log3Info.maxStrings[Log3Info.idxs[k]]
            

        nalphaBeta = get_maxString(Log3Info,frag_seq,k)

        if minSS/max(len(frag_pred),20.0) <= 0.35:
            continue

        x.append(nalphaBeta)
        y.append(nAB)
        l.append(minSS/max(len(frag_pred),20.0))

        if nalphaBeta >=0.18:
            l2.append(minSS/max(len(frag_pred),20.0))
            y2.append(nAB)
            lengths.append(len(frag_pred))
            j += 1

        #Print predicted Fold Switchers
        #if x[-1] >= 0.2:
            #print ID, l[-1],y[-1],x[-1], frag_pred, wp

    plt.hold(False)

    #Information about fold switchers
    #added after 2018.  First 3 are KaiB variants
    #with no structures.  Ratios were determined
    #using JPred and DSSP
    #x and x2: experimentally observed (or estimated for KaiB) SS
    #l and l2: SS:coil ratio
    #y and y2: %predicted helix <-> strand discrepancies
    #lengths: length of fragment used for FSR prediction
    #4KSO
    x.append(0.31)
    y.append(0.257)
    l.append(0.37)
    y2.append(0.257)
    x2.append(0.31)
    l2.append(0.37)
    lengths.append(35)
    #1WWJ
    x.append(0.31)
    y.append(0.257)
    l.append(0.37)
    y2.append(0.257)
    x2.append(0.31)
    l2.append(0.37)
    lengths.append(35)
    #1R5P
    x.append(0.31)
    y.append(0.23)
    l.append(0.37)
    y2.append(0.23)
    x2.append(0.31)
    l2.append(0.37)
    lengths.append(35)
    #ORF9b
    x.append(0.5)
    y.append(0.25)
    l.append(0.5)
    y2.append(0.25)
    x2.append(0.5)
    l2.append(0.5)
    lengths.append(20)
    #BAX
    x.append(0.81)
    y.append(0.76)
    l.append(0.5)
    y2.append(0.76)
    x2.append(0.81)
    l2.append(0.5)
    lengths.append(22)
    #minE
    x.append(0.45)
    y.append(0.35)
    l.append(0.5)
    y2.append(0.35)
    x2.append(0.45)
    l2.append(0.5)
    lengths.append(20)

    of.write('4KSO, , LAKVLPLPVRRIIGDLSDREKVLIGLDLLYGELQD\n')
    of.write('1WWJ, , LAKILPPPVRKIIGDLSDREKVLIGLDLLYDEIRE\n')
    of.write('1R5P, , LSKILPPPVRKIIGDLSDRERVLIGLDLLYEELTE\n')
    of.write('1F16, , WDGLLSYFGTPTWQTVTIFVAG\n')
    of.write('2KXOA, 3RJ9C, VARDRLQIIIAQERAQEGQT\n')
    of.write('6Z4UA, 7KDTB, KTLNSLEDKAFQLTPIAVQMTK\n')

    #End add

    #Close Table S1
    of.close()
    
    plt.hold(False)
    npoints = []
    ratio = []
    cutoff = []
    dTPdFN = []
    TPFN = []
    
    for i in np.linspace(0,0.8,num=81):

        npoints.append(len([n for n in x if n >= i]))
        nTP = float(len([y[n] for n in xrange(len(x)) if y[n] > 0 \
                         and x[n] >=i]))
        nFN = len([n for n in x if n >= i])-nTP
        
        ratio.append(float(len([y[n] for n in xrange(len(x)) if y[n] > 0 \
                            and x[n] >=i]))**2/len([n for n in x if n >= i]))
        cutoff.append(i)


    #Generate Figure S1
    plt.plot(cutoff,ratio,'k^-',markersize=12)

    m = max(ratio)
    midx = [n for n in xrange(len(ratio)) if ratio[n]==m][0]
    mx = cutoff[midx]

    plt.hold(True)
    plt.plot(mx,m,'r^',markersize=12)
    
    plt.ylabel('Number of True Positives * %True Positives')
    plt.xlabel('%Observed Helix <-> Strand Discrepancies')

    plt.savefig('fig_S1.png')

    #Figure S1 finished

    plt.hold(False)


    #Generate Figure 1
    plt.figure(figsize=(5,5), dpi=300)
    
    for i in xrange(len(x)):
        if x[i] >=0.18:
            plt.plot(x[i],y[i],'ro')
            plt.hold(True)

        else:
            plt.plot(x[i],y[i],'ko')
            plt.hold(True)
    

    plt.xlabel(r'Fraction $\alpha \longleftrightarrow \beta$ discrepancies (Observed)')
    plt.ylabel(r'Fraction $\alpha \longleftrightarrow \beta$ discrepancies (Predicted)')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    plt.hold(True)
    plt.plot([0,1.37],[slope*0+intercept,slope*1.37+intercept],'k--')
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.savefig('fig1.png',dpi=600)

    #Figure 1 finished

    plt.hold(False)

    #Output 

    of = open('fs_xy.pik','w')

    #Append amount of predicted secondary structure (l2) and
    #% helix <-> strand discrepancies (y2) to pickled file for
    #use in 
    cPickle.dump(l2,of)
    cPickle.dump(y2,of)

    of.close()

    #Print correlation

    print 'Correlation coefficient: %f' %(np.corrcoef(x,y)[0][1])


 
        

    
