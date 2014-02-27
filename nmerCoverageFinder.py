__author__ = 'mitch'


import sys, string, gzip


transtab = string.maketrans('atcgATCG', 'tagcTAGC')

class contig:
    def __init__(self, name, forseq):
        self.name = name
        totalcov = 0
        covlist = []
        for i in range(0, len(forseq) - nmersize):
            nmer = forseq[i:i + nmersize]
            if nmer[nmersize / 2] == 'a' or nmer[nmersize / 2] == 'c':
                nmer = nmer[::-1]
                nmer = nmer.translate(transtab)
            try:
                totalcov += nmerdict[nmer]
                covlist.append(nmerdict[nmer])
            except:
                covlist.append(0)
        covlist.sort()
        self.coverage = covlist[len(covlist)/2]
        self.to = []
        self.fr = []
        self.forseq = forseq
        revseq = forseq[::-1]
        revseq = revseq.translate(transtab)
        self.revseq = revseq
        self.blast = []
        self.blast2 = set()
        self.placed = False
        self.cand = False
    def length(self):
        return len(self.forseq)

    def repeat(self):
        if len(self.to) > 1 or len(self.fr) > 1:
            return True

    def revseq(self):
        temp = self.seq[::-1]
        return temp.translate(transtab)

    def gccont(self):
        return float(self.forseq.count('a') + self.forseq.count('t')) / len(self.forseq) * 100


def getnmers(reads, prefix):
    global nmersize, nmercut
    nucl = set('atcg')
    nmerdict = {}
    if reads[-3:] == '.gz':
        readfile = gzip.open(reads)
    else:
        readfile = open(reads)
    getseq = False
    getfaseq = False
    for line in readfile:
        if line.startswith('@'):
            getseq = True
        elif line.startswith('>'):
            if getfaseq:
                for i in range(0, len(seq) - nmersize):
                    nmer = seq[i:i+nmersize]
                    if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                        nmer = nmer[::-1]
                        nmer = nmer.translate(transtab)
                    if nmer in nmerdict:
                        nmerdict[nmer] += 1
                    else:
                        if set(nmer).issubset(nucl):
                            nmerdict[nmer] = 1
            getfaseq = True
            seq = ''
        elif getfaseq:
            seq += line.rstrip().lower()
        elif getseq:
            getseq = False
            seq = line.rstrip().lower()
            for i in range(0, len(seq) - nmersize + 1):
                nmer = seq[i:i+nmersize]
                if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                    nmer = nmer[::-1]
                    nmer = nmer.translate(transtab)
                if nmer in nmerdict:
                    nmerdict[nmer] += 1
                else:
                    if set(nmer).issubset(nucl):
                        nmerdict[nmer] = 1
    if getfaseq:
        for i in range(0, len(seq) - nmersize + 1):
            nmer = seq[i:i+nmersize]
            if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                nmer = nmer[::-1]
                nmer = nmer.translate(transtab)
            if nmer in nmerdict:
                nmerdict[nmer] += 1
            else:
                if set(nmer).issubset(nucl):
                    nmerdict[nmer] = 1
    out = open(prefix + '_nmers', 'w')
    for i in nmerdict:
        out.write(i + '\t' + str(nmerdict[i]) + '\n')
    out.close()

def readnmerfile(prefix):
    global nmerdict, nmercut
    nmers = open(prefix + '_nmers')
    nmerdict = {}
    for line in nmers:
        nmer, freq = line.split()
        freq = int(freq)
        if freq >= nmercut:
            nmerdict[nmer] = freq

try:
     nmersize = int(sys.argv[4])
except:
    print '''nmerCoverageFinder.py
USE: Find the median nmer coverage of contigs
USAGE: python nmerCoverageFinder.py contigfile.fa readfile.fa output_prefix nmer_size
outputs nmer counts to output_prefix_nmer and coverage csv to output_prefix_coverage.csv'''
nmercut = 0
getnmers(sys.argv[2], sys.argv[3])
readnmerfile(sys.argv[3])
contigfile = open(sys.argv[1])
first = True
contiglist = []
for line in contigfile:
    if line.startswith('>'):
        if first:
            first = False
        else:
            aninstance = contig(name, seq.lower())
            contiglist.append(aninstance)
        name = line.rstrip()[1:]
        seq = ''
    else:
        seq += line.rstrip()
aninstance = contig(name, seq.lower())
contiglist.append(aninstance)
out = open(sys.argv[3] + '_coverage.csv', 'w')
for i in contiglist:
    out.write(i.name + '\t' + str(i.gccont()) + '\t' + str(i.coverage) + '\n')


