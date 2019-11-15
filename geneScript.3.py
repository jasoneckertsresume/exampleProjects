# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

author: Jason Eckert

01001010
01100001
01110011
01101111
01101110

Classes:
    genSeq - Object to handle genetic sequences. Can be defined by passing either a RNA or 
        Attributes:
                seq = original sequence
                seqDna = DNA seqquence with the same coding values as seq (if seq isDna then seqDna = seq)
                seqRna = RNA seqquence with the same coding values as seq (if seq isDna then seqRna = seq)
                opens = one dimensional list with the locations of all valid start codons (validity defined as the start and stop codon existing within the same reading frame)
                closes = one dimensional list with the locations of all valid stop codons (validity defined as the start and stop codon existing within the same reading frame)
                orfs = list containing DNA nucleotide sequences existing between each pair of valid start and stop codons 
                aaSeqs = list containing each element of orfs translated into 1-letter amino acid codes
        Functions:
            __init__ = constructor
            orfMatch = sets the values of opens and closes to contain only stop/start loci pairs that are in the same reading frame
            getOrfs = returns a list with the DNA nucleotide sequences between valid stop/start loci pairs
            translate = translates the DNA sequences contained within orfs into amino acid codes
            print = function to print a genSeq object

Common Parameters:
    seq = any nucleotide sequence DNA or RNA
    aaSeq = any sequence of amino acids
    anySeq = any molecular sequence
    seqDna = any DNA sequence
    seqRna = any RNA sequence
    
Functions:
    threeFrames = returns a 3 element list with all three reading frames of the sequence
    find_substring = returns a list with the indexes of the first character of all occurances of a substring within a string
    openSearch = returns the index of the first base of every "aug" codon in an RNA sequence
    closeSearch = returns the index of the first base of every STOP codon ("uaa", "uga", "uag")
    orfMatchStandalone = returns a list of the indices of the first character of start/stop codon pairs in the same reading frame
    orfSearch = returns an object with indices of valid start/stop codon pairs
    isDna = evaluates a nucleotide sequence and returns 1 if the sequence is DNA and 0 is the sequence is RNA (checks for presence of uracil)
    baseTranslate = accepts a nucleotide sequence (DNA or RNA) and converts it to its complement of the opposite type
    tDimers = returns a list of the indices of thiamine dimers (allows for overlapping dimers)
    mutateDna = accepts a DNA sequence and performs probablistic mutations at each locus
    mutateRna = accepts a RNA sequence and performs probabilistic mutations at each locus
    molecWeight = accepts an Amino Acid sequence and returns the total molecular weight of the sequence
    aveMolecWeight = accepts an Amino Acid sequence and returns the average molecular weight of an amino acid in that sequence
    getFreq = gets the total counts and relative frequencies of the target molecule(s)
    getNucFreq = gets the total counts and relative frequencies of nucleotides in a nucleotide sequence
    getAaTypes = gets the type of the amino acid
    getMutLocs = compares each base of 2 strings to search for mutations
"""

import random
import pandas as pd

class genSeq :
    seq = ""
    seqDna = ""
    seqRna = ""
    opens = []
    closes = []
    orfs = []
    aaSeqs = []
    def __init__ (self, seq):
        self.seq = seq.lower()
        if isDna(self.seq) == 0:
            self.seqRna = self.seq
            self.seqDna = baseTranslate(self.seq)
        if isDna(self.seq) == 1:
            self.seqRna = baseTranslate(self.seq)
            self.seqDna = self.seq
        self.opens = openSearch(self.seqDna)
        self.closes = closeSearch(self.seqDna)
        self.orfMatch()
        self.getOrfs()
        self.translate()
    def orfMatch(self):
        openLocs = []
        closeLocs = []
        checkList = []
        for openLocus in self.opens:
            for closeLocus in self.closes:
                if openLocus < closeLocus and (closeLocus - openLocus)%3 == 0:
                    openLocs.append(openLocus)
                    closeLocs.append(closeLocus)
                    checkList.append(0)
        if sum(checkList) == 0:
           self.opens = openLocs
           self.closes = closeLocs
        else:
           return "ERROR"
    def getOrfs(self):
        i = 0
        for openLoc in self.opens:
            self.orfs.append(self.seqDna[openLoc:(self.closes[i]+3)]) 
            i = i+1
#DEAD FUNCTION   
#   def getCodons(self):
#        for orf in self.orfs:
#            i = 0
#            while i < len(orf):
#                self.codons.append(orf[i:i+3])
#                i = i +3
#DEAD FUNCTION
#    def translate(self):
#        for orf in self.orfs:
#            i = 0
#            j = 0
#            while i < len(orf):
#                codons = []
#                codons.append(orf[i:i+3])
#                i = i +3
#                for codon in codons:
#                    self.aaSeqs.append(codonTable.get(codon.upper()))
    def translate(self):
        for orf in self.orfs:
            i = 0
            aaString = ""
            while i < len(orf):
                codons = []
                codons.append(orf[i:i+3])
                i = i +3
                for codon in codons:
                    aaString = aaString + str(codonTable.get(codon.upper()))
            self.aaSeqs.append(aaString)
                
    def print(self):
        print(self.opens)
        print(self.closes)
       # print(self.codons)
        print(self.orfs)
        print(self.aaSeqs)
    

def threeFrames(seq):
    seq1 = seq
    seq2 = seq[1:]
    seq3 = seq[2:]
    return seq1, seq2, seq3

def find_substring(substring, string):
    indices = []
    index = -1  
    while True:
        index = string.find(substring, index + 1)
        if index == -1:  
            break
        indices.append(index)
    return indices

def openSearch(seq):
    frameOpens = find_substring(baseTranslate("aug"), seq)
    return frameOpens

def closeSearch(seq):
    frameClose1 = find_substring(baseTranslate("uaa"), seq)
    frameClose2 = find_substring(baseTranslate("uag"), seq)
    frameClose3 = find_substring(baseTranslate("uga"), seq)
    frameCloses = frameClose1 + frameClose2 + frameClose3
    return frameCloses

def orfMatchStandAlone(frameOpens, frameCloses):
    opens = []
    closes = []
    checkList = []
    for openLocus in frameOpens:
        for closeLocus in frameCloses:
            if (closeLocus - openLocus)%3 == 0:
                opens.append(openLocus)
                closes.append(closeLocus)
                checkList.append(0)
    if sum(checkList) == 0:
       return opens, closes
    else:
       return "ERROR"


def orfSearch(seq):
    return orfMatchStandAlone(openSearch(seq), closeSearch(seq))

#Returns 0 if the sequence in RNA and 1 if the sequence is DNA
#Based on detection of Uracil, could be incorrect if RNA sequence
#Contains no Uracil bases(Unlikely)
def isDna(seq):
    if (seq.find("u") != -1 or seq.find("U") != -1):
        return 0
    else:
        return 1

#Detects whether a sequence is an RNA or DNA sequence and
#Translates it to the other nucleotypye 
#Translates to retain coding functionality
def baseTranslate(seq):
    dnaFlag = isDna(seq)
    seq = seq.lower()
    if dnaFlag == 0:
        seq = seq.replace("a", "t")
        seq = seq.replace("u", "a")
        seq = seq.replace("c", "1")
        seq = seq.replace("g", "c")
        seq = seq.replace("1", "g")
        return seq
    if dnaFlag == 1:
        seq = seq.replace("a", "u")
        seq = seq.replace("t", "a")
        seq = seq.replace("c", "1")
        seq = seq.replace("g", "c")
        seq = seq.replace("1", "g")
        return seq
    else:
        return "ERROR"

def tDimers(seq):
    dimerList = []
    dimerList = find_substring("tt", seq)
    return dimerList

#Probabilities to be entered as values between 0 and 1
#Probabilities refer to the likelihood that any individual base/dimer becomes mutant
#Max Precision for *Probs's e-7
def mutateDna(seq, tDimerProb, ranSwapProb, ranDelProb, ranAddProb):
    tDimerLocs = tDimers(seq)
    tDimerProb = 10000000*tDimerProb
    ranSwapProb = 10000000*ranSwapProb
    ranDelProb = 10000000*ranDelProb
    ranAddProb = 10000000*ranAddProb
    for tDimer in tDimerLocs:
        if random.randint(10000000) <= tDimerProb:
            seq = seq[:tDimer] + seq[tDimer+1:]
    for i in seq.len():
        if random.randint(10000000) <= ranDelProb:
            seq = seq[:i] + seq[i+1:]
    for j in seq.len:
        if random.ranint(10000000) <= ranSwapProb:
            baseChoice = random.ranint(4) 
            baseList = ["a","t","c","g"]
            seq = seq[:i-1] + baseList[baseChoice] + seq[i+1:]
    for k in seq.len():
        if random.ranint(10000000) <= ranAddProb:
            baseChoice = random.ranint(4)
            baseList = ["a", "t", "c", "g"]
            seq = seq[:i] + baseList[baseChoice] + seq[i+1:]
    return seq

def mutateRna(seq, tDimerProb, ranSwapProb, ranDelProb, ranAddProb):
    tDimerLocs = tDimers(seq)
    tDimerProb = 10000000*tDimerProb
    ranSwapProb = 10000000*ranSwapProb
    ranDelProb = 10000000*ranDelProb
    ranAddProb = 10000000*ranAddProb
    for tDimer in tDimerLocs:
        if random.randint(10000000) <= tDimerProb:
            seq = seq[:tDimer] + seq[tDimer+1:]
    for i in seq.len():
        if random.randint(10000000) <= ranDelProb:
            seq = seq[:i] + seq[i+1:]
    for j in seq.len:
        if random.ranint(10000000) <= ranSwapProb:
            baseChoice = random.ranint(4) 
            baseList = ["a","u","c","g"]
            seq = seq[:i-1] + baseList[baseChoice] + seq[i+1:]
    for k in seq.len():
        if random.ranint(10000000) <= ranAddProb:
            baseChoice = random.ranint(4)
            baseList = ["a", "u", "c", "g"]
            seq = seq[:i] + baseList[baseChoice] + seq[i+1:]
    return seq

def molecWeight(aaSeq):
    mw = 0
    for aa in range(len(aaSeq)):
        if aaSeq[aa] in list(aaTabGen.index):
            mw = mw + aaTabGen['MW (weight)'].loc[aaSeq[aa]]
    return mw

def aveMolecWeight(aaSeq):
    amw = molecWeight/len(aaSeq)
    return amw

def getFreq(molecs, anySeq):
    mols = []
    counts = []
    freqs = []
    molecs = list(molecs)
    for molec in molecs:
        count = len(find_substring(molec,anySeq))
        mols.append(molec)
        counts.append(count)
        freqs.append(count/len(anySeq))
    table = {"mol": mols, "count": counts, "freq": freqs}
    table = pd.DataFrame(table)
    return table

def getNucFreq (seq):
    if isDna(seq) == 0:
        return getFreq(rnaMolecs, seq)
    else:
        return getFreq(dnaMolecs, seq)
    
def getAaFreq (aaSeq):
    return getFreq(aaMolecs, aaSeq)

def getAaTypes (aaSeq):
    typeSeq = ""
    for aa in range(len(aaSeq)):
        k = 0
        for i in aaTabTypes["Amino Acids Included"]:
            if aaSeq[aa] in i:
                typeSeq = typeSeq + aaTypes[k]
            k = k + 1
    return typeSeq

def getMutLocs (seq1, seq2):
    locs = []
    for i in range(max(len(seq1), len(seq2))):
        if seq1[i] != seq2[i]:
            locs.append(i)
    return locs

codonTable = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 

aaTables = pd.read_html("https://en.wikipedia.org/wiki/Amino_acid")
aaTabGen = aaTables[2]
aaTabGen = aaTabGen.set_index('1-letter[130]')
aaTabSpec = aaTables[3]
aaTabAmbig = aaTables[4]
aaTabAmbig["Amino Acids Included"] = aaTabAmbig["Amino Acids Included"].str.split(", ")
aaTabTypes = aaTabAmbig[4:10]

secStructs = pd.read_html("http://www.bmrb.wisc.edu/referenc/choufas.shtml")
secStructs = secStructs[2]
secStructs = secStructs[0, 3, 5:]

dnaMolecs = ["a", "t", "g", "c"]
rnaMolecs = ["a", "u", "g", "c"]
aaMolecs = list(aaTabGen.index)
aaTypes = list(aaTabTypes["1-letter"])
alphas = ["A", "L", "M", "F", "E", "Q", "H", "K", "R"]
betas = ["Y", "W", "I", "V", "T", "C"]
breakers = ["G", "P", "S", "D", "N", "S"]


testSeq = "UGUCAUUACUACAUUGCGUACGCUUAGGGUAUCUAGACCAAUAUGCAUUUGCUCGAGUUUUCGAAAAUAAUCCGGUACAGCAAAACGCCCAAGCUCGUUCUGGACGGCACGUGACAUCCAGUCAAAGCUGCCCUUUAAGCCUUCCGGACGGCCGGCGACCAACGCUGUGAGAGUUACCGGUUAACCAUUUUUCGACAGUGAAUUGGGAGCAGAAAGAAGGAACAACUAAGCGAAAUGGCCCGGGGAAAAUUAACCUCAAAACGAUGCAAUUAGGCGAAGCGAGCUCAAUAGUCAGGAUAGUCACCUCCAUAUCAAUCCCGUAGCUUCACCUUCAGGACGAUUCAGCAUCAGGGAACGACCGAGAUUAUCAGGGAUUCAGGUGAUGAAUAGCCCUGUAAAUCUGCCGAGACAACGCUAUGCCCUUUUAUGCAACCAGGUCCAAGUUUUACCUCUAGGCCUCGGGUCAGGCAACCUGGGCCUUUCAUCGACGAUACCAAGGUAUGGGGAAAACUAAUUAUCCUCUGUCUCACAUUGGCCUGACUGGGUAGAUAAAGGGAAGUAGCAACUGUAAUUCUGUCUGGGAGAAGACCAACAACUACACACAAUGAAAGCUUAACGUCACUGAGCAGGCCUAAAUGUUAUUGCUCAGCGACAUUGAUGCCUUGAUCAAUUGCGCACACUGUCUUAAAAAUGAUUAAUUCUGUGUCAGAGGACCCGUUAACGCCCGACGCUAAAUCCGCACCAAUCACAACCCUUAGCUAUCUGACAGUGGACAUAUUCUAGGAGAUUUAGCGGCUCUUACACUGUUCACUGAUUAUCCUUGCCUAAACGGAUGGGAUGUUGACCAUCAGCGCUAAACGAUCAAGCGGGCCACAGGGCGAUCGCCAAUUAAGGUGUAGGUGACGGUCUAUUAUUGCGACGGUCUCGCACGUGCUUUAACGGCCGACAGGCUUCGUCUGGCGUCAGGCUCUUGAUCCAGCUGCUAGAAAAUUAAUGGCUG"
testSeq = testSeq.lower()
testSeq1 = "AUGC"
testSeq2 = "ACGT"
ts3 = "aaaaaaaau"
ts4 = "aaaaaaaaa"
a = threeFrames(testSeq)
b = openSearch(testSeq)
c = closeSearch(testSeq)
d = orfMatchStandAlone(b,c)
e = orfSearch(testSeq)
f = genSeq(testSeq)
g = tDimers(baseTranslate(testSeq))
#f.orfMatch()
#f.getOrfs()
#f.translate()
print(isDna(testSeq1), isDna(testSeq2))
print(baseTranslate(testSeq))
#f.getCodons()
#f.print()
print(g)

#print(orfSearch(testSeq))