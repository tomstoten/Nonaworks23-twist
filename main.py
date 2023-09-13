from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Restriction.Restriction import RestrictionBatch
import re
import Gibson
from Bio.SeqUtils import MeltingTemp

gene = Seq(
    "TCCCTGGGCTCTTTTAGTGGACGGAGACCCAGCTGTCAGTTTGTTGTAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACCACTTACTCCATGGAGGGATTAGATGATTCACGGTAGGCTTGGGCAG"
)

oligo1 = Seq(
    "TCCCTGGGCTCTTTTAGTGGACGGAGACCCAGCTGTCAGTTTGTTGTAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
)

oligo2 = Seq(
    "CTGCCCAAGCCTACCGTGAATCATCTAATCCCTCCATGGAGTAAGTGGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
)

nucleo_reference = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "U",
    "R": "GA",
    "Y": "CT",
    "K": "GT",
    "M": "AC",
    "S": "GC",
    "W": "AT",
    "B": "GTC",
    "D": "GAT",
    "H": "ACT",
    "V": "GCA",
    "N": "AGCT"
}


# Formula from http://biotools.nubic.northwestern.edu/OligoCalc.html
def Tm(seq):
  wA = seq.count("A")
  xT = seq.count("T")
  yG = seq.count("G")
  zC = seq.count("C")
  if len(seq) < 14:
    tm = (wA+xT) * 2 + (yG+zC) * 4
  else:
    tm = 64.9 + 41 * (yG+zC - 16.4)/(wA+xT+yG+zC)
  return tm
  


def specifyOverhang(overhangDict, overhang, indices, seq, topOverhang):
  checkOverhang = False
  for c in overhang:
    if c not in 'ATCG':
      checkOverhang = True
      break
    #if nucleotide is not ATCG, we want to go find the index in oligo, and set that range as the new key for the overhang.
  if checkOverhang:
    length = len(overhang)
    # print(overhang)
    for i in indices:
      realOverhang = seq[i - 1:i - 1 + length]
      # if overhang sequence already in overhang, append new indeces and
      if realOverhang in overhangDict:
        overhangDict[realOverhang][i] = topOverhang
        # overhangDict[realOverhang]
      else:
        overhangDict[realOverhang] = {i: topOverhang}


def getIndexAndOverhang(seq):
  length = len(seq)
  a = Analysis(AllEnzymes, seq).with_sites()
  # print(a)
  overhangDict = {}
  #find cut sites, and overhangs for each enzyme.
  for k, v in dict(a).items():
    # print(str(k) + ", " + k.site + ": " + k.elucidate())
    elu = k.elucidate()
    site = elu.replace('^', '_')
    site_sep = site.split('_')
    # if _ before ^  top over hang, if ^_ = bottom overhang
    #overhangFlags False = bottom , True = top
    if len(site_sep) > 2:
      topOverhang = False
      if elu.index('^') > elu.index('_'):
        topOverhang = True
      overhang = site_sep[1]
      if overhang != '':
        specifyOverhang(overhangDict, overhang, v, seq, topOverhang)

  return overhangDict


def overhangsCompatible(str1, str2, f1, f2):
  # if overhang in one oligo == overhang in other oligo AND they are opposite overhang flags
  return str1 == str2 and f1 != f2


def getValidCombinations(oligo1, oligo2):
  #get dictionary of cut sites for each oligo and its reverse complement
  o1 = getIndexAndOverhang(oligo1)
  o2 = getIndexAndOverhang(oligo2)

  combos = {}

  for key1 in o1:
    if key1 in o2:
      for k1, b1 in o1[key1].items():
        for k2, b2 in o2[key1].items():
          if b1 != b2:
            ss1 = oligo1[:k1]
            ss2 = oligo2[k2:]
            avg_temp = (Tm(ss1) + Tm(ss2))/2
            combos[ss1 + ss2] = avg_temp

  return combos


seqs = getValidCombinations(oligo1, oligo2.reverse_complement())
sorted_seqs = sorted(seqs.items(), key=lambda x: x[1])
#print(sorted_seqs)
print("============ Restriction =============")
for s in sorted_seqs:
  print(s[0], s[1])