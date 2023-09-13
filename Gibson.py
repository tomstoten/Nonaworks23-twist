from Bio.Seq import Seq
from difflib import SequenceMatcher

gene = Seq(
    "TCCCTGGGCTCTTTTAGTGGACGGAGACCCAGCTGTCAGTTTGTTGTAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACCACTTACTCCATGGAGGGATTAGATGATTCACGGTAGGCTTGGGCAG"
)

oligo1 = Seq(
    "TCCCTGGGCTCTTTTAGTGGACGGAGACCCAGCTGTCAGTTTGTTGTAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
)

oligo2 = Seq(
    "CTGCCCAAGCCTACCGTGAATCATCTAATCCCTCCATGGAGTAAGTGGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
)

def findLongestMatch(seq1, seq2):
# find the complement sequence of 1 oligo.
  comp = seq2.reverse_complement()
# find longest match from both oligos
  match = SequenceMatcher(None, seq1, comp).find_longest_match()
  
  return match, comp

#combine two oligos at their longest match locations
def combineFromLongestMatch(oligo1, oligo2):
  match, comp = findLongestMatch(oligo1, oligo2)
  output = oligo1[:match.a] + comp[match.b:]

  return output

bestGibsonOutput = combineFromLongestMatch(oligo1, oligo2)
print("=========== Gibson Assembly ===========")
print("Target gene", gene)
print("Best match:", bestGibsonOutput)
print("Match == Gene? -->", gene == bestGibsonOutput)