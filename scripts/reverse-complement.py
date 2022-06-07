#!/usr/bin/python
import sys
def complement(seq):
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    return complement[seq]
def reverse_complement(seq):
    # seq[::-1]  is the string seq, in reverse order. 
    reverseComplement = ""
    for residue in seq[::-1]           :
        reverseComplement += complement(residue)
    print(reverseComplement)
print("Usage: reverse-complement.py <RNA sequence>")
print ("Reverse Complement of : ",sys.argv[1], " is : ")
#reverse_complement("UCGGGCCCC")
reverse_complement(sys.argv[1])
