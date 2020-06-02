#!/usr/bin/env python3
# Max Genetti (mgenetti)

'''
Find substrings of a genome encoding a given amino acid sequence.

Main expects a file formatted for rosalind problem 24 handled by a general file
reader class FileReader.

Given: A DNA string Text and an amino acid string Peptide.

Return: All substrings of Text encoding Peptide (if any such substrings exist).
'''

import sys
import numpy as np



class FileReader :
   '''
   Note:
      Edited from class FastaReader provided by David Bernick
   Args:
      fname: a file formatted with a string in each line
      
   Returns:

   '''
   
   def __init__ (self, fname='') :
      '''constructor: saves attribute fname'''
      self.fname = fname
      self.lines = []

   def doOpen (self) :
      '''opens the file'''
      if self.fname is '':
         return sys.stdin
      else:
         return open(self.fname)
 
   def readFile (self) :
      '''parses the file'''
      sequence = ''
      peptide = ''
      with self.doOpen() as fileH:
         sequence = fileH.readline()
         peptide = fileH.readline()
      return sequence, peptide


class PeptideSearch:
   '''
   Note:
      This class finds the sequences in a DNA strand coding for a given peptide
   Args:
      sequence: a string of A DNA sequence and
      peptide: a string of a peptide sequence
   Returns:
      The DNA subsequences coding for the amino acid
   '''
   dnaCodonTable = {
      # DNA codon table
      # T
      'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', # TxT
      'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', # TxC
      'TTA': 'L', 'TCA': 'S', 'TAA': '-', 'TGA': '-', # TxA
      'TTG': 'L', 'TCG': 'S', 'TAG': '-', 'TGG': 'W', # TxG
      # C
      'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', # CxT
      'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', # CxC
      'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', # CxA
      'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', # CxG
      # A
      'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', # AxT
      'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', # AxC
      'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', # AxA
      'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', # AxG
      # G
      'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', # GxT
      'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', # GxC
      'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', # GxA
      'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
      }
      
   dnaCodonTableCompliment = {
      # DNA codon table of Reverse Compliments
      
      'AAA': 'F', 'AGA': 'S', 'ATA': 'Y', 'ACA': 'C',
      'GAA': 'F', 'GGA': 'S', 'GTA': 'Y', 'GCA': 'C',
      'TAA': 'L', 'TGA': 'S', 'TTA': '-', 'TCA': '-',
      'CAA': 'L', 'CGA': 'S', 'CTA': '-', 'CCA': 'W',
      
      'AAG': 'L', 'AGG': 'P', 'ATG': 'H', 'ACG': 'R',
      'GAG': 'L', 'GGG': 'P', 'GTG': 'H', 'GCG': 'R',
      'TAG': 'L', 'TGG': 'P', 'TTG': 'Q', 'TCG': 'R',
      'CAG': 'L', 'CGG': 'P', 'CTG': 'Q', 'CCG': 'R',
      
      'AAT': 'I', 'AGT': 'T', 'ATT': 'N', 'ACT': 'S',
      'GAT': 'I', 'GGT': 'T', 'GTT': 'N', 'GCT': 'S',
      'TAT': 'I', 'TGT': 'T', 'TTT': 'K', 'TCT': 'R',
      'CAT': 'M', 'CGT': 'T', 'CTT': 'K', 'CCT': 'R',
      
      'AAC': 'V', 'AGC': 'A', 'ATC': 'D', 'ACC': 'G',
      'GAC': 'V', 'GGC': 'A', 'GTC': 'D', 'GCC': 'G',
      'TAC': 'V', 'TGC': 'A', 'TTC': 'E', 'TCC': 'G',
      'CAC': 'V', 'CGC': 'A', 'CTC': 'E', 'CCC': 'G'
      }

   def __init__(self, sequence, peptide) :
      '''
      Stores sequence and peptide
      '''
      self.sequence = sequence
      self.peptide = peptide


   def parseSequence(self):
      '''
      Note:
         This method returns the sequences coding for the peptide in both 
         the forward and reverse strands.
      '''

      p = len(self.peptide)

      for i in range(len(self.sequence)-(p*3-1)): #iterate through sequence
         
         seq = ''
         for j in range(0, p): #Check Forward Strand
            if(PeptideSearch.dnaCodonTable[self.sequence[i+(3*j):i+(3*j)+3]] == self.peptide[j]):
               #checks forward strand
               seq += self.sequence[i+(3*j):i+(3*j)+3]
               if(j == p-1):
                  print(seq)
            else:
               break

         seq = ''
         for j in range(0, p): #Check Reverse Compliment
            if(PeptideSearch.dnaCodonTableCompliment[self.sequence[i+(3*j):i+(3*j)+3]] == self.peptide[-(j+1)]):
               #checks peptide in reverse
               seq += self.sequence[i+(3*j):i+(3*j)+3]
               if(j == p-1):
                  print(seq)
            else:
               break

def main() :

   # read file and store as matrices
   fileReader = FileReader('')
   sequence, peptide = fileReader.readFile()

   # initialize hmm object and print the probability
   seqs = PeptideSearch(sequence, peptide)
   seqs.parseSequence()


if __name__ == "__main__":
   main()
