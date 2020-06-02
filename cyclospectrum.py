#!/usr/bin/env python3
# Max Genetti (mgenetti)

'''
This program generates the theoretical spectrum of a cyclic peptide.

Main expects a file formatted for rosalind problem 25 handled by a general file
reader class FileReader.

Given: An amino acid string Peptide.

Return: Cyclospectrum of the Peptide.
'''

import sys
import numpy as np



class FileReader :
   '''
   Note:
      Edited from class FastaReader provided by David Bernick
   Args:
      fname: a file with a string of amino acids
      
   Returns:
      peptide: the string of amino acids
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

      peptide = ''
      with self.doOpen() as fileH:
         peptide = fileH.readline()
      return peptide


class CyclicPeptide:
   '''
   Note:
      This class generates the spectrum of a cyclic peptide
   Args:
      peptide: a peptide sequence
   Returns:
      possible dna sequences coding this peptide
   '''

   aa2mw = {
      'A': 71,  'G': 57,  'M': 131, 'S': 87, 'C': 103,
      'H': 137, 'N': 114, 'T': 101, 'D': 115, 'I': 113,
      'P': 97, 'V': 99, 'E': 129, 'K': 128, 'Q': 128,
      'W': 186, 'F': 147, 'L': 113, 'R': 156, 'Y': 163
      }

   def __init__(self, peptide) :
      '''
      Stores sequence and peptide
      '''
      self.peptide = peptide
      

   def printSpectrum(self) :
      '''
      Note:
         This method prints the spectrum of a cyclic peptide.
         In order of smallest to largest spectra.
      '''


      pep = self.peptide + self.peptide[:-2]
      #linear representation of clyclic relationship
      profile = [0, 0] #initialize profile

      for aa in self.peptide: #calculate molecular weight of cyclic peptide
         profile[1] += CyclicPeptide.aa2mw[aa]

      for i in range(len(self.peptide)):
         #calculate molecular weights of all fragments
         mw = 0
         for j in range(len(self.peptide)-1):
            mw += CyclicPeptide.aa2mw[pep[i+j]]
            profile.append(mw)

      profile.sort() #sort the string

      print(' '.join(str(x) for x in profile))


def main() :

   # read file and store as matrices
   fileReader = FileReader('')
   peptide = fileReader.readFile()

   # initialize object and print
   pep = CyclicPeptide(peptide)
   pep.printSpectrum()


if __name__ == "__main__":
   main()
