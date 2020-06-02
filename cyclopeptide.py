#!/usr/bin/env python3
# Max Genetti (mgenetti)

'''
This program can determine the possible sequences of a cyclopeptide given an
ideal spectrum (where all subsets of the full length sequence are included).

Main expects a file formatted for rosalind problem 26 handled by a general file
reader class FileReader.

Given: A spectrum corresonding to an cyclopeptide of unknown sequence

Return: All possible sequences of Peptide, if any such exist.
'''

import sys
from collections import defaultdict


class FileReader :
   '''
   Note:
      This class expects a string of ints consisting of the ideal spectrum of a
      cyclopeptide.
   Args:
      fname: a file formatted with numbers seperated by white space
   Returns:
      A list of the numbers
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
      spectrum = []
      with self.doOpen() as fileH:
         spectrum = list(int(x) for x in fileH.readline().split())
         
      return spectrum


class CycloPeptideSequencing:
   '''
      Note:
         This class takes in a spectrum of a cyclic peptide and finds all of
         the possible cyclopeptide mass sequences to generate the spectrum
      Args:
         spectrum: the ideal spectrum of the cycliceptide
      Returns:
         Prints the cyclopeptides found matching the spectrum
      ''' 

   def __init__(self, spectrum) :
      self.cyclopeptideSequencing(spectrum)


   def expand(self, peptides) :
      '''
         Note:
            This method expands the peptide list by adding each possible
            amino acid spectrum, some of which include multiple amino acids.
         Args:
            peptides: a list of the current possible cyclopeptides returned from
               cyclopeptidesequencing.
         Returns:
            newPeptides: The updated peptide list with each amino acid spectrum
               added to each peptide in pepList.
      '''
      aaList = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
      #all amino acid spectrum
      newPeptides = []
      for peptide in peptides:
      #iterate through all peptides
         for aa in aaList:
         #iterate through all amino acids
            newPeptide = peptide.copy()
            newPeptide.append(aa)
            newPeptides.append(newPeptide)
      return newPeptides


   def cyclospectrum(self, peptide, cyclic) :
      '''
      Note:
         This method calculates the spectrum of the peptide to see if it is
         identical to the given spectrum in cyclopeptideSequencing.
      Args:
         peptide: the current peptide being checked in cyclopeptideSequencing
      Returns:
         spectrum: the spectrum of the peptide
      '''
      spectrum = [0] #start 0 since it is included in spectrum
      for i in range(1, len(peptide)):
      #for all masses starting after the first

         for x in range(len(peptide)):
         #for all masses
            if i+x <= len(peptide):
            #check if in range
               spectrum.append(sum(peptide[x:x+i]))  #add linear masses
            
            elif(cyclic): #if final check
               y = i+x-len(peptide)
               #update range
               spectrum.append(sum(peptide[x:])+sum(peptide[:y])) #add cyclic masses

      spectrum.append(sum(peptide))
      spectrum.sort() #sort the spectrum
      return spectrum


   def cyclopeptideSequencing(self, spectrum) :
      '''
      Note:
         Finds valid sequences of amino acid weights to match the spectrum
      Args:
         spectrum: a list of the spectrum of the peptide
      Returns:
         Calls printPeptide to output peptides that match the spectrum.
      '''

      spectrumDict = defaultdict(int)
      #A dict of all the mass occurences in the spectrum
      for mass in spectrum:
         spectrumDict[mass] += 1

      peptides = [[]] #a set containing only the empty peptide
      while peptides != []:
      #while Peptides is nonempty
         
         peptides = self.expand(peptides) #Peptides ← Expand(Peptides)
         newPeptides = [] #emtpy for new peptides that pass checks
         
         for peptide in peptides[::-1]:
         #for each peptide Peptide in Peptides in reverse for faster removal

            if sum(peptide) == spectrum[-1]:
            #if Mass(Peptide) = ParentMass(Spectrum)

               spec = self.cyclospectrum(peptide, True) #calculates spectrum
               
               if spec == spectrum:
               #if Cyclospectrum(Peptide) = Spectrum

                  self.printPeptide(peptide) #output the peptide

            elif sum(peptide) < spectrum[-1]:  #check if under maximum mass

               spec = self.cyclospectrum(peptide, False) #calculates spectrum
               match = True  #assume it matches
               for mass in spec:
               
                  if spec.count(mass) > spectrumDict[mass]:
                  #check if each mass is valid

                     match = False #if not match is false

               if match: #if match is true add to newPepides
                  newPeptides.append(peptide)

            peptides = newPeptides   #update peptides with newPeptides



   def printPeptide(self, peptide) :
      '''
      Note:
         Prints the peptides matching the spectrum in the required format
      Args:
         peptide: a list of an amino acid spectrum 
      Returns:
         Prints the peptides in a string seperated by '-'
      '''
      pep = ""
      for aa in peptide:
         pep += str(aa) + '-'
      print(pep[:-1])


def main() :

   # read file and store as matrices
   fileReader = FileReader('')
   spectrum = fileReader.readFile()

   # initialize hmm object and print the probability
   CycloPeptideSequencing(spectrum)


if __name__ == "__main__":
   main()
