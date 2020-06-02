#!/usr/bin/env python3
# Max Genetti (mgenetti)

'''
This program generates a DeBruijn graph using specified lengths of a given string

Main expects a file formatted for rosalind problem 11 handled by a general file
reader class FileReader.

Given: An integer k and a string Text.

Return: DeBruijnk(Text), in the form of an adjacency list.
'''

import sys

class FileReader :
   '''
   Args:
      fname: a file formatted with an int in the top line and a string in the
      second line
      
   Returns:
      The int in the first line and the sequence in the second line
   '''
   
   def __init__ (self, fname='') :
      '''constructor: saves attribute fname'''
      self.fname = fname
            
   def doOpen (self) :
      '''opens the file'''
      if self.fname is '':
         return sys.stdin
      else:
         return open(self.fname)
 
   def readFile (self) :
      '''
      Parses the file and returns the int in the first line and sequence in
      the second
      '''
      kmer = 0
      sequence = ''
      with self.doOpen() as fileH:
         for line in fileH:
            if kmer == 0:
               kmer = line.rstrip()
            else:
               sequence += ''.join(line.rstrip().split()).upper()
      yield kmer, sequence

class sequenceParser :
   """
   Generates the k-mer composition of a string.

   Args:
      entry: An indexable object where entry[0] is int and entry[1] is string
   
   Returns:
      A DeBruijn in the form of an adjacency list of the sequences of length
      entry[0] from string entry[1]
   """
   
   def __init__ (self, entry) :
      '''constructor: saves attribute kmer and sequence'''   
      self.kmer = int(entry[0])
      self.sequence = entry[1]

   def parse(self) :
      '''Returns parsed sequence in lexigraphic order'''
      kList = []
      for i in range(0, len(self.sequence)  - self.kmer +1):
         kList.append(self.sequence[i:i+self.kmer])
      kList.sort()
      return kList

   def pairItUp(self) :
      '''Returns paired sequence dictionary'''
      seqList = self.parse()
      pairD= {}

      for item in seqList:
         if item[0:-1] not in pairD:
            pairD[item[0:-1]] = (item[1:])
         else: pairD[item[0:-1]] += (', '+item[1:])
      return pairD
  

def main():
   
   '''Reads the file given in systdin'''
   fileReader = FileReader()
   seqFile = fileReader.readFile()

   '''generates a kmer and sequence in the file'''
   for entry in seqFile:
      x = entry

   '''Initializes the class object sequenceParser and calls pairItUp'''
   seqs = sequenceParser(x)
   pairDict = seqs.pairItUp()

   '''prints all the elements in the dictionary as: key -> value'''
   for key in pairDict:
      print(key,' -> ', pairDict[key])

if __name__ == "__main__":
   main()
