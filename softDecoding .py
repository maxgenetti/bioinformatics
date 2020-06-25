#!/usr/bin/env python3
# Max Genetti (mgenetti)

'''
This program generates the probability of a sequence being in a possible state
given its emission sequence, emission matrices and transition matrices.

Main expects a file formatted for rosalind problem 22 handled by a general file
reader class FileReader.

Given: A string x, followed by the alphabet Σ from which x was constructed,
followed by the states States, transition matrix Transition, and emission
matrix Emission of an HMM (Σ, States, Transition, Emission).

Return: The probability Pr(πi = k|x) that the HMM was in state k at step i
(for each state k and each step i). 
'''

import sys
import numpy as np


class FileReader :
   '''

   Args:
      fname: a file formatted with a string in each line
      
   Returns:
      The list parsed by lines starting with - into a lists containing each line between.
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
      matrices = []
      matrix = []
      with self.doOpen() as fileH:
         for line in fileH:
            if line.startswith('-') : #seperates lists on lines starting with -
               matrices.append(matrix)
               matrix = []
               continue
            else:
               line = line.split()  #saves each line as a list
               matrix.append(line) #saves lines grouped together as a list of lists
         matrices.append(matrix)
      return matrices


class HMM:
   '''
   Note:
      Recieved help from Alex Bagi on Output Method
   Args:
      observation: a string of the observation
      alphabet: a list of the alphabet in a hmm
      states: a list of the states in a hmm
      transitionList: a list of lists containing the transition matrix
      emissionList: a list of lists containing the emission matrix
   Returns:
      A list of the probability of the states at each position in the observation
   '''

   def __init__(self, observation, alphabet, states, transitionMatrixList, emissionMatrixList) :
      '''
      Stores observation, alphabet, states and then builds the matrices
      '''
      self.observation = observation
      self.alphabet = alphabet
      self.states = states
      self.transition = self.buildMatrix(transitionMatrixList)
      self.emission = self.buildMatrix(emissionMatrixList)
      self.observationIndex = self.convertObservation()
      self.len = len(self.observation)


   def buildMatrix(self, matrixList) :
      '''
      creates in the numpy array of the matrix
      '''
      m = []
      #iterate through, skipping the list and first item in each following list
      for i in range(1,len(matrixList)):
         m.append(matrixList[i][1:])
      #append new lines for each row of the matrix
      m = np.array(m)
      #set alll the values to float
      matrix = m.astype(np.float)
      return matrix


   def convertObservation (self) :
      '''
      iterate through observation storing alphabet index in list
      Returns the numeric equivalent of the observation
      '''
      num = []
      for i in self.observation:
         num.append(self.alphabet.index(i))
      return num



   def forward(self):
      '''
      calculates the probability of the sequence
      '''
      observed = self.convertObservation()   #convert observed into numberic representation
      prior = 1/len(self.states)   #starting probability of that state
      
      #create empty array to hold calculations for each poisition in the observation
      fwd=np.zeros((len(self.observation),len(self.states)))   
      
      #initializes with the starting probability prior
      fwd[0]=prior*self.emission[:,observed[0]]   
     
      #Iterate through multiplying the probability for each position
      for t in range(len(self.observation)-1):   
         fwd[t+1]=np.dot(fwd[t],self.transition)*self.emission[:,observed[t+1]]
      #returns the sum
      return fwd, fwd[-1].sum()


   def backward(self):
      '''
      Calculates the backwords probability from a given observation using the
      emission and transition matrices
      '''
      obsArray = np.array(self.observationIndex) #converts the observation into an Array
      matrix = np.zeros((self.len, len(self.states))) #creates an empty matrix
      for i in range(len(self.states)): #iterates through states storing positions in matrix
         matrix[self.len - 1, i] = 1
      for i in reversed(range(self.len - 1)): #iterates through the observation in reverse
         matrix[i] = np.dot(matrix[i + 1] * self.emission[:, obsArray[i + 1]], self.transition.T)
      return matrix


   def output(self):
      '''
      Calculates the probability for the states at each position 
      Note: Recieved help from Alex Bagi on this method
      '''
      fwd = self.forward() #calculate forward probability
      bwd = self.backward() #calculate reverse probability

      prob = (fwd[0] * bwd) / fwd[1] #calculates the initial probability

      print(*self.states, sep='\t') #prints states
      for s in range(0, len(prob)): #prints the probability of each state
         print('{0}'.format('\t'.join([str(number) for number in prob[s, :]])))


def main() :

   # read file and store as matrices
   fileReader = FileReader('')
   matrices = fileReader.readFile()

   # initialize hmm object and print the probability
   hmm = HMM(matrices[0][0][0], matrices[1][0], matrices[2][0], matrices[3], matrices[4])
   hmm.output()
   


if __name__ == "__main__":
   main()
