#!/usr/bin/env python3
# Max Genetti (mgenetti)

'''
This program implements Viterbi learning to generate the transition and
emission probabilities from a sequence.

Main expects a file formatted for rosalind problem 21 handled by a general file
reader class FileReader.

Given: A sequence of emitted symbols x = x1 ... xn in an alphabet A, generated
by a k-state HMM with unknown transition and emission probabilities, initial
Transition and Emission matrices and a number of iterations i.

Return: A matrix of transition probabilities Transition and a matrix of
emission probabilities Emission that maximizes Pr(x, π) over all possible
transition and emission matrices and over all hidden paths π. 
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
   Args:
      iteration: an int for the number of iteraions
      observation: a string of the observation
      alphabet: a list of the alphabet in a hmm
      states: a list of the states in a hmm
      transitionMatrixList: a list of lists containing the transition matrix
      emissionMtrixList: a list of lists containing the emission matrix
   Returns:
      A matrix of transition probabilities and a matrix of emission
      probabilities
   '''

   def __init__(self, iteration, observation, alphabet, states, transitionMatrixList, emissionMatrixList) :
      '''Stores iteration, observation, alphabet, states and then builds the matrices'''
      self.iteration = iteration
      self.observation = observation
      self.alphabet = alphabet
      self.states = states
      self.path = ''
      self.transition = self.buildMatrix(transitionMatrixList)
      self.emission = self.buildMatrix(emissionMatrixList)


   def iterator(self) :
      '''
      Iterates through rebuilding the emission and transition matrices based on viterbi
      '''
      for i in self.iteration:
         self.viterbi()
         self.buildEmission()
         self.buildTransition()
      self.printMatrices()


   def buildMatrix(self, matrixList) :
      '''
      creates the numpy array of the matrix
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


   def convertOutput (self, outputIndex) :
      '''
      iterate through output converting number to state based on index
      Returns the viterbi output converted into string of state charecters
      '''
      output = ''
      for i in outputIndex:
         output += self.states[i]
      return output


   def viterbi(self):
      '''
      Returns the path of a given observation sequence using the viterbi algorithm
      '''
      #convert observed sequence into an array of the indexed values
      obsArray = np.array(self.convertObservation())

      # sets prior probability as equal for all possible starting observations
      prior = 1 / len(self.states)
      #set to length of the obersation obsLen
      obsLen = len(self.observation)

      #Initialize tracking tables of the needed dimensions
      track1 = np.empty((len(self.states), obsLen))
      track2 = np.empty((len(self.states), obsLen))

      #initialize tracking table1 with the prior probability
      track1[:, 0] = prior * self.emission[:, obsArray[0]]
      track2[:, 0] = 0
      # Iterate throught the observations updating the tracking tables
      for i in range(1, obsLen): #starting at 1 to avoid overindexing
         track1[:, i] = np.max(track1[:, i - 1] * self.transition.T * self.emission[np.newaxis, :, obsArray[i]].T, 1)
         track2[:, i] = np.argmax(track1[:, i - 1] * self.transition.T, 1)

      # initialize the output array
      x = np.empty(obsLen, 'B')
      x[-1] = np.argmax(track1[:, obsLen - 1])
      for i in reversed(range(1, obsLen)):
         x[i - 1] = track2[x[i], i]

      self.path = self.convertOutput(x)  # converts from numbers back to state chars



   def convertObservation (self) :
      '''
      iterate through observation storing alphabet index in list
      Returns the numeric equivalent of the observation
      '''
      num = []
      for i in self.observation:
         num.append(self.alphabet.index(i))
      return num


   def convertPath (self) :
      '''
      iterate through observation storing alphabet index in list
      Returns the numeric equivalent of the observation
      '''
      num = []
      for i in self.path:
         num.append(self.states.index(i))
      return num


   def buildTransition (self) :
      '''
      builds the transition matrix based on the path and possible states
      '''
      k =len(self.states)
      m = [[0]*k]*k   #make list of lists from number of states
      m = np.array(m)   #build np array from list of lists to hold counts
      matrix = m.astype(np.float) #make matrix for holding probabilities
      seq = self.convertPath()
      
      for x in range(len(seq)-1): #count each transition
         m[seq[x]][seq[x+1]] += 1
      
      s = np.sum(m, axis = 1) #sum the counts for each row
  
      for i in range(k): #converts counts into probability
         if s[i] == 0:
            matrix[i] = [1/k] * k
         else:
            for j in range(k):
               matrix[i, j] = m[i, j]/s[i]

      self.transition = matrix


   def buildEmission (self) :
      '''
      builds the emission matrix based on the path and possible states
      '''
      k =len(self.states)
      l =len(self.alphabet)

      m = [[0]*l]*k   #make list of lists from number of states
      m = np.array(m)   #build np array from list of lists to hold counts
      matrix = m.astype(np.float) #make matrix for holding probabilities
      path = self.convertPath()
      obs = self.convertObservation()
      
      for x in range(len(path)): #count each emission
         m[path[x]][obs[x]] += 1

      s = np.sum(m, axis = 1) #sum the counts for each row

      for i in range(k): #converts counts into probability
         if s[i] == 0:
            matrix[i] = [1/k] * l
         else:
            for j in range(l):
               matrix[i,j] = m[ i, j]/s[i]

      self.emission = matrix


   def printMatrices (self) :
      '''
      Prints the emission and transmission matrices formatted
      '''
      e = self.emission.astype(str)
      t = self.transition.astype(str)
      print(" " +" ".join(self.states)) #prints transition matrix
      for i in range(len(t)):
         print(self.states[i]+" "+" ".join(t[i]))
      print('--------')
      print(" "+" ".join(self.alphabet)) #prints emission matrix
      for i in range(len(e)):
         print(self.states[i]+" "+" ".join(e[i]))




def main() :

   # read file and store as matrices
   fileReader = FileReader('')
   matrices = fileReader.readFile()

   # initialize hmm object and print the probability
   hmm = HMM(matrices[0][0][0], matrices[1][0][0], matrices[2][0], matrices[3][0], matrices[4], matrices[5])
   hmm.iterator()


if __name__ == "__main__":
   main()
