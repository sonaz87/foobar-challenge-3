# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 10:48:11 2021

@author: lor

Doomsday Fuel
=============

Making fuel for the LAMBCHOP's reactor core is a tricky process because of the exotic matter involved. 
It starts as raw ore, then during processing, begins randomly changing between forms, eventually reaching a stable form. 
There may be multiple stable forms that a sample could ultimately reach, not all of which are useful as fuel. 

Commander Lambda has tasked you to help the scientists increase fuel creation efficiency by predicting the end state of a given ore sample. 
You have carefully studied the different structures that the ore can take and which transitions it undergoes. 
It appears that, while random, the probability of each structure transforming is fixed. That is, each time the ore is in 1 state, 
it has the same probabilities of entering the next state (which might be the same state).  
You have recorded the observed transitions in a matrix. The others in the lab have hypothesized more exotic forms that the ore can become, 
but you haven't seen all of them.

Write a function solution(m) that takes an array of array of nonnegative ints representing how many times that state has gone to the next state 
and return an array of ints for each terminal state giving the exact probabilities of each terminal state, 
represented as the numerator for each state, then the denominator for all of them at the end and in simplest form. 
The matrix is at most 10 by 10. It is guaranteed that no matter which state the ore is in, there is a path from that state to a terminal state. 
That is, the processing will always eventually end in a stable state. The ore starts in state 0. 
The denominator will fit within a signed 32-bit integer during the calculation, as long as the fraction is simplified regularly. 

For example, consider the matrix m:
[
  [0,1,0,0,0,1],  # s0, the initial state, goes to s1 and s5 with equal probability
  [4,0,0,3,2,0],  # s1 can become s0, s3, or s4, but with different probabilities
  [0,0,0,0,0,0],  # s2 is terminal, and unreachable (never observed in practice)
  [0,0,0,0,0,0],  # s3 is terminal
  [0,0,0,0,0,0],  # s4 is terminal
  [0,0,0,0,0,0],  # s5 is terminal
]
So, we can consider different paths to terminal states, such as:
s0 -> s1 -> s3
s0 -> s1 -> s0 -> s1 -> s0 -> s1 -> s4
s0 -> s1 -> s0 -> s5
Tracing the probabilities of each, we find that
s2 has probability 0
s3 has probability 3/14
s4 has probability 1/7
s5 has probability 9/14
So, putting that together, and making a common denominator, gives an answer in the form of
[s2.numerator, s3.numerator, s4.numerator, s5.numerator, denominator] which is
[0, 3, 2, 9, 14].

Languages
=========

To provide a Java solution, edit Solution.java
To provide a Python solution, edit solution.py

Test cases
==========
Your code should pass the following test cases.
Note that it may also be run against hidden test cases not shown here.

-- Java cases --
Input:
Solution.solution({{0, 2, 1, 0, 0}, {0, 0, 0, 3, 4}, {0, 0, 0, 0, 0}, {0, 0, 0, 0,0}, {0, 0, 0, 0, 0}})
Output:
    [7, 6, 8, 21]

Input:
Solution.solution({{0, 1, 0, 0, 0, 1}, {4, 0, 0, 3, 2, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}})
Output:
    [0, 3, 2, 9, 14]

-- Python cases --
Input:
solution.solution([[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0,0], [0, 0, 0, 0, 0]])
Output:
    [7, 6, 8, 21]

Input:
solution.solution([[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])
Output:
    [0, 3, 2, 9, 14]


"""

def solution(m):
    
    import fractions
    import numpy as np

    def invert(A, Identity):
        A = np.array(A)
        m,n = A.shape
        if m != n:
            raise Exception("The input matrix should be a square")
        
        # create an Identity matrix of size mxm
        I = Identity
    
        # Augment matrix A
        aug_A = np.c_[A,I]
    
        #Gaussian Elimination to generate an upper triangular matrix
        j = 0
        for i in range(m-1):
            pivot = aug_A[i][j]
            if pivot == 0:
            
                found = False
                for k in range(i+1,m):
                    if aug_A[k][j] != 0:
                        temp = aug_A[k].copy()
                        aug_A[k] = aug_A[i].copy()
                        aug_A[i] = temp.copy()
                        found = True 
                        break
                if found == False:
                    raise Exception("The matrix is singular and hence cannot be inverted")
                else:
                    pivot = aug_A[i][j]
            for k in range(i+1, m):
                target = aug_A[k][j]
                multiplier = target / pivot
                aug_A[k] = aug_A[k] - multiplier*aug_A[i]
            j += 1
    
        #Generate 0s above the pivot and create a diagonal matrix
        j = m-1
        for i in range(m-1,0,-1):
            pivot = aug_A[i][j]
            for k in range(i-1,-1,-1):
                target = aug_A[k][j]
                multiplier = target / pivot
                aug_A[k] = aug_A[k] - multiplier*aug_A[i]
            j -= 1
    
        for i in range(m):
            aug_A[i] /= aug_A[i][i]
    
        return aug_A[:,m:]

    def least_common_multiplier(rationals):
        
        denominators = [fractions.Fraction(r).denominator for r in rationals]
        
        lcm = denominators[0]
        for d in denominators[1:]:
            lcm = lcm // fractions.gcd(lcm, d) * d
        
        return lcm    
    
    # creating fractional probability matrix        
    
    matrix = []
    for line in m:
        newline = []
        counter = 0
        for element in line:
            counter += element
        for element in line:
            if counter != 0:
                newelement = fractions.Fraction(element, counter)
            else:
                newelement = fractions.Fraction(0)
            newline.append(newelement)
        matrix.append(newline)
        
        
    a = np.array([[matrix[i][j] for j in range(len(m))] for i in range(len(m))])
    
    print(a)
      
    # solving the probability matrix

    new_order = []
    terminals = []
    for i in range(len(a)):
    
        terminal = True
        for j in a[i]:
            if j != 0:
                terminal = False
        if terminal == True:
            a[i][i] = fractions.Fraction(1)
            terminals.append(i)
            new_order.append(i)
    
    for i in range(len(a)):
        if i not in new_order:
            new_order.append(i)

    
    A1 = np.array([[a[i][j] for j in new_order] for i in new_order])

      
    # print(b)
    
    
    # print(terminals)
    
    I = np.array([[A1[i][j] for j in range(len(terminals))] for i in range(len(terminals))])
    O = np.array([[A1[i][j] for j in range(len(terminals), len(A1))] for i in range(len(terminals))])
    R = np.array([[A1[i][j] for j in range(len(terminals))] for i in range(len(terminals), len(A1))])
    Q = np.array([[A1[i][j] for j in range(len(terminals), len(A1))] for i in range(len(terminals), len(A1))])
    I2 = np.array([[A1[i][j] for j in range(len(Q))] for i in range(len(Q[0]))])
    
    # print("I")
    # print(I)
    # print("O")
    # print(O)
    print("R")
    print(R)
    print("Q")
    print(Q)
    # print("I2")
    # print(I2)
    
    F = (I2-Q)
    print("i-q")
    print(F)
    # TODO! create inverted matrix of F with gauss elimination
    F = invert(F, I2)
    print("F")
    print(F)
    # print("F")
    # print(F)
    result = (np.matmul(F,R))
    print("result")
    print(result)
    
    solution = result[0]
    denom = least_common_multiplier(solution)
    # print(denom)
    final = []
    for i in solution:
        final.append(int(i.numerator * (denom/i.denominator)))
        
    final.append(denom)
    return final

    


        


test = [[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0,0], [0, 0, 0, 0, 0]]
test = [[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], [0, 0, 0, 0, 0], [0, 0, 0, 0,0], [0, 0, 0, 0, 0]]  # s5 is terminal
test6 = [
        [0, 7, 0, 17, 0, 1, 0, 5, 0, 2],
        [0, 0, 29, 0, 28, 0, 3, 0, 16, 0],
        [0, 3, 0, 0, 0, 1, 0, 0, 0, 0],
        [48, 0, 3, 0, 0, 0, 17, 0, 0, 0],
        [0, 6, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]

test9 = [
        [0, 86, 61, 189, 0, 18, 12, 33, 66, 39],
        [0, 0, 2, 0, 0, 1, 0, 0, 0, 0],
        [15, 187, 0, 0, 18, 23, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]    

test10 = [
        [0, 0, 0, 0, 3, 5, 0, 0, 0, 2],
        [0, 0, 4, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 4, 4, 0, 0, 0, 1, 1],
        [13, 0, 0, 0, 0, 0, 2, 0, 0, 0],
        [0, 1, 8, 7, 0, 0, 0, 1, 3, 0],
        [1, 7, 0, 0, 0, 0, 0, 2, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]
        
test2 = [[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]

result = solution(test2) 
print(result)

    
    