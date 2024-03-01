# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 10:22:38 2024

@author: pbrown2
"""

from itertools import chain, combinations
import numpy as np
from scipy import optimize as opt

def needed_subsets(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)))



m = 4 # number of row actions
n = 3 # number of column actions

row_action_names = {idx:'a'+str(idx) for idx in range(m)}
row_action_indices = {row_action_names[name]:name for name in row_action_names}

col_action_names = {idx:'b'+str(idx) for idx in range(n)} 
col_action_indices = {col_action_names[name]:name for name in col_action_names}

def generate_decision_variables() :
    # generates all the info for the problem's main decision variables (other than t)
    A = {idx:set() for idx in row_action_names}
    B = {idx:set() for idx in col_action_names}
    D = {}
    R = {}
    k = -1 # the -1 is for (),() action
    row_subsets = needed_subsets(row_action_names)
    col_subsets = list(needed_subsets(col_action_names))
    for row_resource in row_subsets : # row_resource is a tuple of row actions containing this resource
        for col_resource in col_subsets : # col_resource is a tuple of col actions containing this resource 
            # these resources specify the kth decision variable
            for action in row_resource :
                A[action].add(k)
            for action in col_resource :
                B[action].add(k)
            resource_description = ''.join(row_action_names[i] for i in row_resource)
            resource_description += ''.join(col_action_names[i] for i in col_resource)
            D[k] = resource_description
            R[resource_description] = k
            k += 1
    for action in A :
        A[action].discard(-1)
    for action in B :
        B[action].discard(-1)
    del D[-1]
    del R['']
    return A,B,D,R,k
                
def generate_constraint_matrices(Phi,A,B,k) :
    # given m x n array Phi specifying which values of potential function we want,
    # output Aeq beq Aineq bineq arrays specifying the LP constraints.
    # Phi: specified values are > -0.5, unspecified <= -0.5.
    # something will break if Phi is the wrong dimensions
    # k is # of resource constraints
    
    Aeq = []
    beq = []
    Aineq = []
    bineq = []
    
    Cineq = [] # list to store the (i,j) pair corresponding to each ineq constraint
    Ceq = [] # list to store the (i,j) pair corresponding to each eq constraint
    
    for i in range(m) :
        for j in range(n) :
            a = np.zeros(k+1) # initialize constraint matrix row
            if Phi[i,j] > -0.5 :
                a[list(A[i] | B[j]) ] = 1 # add in resource values corresponding to action pair (i,j)
                Aeq.append(a)
                beq.append(Phi[i,j])
                Ceq.append((i,j))
            else : 
                a[list(A[i] | B[j]) ] = 1 # add in resource values corresponding to action pair (i,j)
                a[-1] = -1
                Aineq.append(a)
                bineq.append(0)
                Cineq.append((i,j))
    return np.array(Aeq), np.array(beq), np.array(Aineq), np.array(bineq), Ceq, Cineq
                
if __name__ == "__main__" :
    Phi = [[1,-1,-1],
           [.9,.8,-1],
           [-1,.7,.6],
           [-1,-1,.5]]
    Phi = np.array(Phi)  
    
    # This code needs some stuff:
        # refactor as a class
        # read m,n from Phi
        # have nice printouts to facilitate making games
    
    A,B,D,R,k = generate_decision_variables()
    Aeq,beq,Aineq,bineq,Ceq,Cineq = generate_constraint_matrices(Phi,A,B,k)
    c = np.zeros(k+1)
    c[-1] = 1
    
    answers = opt.linprog(c,
                          A_ub = Aineq,
                          b_ub = bineq,
                          A_eq = Aeq,
                          b_eq = beq,
                          bounds = (0,None))
    print(answers)
    
    
    
