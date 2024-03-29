# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 10:22:38 2024

@author: pbrown2
"""

from potential_game import *
from tabulate import tabulate

if __name__ == "__main__" :
    # Phi = [[10,-1,-1,-1],
    #        [8,6,-1,-1],
    #        [-1,5,5,5],
    #        [-1,-1,5,5]]
    
    Phi = [[10,6,-1],
           [8,6,4],
           [4,4,2]]
           # [-1,-1,5,5]]
    
    Phi = np.array(Phi)
    

    pg = Potential_Game(Phi)
    pg.solve_potential_function()
    print(pg.answers)
    
    pg.print_resource_values()