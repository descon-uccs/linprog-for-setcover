# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 10:22:38 2024

@author: pbrown2
"""

from potential_game import *

if __name__ == "__main__" :
    Phi = [[1,-1,-1],
           [.9,.8,-1],
           [-1,.7,.6],
           [-1,-1,.5]]
    Phi = np.array(Phi)  
    
        
    pg = Potential_Game(Phi)
    pg.solve_potential_function()
    print(pg.answers)
    
    
    
