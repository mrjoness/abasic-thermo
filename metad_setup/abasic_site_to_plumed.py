import os
import sys
import subprocess
import numpy as np
import glob

def get_pairs(num_bp, b_idx, strand=1):
    '''return matching for abasic strands'''
    
    atom_pairs = []
    
    # account for control
    if b_idx is None:
        
        s1 = [i*3 + 1 for i in range(num_bp)]
        s2 = [i*3 + num_bp*3 for i in range(num_bp)]
        
        for idx, (b1, b2) in enumerate(zip(s1, reversed(s2))):
            atom_pairs.append((b1, b2))

    else:
        
        # shift index down for missing bp
        s1 = [i*3 + 1 for i in range(b_idx-1)] + [i*3 for i in range(b_idx-1, num_bp)]
        s2 = [i*3 + num_bp*3 - 1 for i in range(num_bp)]

        for idx, (b1, b2) in enumerate(zip(s1, reversed(s2))):

            # just append same base for missing index
            if idx != b_idx-1:
                atom_pairs.append((b1, b2))
    
    # atom_pairs is list of tuples
    plumed_str = 'd1: DISTANCES'
    
    # add each pair with plumed ormat
    for i, pair in enumerate(atom_pairs):
        plumed_str += f' ATOMS{i}={pair[0]},{pair[1]}'
        
    plumed_str += ' MEAN\n'
    
    return plumed_str
    

atom_pairs = get_pairs(num_bp, b)
