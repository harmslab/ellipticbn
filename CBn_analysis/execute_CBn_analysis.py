#!/usr/bin/env python3
# coding: utf-8

"""
Command line interface to use the CBn_analysis.run_all
"""

from CBn_analysis import run_all
import sys
import os
import glob

__usage__ = "No xyz file or folder of xyz files specified. \n \n Please specify either an xyz file or folder of xyz files"

def main(argv=None):
   
    try:
        file_input = sys.argv[1]
    except IndexError:
        err = __usage__
        
        raise ValueError(err)
        
    if os.path.isfile(file_input):
        file_input = [file_input]
        
    elif os.path.isdir(file_input):
        file_input = glob.glob(os.path.join(file_input,"*.xyz"))
        
    else:
        err = __usage__
        
        raise ValueError(err)
        
    for f in file_input:
        run_all(f)   # execute all functions and analyze xyz file
    
       
if __name__ == "__main__":
    main()