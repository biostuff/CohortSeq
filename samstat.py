#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 15:26:19 2015

@author: wolfthomas
"""

import click
from joblib import Parallel  
from joblib import delayed
import sh
import os
import fnmatch

 #define a function to get bam files from the bam directory
def listdir(dirname, pattern="*"):
    return fnmatch.filter(os.listdir(dirname), pattern)



@click.command()
@click.option('--bamFolder',help='The folder where.')
@click.option('--ncores',default=1,help='How many cores should be used if multicore processing is activated ? (One core for each bam file)')
@click.option('--samstatLocation',default="samstat",help='Where is samstat located')

    
def getsamstat (bamFolder,multicore,ncores,samstatLocation): 
    """Wrapper for the samstat tool. Also adds support for multicore processing. One core for each bam file"""

   
    #get a list of bam files from the bam directory
    bamsToStat = listdir(dirname = bamFolder,pattern = "*.bam")

    #define a samstat python function using sh
    samstat = sh.Command(samstatLocation)

    #define a samstat function that returns the filename of the bam that was statted
    def samstatReturn(bamstat): 
        samstat(bamstat)
        return(bamstat)


    if isinstance(bamsToStat,list) == False:
        bamsToStat = [bamsToStat]
        
    #samstat can also accept multiple bam files at once
    #try instead of a loop    
    if ncores < 2:    
        for i in range(len(bamsToStat)):
            samstat(bamFolder + bamsToStat[i])
    else:
        Parallel(n_jobs=ncores)(delayed(samstatReturn)(bamFolder + i) for i in bamsToStat)

        
#generate the samstat output   
if __name__ == '__main__':
    getsamstat()
