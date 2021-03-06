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

#sequences = "/media/wolfthomas/BackupData/PipelineTest/Input/"
#folder = True
#samstatlocation = "/home/wolfthomas/NGSTools/samstat-1.5/src/samstat"

#getsamstat(sequences = sequences,folder = True,ncores = 15,samstatlocation = samstatlocation)

def samstatReturn(bamstat,samstatlocation): 
        samstat = sh.Command(samstatlocation)
        samstat(bamstat)
        return(bamstat)




@click.command()
@click.option('--sequences',help='The filename of a bam file, ot if --folder a folder containing the bams.')
@click.option('--folder',is_flag=True,help='Should all bams in a folder, as specified by --sequences, be statted.')
@click.option('--ncores',default=1,help='How many cores should be used if ? (Maximum one core for each bam file)')
@click.option('--samstatlocation',default="samstat",help='Where is samstat located ?')



    
def getsamstat (sequences,folder,ncores,samstatlocation): 
    """Wrapper for the samstat tool. Also adds support for multicore processing. One core for each bam file"""

   
    #define a function to get bam files from the bam directory
    def listdir(dirname, pattern="*"):
        return fnmatch.filter(os.listdir(dirname), pattern)

   
   
    #get a list of bam files from the bam directory
    if folder == True:
        bamsToStat = listdir(dirname = sequences,pattern = "*.bam")
        bamsToStat = [sequences + i for i in bamsToStat]
    else:
        bamsToStat = sequences


    

    if isinstance(bamsToStat,list) == False:
        bamsToStat = [bamsToStat]
        
    #samstat can also accept multiple bam files at once
    #try instead of a loop    
    if ncores < 2:    
        for i in range(len(bamsToStat)):
            samstatReturn(bamsToStat[i],samstatlocation)
    else:
        Parallel(n_jobs=ncores)(delayed(samstatReturn)(i,samstatlocation) for i in bamsToStat)

        
#generate the samstat output   
if __name__ == '__main__':
    getsamstat()
