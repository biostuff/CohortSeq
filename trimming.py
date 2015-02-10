# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 11:22:28 2015

@author: wolfthomas
"""


import click
import sh



@click.command()
@click.option('--sequences',help='The bam or fastq files files where the sequences to be aligned are stored. For mutiple bams please seperate them by commas.')
@click.option('--sequencesReverse',help='The fastq files  where the reverse sequences to be aligned are stored. For mutiple bams please seperate them by commas.')
@click.option('--outputNames',help='The bam,sam or fastq files files where the sequences to be aligned are stored. For mutiple bams please seperate them by commas.')
@click.option('--ncores',default=1,help='How many cores should be used if multicore processing is activated ?')
@click.option('--maxmem',default=1,help='How many cores should be used if multicore processing is activated ?')
@click.option('--outoutNamesReverse',help='The fastq files  where the reverse sequences to be aligned are stored. For mutiple bams please seperate them by commas.')
@click.option('--bam',is_flag=True,help='Is the data in bam format ? Currently only supported for single end data.')
@click.option('--paired',is_flag=True,help='Is the data in bam format ? Currently only supported for single end data.')
@click.option('--bbdukLocation',default="bbduk.sh",help='Where is the bbduk sh located.')
@click.option('--trimseqLocation',default = "truseq.fa.gz",help='Where are the bbtools  located.')



def readTrimming (sequences,sequencesReverse,outputNames,outputNamesReverse,ncores,maxmem,outoutNamesReverse,bam,paired,bbdukLocation,trimseqLocation):
    sequences = sequences.split(",")
    outputNames = outputNames.split(",")    
    
    if paired == True:
        sequencesReverse = sequencesReverse.split(",")
        outputNamesReverse = outputNamesReverse.split(",")        
    
    bbduk = sh.Command(bbdukLocation)

    for i in enumerate(sequences):
        if paired == True and bam == False:
            bbduk("-Xmx" + str(maxmem) + "g","in1=" + sequences[i],"in2=" + sequencesReverse[i],"out1=" + outputNames[i],"out2=" + outputNamesReverse[i],"qtrim=" + "rl","trimq=" + str(10),"ktrim=r","k=25","mink=11","ref=" + trimseqLocation,"hdist=1","threads=" + str(ncores))   
        else:
            bbduk("-Xmx" + str(maxmem) + "g","in=" + sequences[i],"out=" + outputNames[i],"qtrim=" + "rl","trimq=" + str(10),"ktrim=r","k=25","mink=11","ref=" + trimseqLocation,"hdist=1","threads=" + str(ncores))   



if __name__ == '__main__':
    readTrimming()
