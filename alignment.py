# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 11:06:20 2015

@author: wolfthomas
"""


import click
import sh
import os
import subprocess
import tempfile
import shutil
import math
from fs.osfs import OSFS
from fs.memoryfs import MemoryFS
from fs.expose import fuse


@click.command()
@click.option('--mapper',default="ngm",help='Which aligner should be used. Currently supports NGM and BWA.')
@click.option('--sequences',help='The bam,sam or fastq files files where the sequences to be aligned are stored. For mutiple bams please seperate them by commas.')
@click.option('--sequencesreverse',help='The fastq files  where the reverse sequences to be aligned are stored. For mutiple bams please seperate them by commas.')
@click.option('--bam',is_flag=True,help='Is the data in bam format ? Currently only supported for single end data.')
@click.option('--samplenames',help='If nor specified the filenames are used')
@click.option('--outputnames',help='The bam,sam or fastq files files where the sequences to be aligned are stored. For mutiple bams please seperate them by commas.')
@click.option('--reference',help='Where is the reference stored.')
@click.option('--paired',is_flag=True,help='Is the data paired end.')
@click.option('--markduplicates',default="none",help='Should PCR duplicates be marked. Either none,sambamba or samblaster. Only use sambalster if data is paired end.')
@click.option('--baq',is_flag=True,help='Should the BAQ tag be added to the bams.')
@click.option('--gpu',is_flag=True,help='Should GPU computing be used to speed up the alignment (currently only supported for NextGenMap(NGM))')
@click.option('--ncores',default=1,help='How many cores should be used if multicore processing is activated ?')
@click.option('--maxmem',default=1,help='How much memory is available for file processing')
@click.option('--baminmemory',is_flag=True,help='If set to True intermediary bams will be stored in memory. Make sure that maxmem + bam size is available.')
@click.option('--samtoolslocation',default="samtools",help='Where is samtools located')
@click.option('--sambambalocation',default="sambamba_latest",help='Where is sambamba located')
@click.option('--samblasterlocation',default="samblaster",help='Where is samblaster located')
@click.option('--bwalocation',default="bwa",help='Where is bwa located')
@click.option('--ngmlocation',default="ngm",help='Where is ngm located')
@click.option('--picardlocation',default="picard.jar",help='Where is the picard jar located.')
@click.option('--tempfolder',default="/tmp/",help='Give the file path to the temporary folder to be used.')
@click.option('--rg-pu',default="picard.jar",help='The sequencing platform unit used.')
@click.option('--rg-pl',default="picard.jar",help='The sequencing platform used.')
@click.option('--rg-id',default="missing",help='The read group id for each sample.')
@click.option('--rg-sm',default="missing",help='The samplename.')




def alignAndSort (mapper,sequences,sequencesreverse,bam,samplenames,outputnames,reference,paired,markduplicates,baq,gpu,ncores,maxmem,baminmemory,samtoolslocation,sambambalocation,samblasterlocation,bwalocation,ngmlocation,picardlocation,tempfolder,rg_pu,rg_pl,rg_id,rg_sm):

    #split the files into a list or if only one file define as a list
   #sequences = sequences.split(",")
   #outputnames = outputnames.split(",")    
    
   if paired == True:
       sequencesreverse = sequencesreverse.split(",")
      
    
   #BWA
   #bwa = sh.Command(bwaLocation)

   #NGM
   ngm  = sh.Command(ngmlocation)


   #Java
   javaString = "/usr/bin/java" + " " + "-XX:ParallelGCThreads=" + str(ncores) + " " + "-Xmx" + str(maxmem)  + "g" + " " + "-jar"
   #-Djava.io.tmpdir=/tmp

   #Picard
   picardString =  javaString  + " " + picardlocation 

     
   #define samtools as a python function
   samtools  = sh.Command(samtoolslocation)

   #sambamba_latest
   sambamba = sh.Command(sambambalocation)


   #generate a name for the temporary directory
   directoryTemp = tempfile.mkdtemp(dir = tempfolder)
  
   #generate a temporary bam either on disk or in memory
   #using fuse

   if baminmemory == True:
       memory_fs_dir = OSFS(directoryTemp)
       memfs = MemoryFS()
       memory_fs_dir.makedir('ramdrive', allow_recreate=True)
       memory_fused = fuse.mount(memfs,memory_fs_dir.getsyspath('ramdrive'))
       tempBam = memory_fused.path + "/" + "temp.bam"
   else:
      tempBam = directoryTemp + "/" + "temp.bam"

  
     
   #if the read group id is missing use the filename instead
   if(rg_sm == "missing"):
       pat_rg_sm = str(sequences) 
   else:
       pat_rg_sm = rg_sm
    
   #if the read group sample name is missing use filename instead
   if(rg_sm == "missing"):
       pat_rg_id = str(sequences) 
   else:
       pat_rg_id = rg_id
    
    
   if mapper == "NGM":   
       if paired == False:
           ngm("--rg-id",pat_rg_id,"--rg-sm",pat_rg_sm,"--rg-pu",rg_pu,"--rg-pl",rg_pl,"--rg-pg","NGM",t=ncores,r = reference,q = sequences,o = tempBam,b = True,g = gpu,paired = paired)
       else: 
           if markduplicates == "samblaster":
               outputter = open(tempBam,mode = "w")
               mapcores = round(ncores/2)
               if gpu == True:
                   gpu_string = "-g"
               else:
                   gpu_string = ""
                       
               p1 = subprocess.Popen([ngmlocation,"--rg-id",str(sequences),"--rg-sm",str(sequences),"--rg-pu",rg_pu,"--rg-pl",rg_pl,"--rg-pg","NGM","-1",sequences,"-2",sequencesreverse,"-r",reference,"-t",str(mapcores),gpu_string,"--paired"],stdout=subprocess.PIPE)
               p2 = subprocess.Popen([samblasterlocation],stdin=p1.stdout,stdout=subprocess.PIPE)
               p1.stdout.close()            
               p3 = subprocess.Popen([samtoolslocation,"view","-S","-b","-1","-@",str(mapcores),"-"], stdin=p2.stdout,stdout=outputter)
               p2.stdout.close()         
               p3.communicate()[0] #run our command
               outputter.close()  
           else:    
               ngm("--rg-id",pat_rg_id,"--rg-sm",pat_rg_sm,"--rg-pu",rg_pu,"--rg-pl",rg_pl,"--rg-pg","NGM","-1",sequences,"-2",sequencesreverse,r = reference,t=ncores,o = tempBam,b = True,strata = False,g = gpu,paired = paired,_no_out=False)    
            
               
           
   if mapper == "BWA": 
       R = "@RG\tID:" + pat_rg_id + "\tSM:" + pat_rg_sm   + "\tPU:" + rg_pu + "\tPL:" + rg_pl + "\tPG:" + "BWAMEM"
       outputter = open(tempBam,mode = "w")
       mapcores = round(ncores/2)
       if paired == False:
           if bam == True:
               p1 = subprocess.Popen([samtoolslocation,"bam2fq",sequences],stdout=subprocess.PIPE)  
               p2 = subprocess.Popen([bwalocation,"mem",reference,"-t",str(mapcores),"-R",R,"-"],stdin = p1.stdout,stdout=subprocess.PIPE)
               p1.stdout.close() #make sure we close the output so p2 doesn't hang waiting for more input
           else:
               p2 = subprocess.Popen([bwalocation,"mem",reference,sequences,"-t",str(mapcores),"-R",R,"-"],stdout=subprocess.PIPE)
       else:
           if markduplicates == "samblaster":
               p2bwa = subprocess.Popen([bwalocation,"mem",reference,sequences,sequencesreverse,"-t",str(mapcores),"-R",R,"-"],stdout=subprocess.PIPE)
               p2 = subprocess.Popen([samblasterlocation],stdin=p2bwa.stdout,stdout=subprocess.PIPE)
               p2bwa.stdout.close()
           else:
               p2 = subprocess.Popen([bwalocation,"mem",reference,sequences,sequencesreverse,"-t",str(mapcores),"-R",R,"-"],stdout=subprocess.PIPE)
       p3 = subprocess.Popen([samtoolslocation,"view","-S","-b","-1","-@",str(mapcores),"-"], stdin=p2.stdout,stdout=outputter)
       p2.stdout.close()
       p3.communicate()[0] #run our commands
       outputter.close()
  
    
       if paired == True:
           args = [picardString,"FixMateInformation","INPUT=",tempBam,"TMP_DIR=",directoryTemp,"VALIDATION_STRINGENCY=SILENT","ASSUME_SORTED=true"]
           subprocess.call([" ".join(args)],shell = True)
        

       #sort the bam file, calculate BAQ and generate bam
   if baq == True:
       if markduplicates == "sambamba":
           outputter = open(tempBam,mode = "r+") 
       else:
           outputter = open(outputnames,mode = "w")
       p1 = subprocess.Popen([samtoolslocation,"sort","-@",str(ncores),"-o","-l","0","-m",str(int(math.floor(maxmem/ncores)))+"G","-o",tempBam,directoryTemp + "/" + "sorttempbam"],stdout=subprocess.PIPE)
       p2 = subprocess.Popen([samtoolslocation,"fillmd","-Erub","-",reference], stdin=p1.stdout,stdout=subprocess.PIPE) #send p1's output to p2
       p1.stdout.close() #make sure we close the output so p2 doesn't hang waiting for more input
       p3 = subprocess.Popen([samtoolslocation,"view","-S","-b","-1","-@",str(mapcores),"-"], stdin=p2.stdout,stdout=outputter)
       p2.stdout.close()
       p3.communicate()[0] #run our command
       outputter.close()  
       if markduplicates != "sambamba":
           sambamba.index(outputnames,t = str(ncores))
         
   else:
       if markduplicates == "sambamba": 
           samtools.sort("-o",tempBam,"-O","bam","-@",str(ncores),"-m",str(int(math.floor(maxmem/ncores)))+"G","-T",directoryTemp + "/" + "sorttempbam/",tempBam)                  
       else:
           sambamba.sort(tempBam,o = outputnames,t = str(ncores),tmpdir = directoryTemp + "/" + "sorttempbam",m = str(maxmem) + "GB") 
        
       
       #mark  the duplicates
   if markduplicates == "sambamba":
       sambamba.markdup(tempBam,outputnames,t = ncores,tmpdir = directoryTemp,r = False)
       
   #close and remove the temporary bam
   os.unlink(tempBam)

 #if an in memory dir was generated remove that one
   if baminmemory == True:
      memory_fused.unmount()
      memory_fs_dir.removedir("ramdrive")


   #remove the temporary directory
   #os.removedirs(directory_temp)
   shutil.rmtree(directoryTemp,ignore_errors=True)


#generate the samstat output   
if __name__ == '__main__':
    alignAndSort()






