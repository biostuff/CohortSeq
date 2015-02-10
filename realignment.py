# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 15:28:47 2015

@author: wolfthomas
"""

region_information = BedTool(directory_temp + "/" + "amplified.bed")
import numpy
realign_chromosomes = numpy.unique([x[0]  for x in region_information])



#how are the splitted bed files called
#realign_regions = listdir(bed_split_dir, pattern="*")



for i in range(len(bams_to_realign)):
    
    #where are the splitted bams to be located
    #bam_split_dir =  tempfile.mkdtemp(dir = directory_temp) 
    
    #where should the realigned bams be stored
    #pyro_split_dir = tempfile.mkdtemp(dir = directory_temp) + "/"
    
    final_realigned_bam_dir =  tempfile.mkdtemp(dir = directory_temp)  
        
    #cur_chromosome =    realign_chromosomes[0]     
        
    bam_to_realign = bams_to_realign[i]        
      
    #where should the realigned bams be stored
    pyro_split_dir = tempfile.mkdtemp(dir = directory_temp)     
      
    #final_realigned_bam_dir =  tempfile.mkdtemp(dir = directory_temp)           
    
    
    def pyroMulti(input_folder,bam_to_realign,realign_chromosome,output_folder,maxmem,ncores,pyromodel,directory_temp):       
                   
        #where are the splitted bams to be located
        bam_split_dir =  tempfile.mkdtemp(dir = directory_temp) 
    
           
        #realign_chromosome = realign_chromosomes[0] 
    
        cur_bam = bam_split_dir + "/" +  str(realign_chromosome) + "_"  + "filtered.bam"
        
        sambamba.slice("-o",cur_bam,input_folder + bam_to_realign,realign_chromosome)        
        #sambamba.index(cur_bam)
        cur_tempbampyro = output_folder + "/" + str(realign_chromosome)  + "_" + "temppyro.bam"           
        cur_regions = realign_chromosome + ":0"
        pyroTools.GraphReAlign("-realign","-bam",cur_bam,"-ref",reference,"-c","-out",cur_tempbampyro,"-config",pyromodel,"-large-homopolymer","-update-graph-by-read","-amplicon","-keep-duplicate","-min-read-len",str(70),"-region",cur_regions,"-verbosity",5,_ok_code=[0,1,3,5])       
        samtools.sort("-o",cur_tempbampyro,"-O","bam","-m",str(int(math.floor(maxmem/ncores)))+"G","-T",directory_temp + "/" + "sorttempbam/",cur_tempbampyro)                
    
    
    #pyroMulti(output_folder,bams_to_realign[1],realign_chromosomes[0],pyro_split_dir,maxmem,ncores,pyromodel)      
    
    Parallel(n_jobs=ncores)(delayed(pyroMulti)(output_folder,bams_to_realign[i],j,pyro_split_dir,maxmem,ncores,pyromodel,directory_temp) for j in realign_chromosomes)    
        
    tempbampyro = final_realigned_bam_dir  + "/" + "temppyro.bam"
 
    #get all the realigned bams 
    realigned_bams = listdir(pyro_split_dir, pattern="*.bam")     
    realigned_bams = [pyro_split_dir + "/" + x for x in realigned_bams]     
     
    
    picard_input = " ".join(["INPUT=" + x for x in realigned_bams])  
    args = [picard_string,"MergeSamFiles",picard_input,"OUTPUT=",tempbampyro,"TMP_DIR=",directory_temp,"VALIDATION_STRINGENCY=SILENT"," CREATE_INDEX=","true","USE_THREADING=","true"]
    subprocess.call([" ".join(args)],shell = True)
    
    #samtools.sort(tempbampyro,tempbampyro,"-@",str(ncores),"-m",str(int(math.floor(maxmem/ncores)))+"G")        
    
    
    #sambamba.index(tempbampyro,"-t",ncores)
    
    if paired == True:
        samtools.sort("-n","-o",tempbampyro,"-O","bam","-@",str(ncores),"-m",str(int(math.floor(maxmem/ncores)))+"G","-T",directory_temp + "/" + "sorttempbam/",tempbampyro)  
        args = [picard_string,"FixMateInformation","INPUT=",tempbampyro,"TMP_DIR=",directory_temp,"VALIDATION_STRINGENCY=SILENT"," CREATE_INDEX=","false","ASSUME_SORTED=true"]
        subprocess.call([" ".join(args)],shell = True)
        samtools.sort("-o",tempbampyro,"-O","bam","-@",str(ncores),"-m",str(int(math.floor(maxmem/ncores)))+"G","-T",directory_temp + "/" + "sorttempbam/",tempbampyro)  
   
    sambamba.index(tempbampyro,"-t",ncores)   
   
    #left align the indels          
    args = [gatk_string," -T LeftAlignIndels "," -I ",tempbampyro, " -R ",reference," -o ",output_folder_postprocessed + bams_to_realign[i]]
    subprocess.call(["".join(args)],shell = True)
    sambamba.index(output_folder_postprocessed + bams_to_realign[i],"-t",ncores)
    
    os.unlink(tempbampyro)
os.removedirs(directory_temp)
    