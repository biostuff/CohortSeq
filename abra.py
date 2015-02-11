# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 22:38:24 2015

@author: wolfthomas
"""


def listdir(dirname, pattern="*"):
    return fnmatch.filter(os.listdir(dirname), pattern)


bams_to_realign = listdir(dirname = output_folder,pattern = "*.bam")

directory_temp = tempfile.mkdtemp(dir = temp_folder)     

merge_bedfile(bedfile,"amplified.bed",directory_temp,chrrem =  True)  


#identity =  [1,2,3]
identity = range(len(bams_to_realign))

identity_levels = unique(identity)
for i in range(identity_levels):
    directory_temp_abra = tempfile.mkdtemp(dir = directory_temp)     
    
    cur_bam_pos = [j for j, k in enumerate(identity) if k == identity_levels[i]]
    cur_bams = [bams_to_realign[j] for j in cur_bam_pos]
    output_bams = [directory_temp + cur_bams[j] for j in range(len(cur_bams))]
    output_bams_realigned = [output_folder_postprocessed + cur_bams[j] for j in range(len(cur_bams))]
    cur_bams = [output_folder + j for j in cur_bams]   
   
   
    cur_bams_string = ",".join(cur_bams)     
    output_bams_string = ",".join(output_bams)    
   
    args = [abra_string,"--in",cur_bams_string,"--out", output_bams_string,"--ref",reference,"--targets",directory_temp + "/" + "amplified.bed","--threads",str(ncores),"--working",directory_temp_abra + "test","--mad ",str(1000000)]
  
    subprocess.call([" ".join(args)],shell = True) 
    for k in range(len(cur_bam_pos)):
        if paired == True:
            samtools.sort("-n","-o",output_bams[k],"-O","bam","-@",str(ncores),"-m",str(int(math.floor(maxmem/ncores)))+"G","-T",directory_temp + "/" + "sorttempbam/",output_bams[k])  
            args = [picard_string,"FixMateInformation","INPUT=",output_bams[k],"TMP_DIR=",directory_temp,"VALIDATION_STRINGENCY=SILENT","ASSUME_SORTED=true"]
            subprocess.call([" ".join(args)],shell = True)  
        sambamba.sort(output_bams[k],o =  output_bams_realigned[k],t = str(ncores),tmpdir = directory_temp + "/" + "sorttempbam",m = str(maxmem) + "GB")   