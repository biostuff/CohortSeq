# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 22:51:27 2015

@author: wolfthomas
"""


#where are the bam files located you want to use for calling
bams_to_call_folder = output_folder_postprocessed

#get the names of the bam files for variant calling

def listdir(dirname, pattern="*"):
    return fnmatch.filter(os.listdir(dirname), pattern)

bams_to_call = listdir(dirname = bams_to_call_folder,pattern = "*.bam")


#define the class for each sample 0=normal 1 = tumor
origin = [1] * len(bams_to_call)

#define a patient idetifier important for tumor against matched normal
identity =  [1] * len(bams_to_call)

#which of the bams are normal
normal_pos =  [i for i, j in enumerate(origin) if j == 0]
#which of the bams are tumor
tumor_pos =  [i for i, j in enumerate(origin) if j == 1]


#generate a temp folder
directory_temp = tempfile.mkdtemp(dir = temp_folder)


#make sure that bams to call is a list
if isinstance(bams_to_call,list) == False:
    bams_to_call = [bams_to_call]



#get the names of the vcf files that should be output in the temp folder
vcf_output_temp = [directory_temp  + "/" + x.rstrip(".bam") + ".vcf" for x in bams_to_call] 





#transform the bed file into a platypus region string
def get_platypus_regions(bedfile,directory):
    #read the bed file used for calling
    region_information = BedTool(directory + "/" + bedfile)

    platypus_regions = ","
    platypus_regions = platypus_regions.join([x[0] + ":"+ x[1] + "-" + x[2] for x in region_information])

    return(platypus_regions)
    

merge_bedfile(bedfile,"amplified.bed",directory_temp,chrrem =  True)    

#get the region information that can be used for calling 
platypus_regions = get_platypus_regions("amplified.bed",directory_temp)


#run the variant calling
for i in range(len(bams_to_call)):
    if (origin[i]  == 1 or somatic == False):
        cur_bam_to_call = bams_to_call_folder + bams_to_call[i]
        if somatic == True and normal == True:
            normal_bam_pos = [j for j, k in enumerate(identity) if k == identity[i] and origin[i] == 0][0]
            matched_normal_bam =  bams_to_call_folder + bams_to_call[normal_bam_pos] 
            cur_bam_to_call = normal_bam_pos + "," + cur_bam_to_call
        platypus.callVariants(refFile = reference,bamFiles = cur_bam_to_call,output= vcf_output_temp[i],regions=platypus_regions,nCPU = ncores,filterDuplicates=0,genIndels = 1,minFlank = 3,mergeClusteredVariants = 0,maxReadLength=500,minGoodQualBases=70,countOnlyExactIndelMatches = 1,badReadsThreshold = 10)
     
      