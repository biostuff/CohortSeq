# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 22:51:56 2015

@author: wolfthomas
"""


#define a function to annotate the vcfs with vep
def annotatevcf(vcfins,vcfouts,vep_loc,reference,directory_temp,ncores): 
   
   if isinstance(vcfins,list) == False:
       vcfins = [vcfins]
       vcfouts = [vcfouts]
        
   tempvcf = directory_temp + "/" + "temp.vcf"    
        
   for i in range(len(vcfins)):  
       variant_effect_predictor("--input_file",vcfins[i],"--output_file",tempvcf,"--offline","--cache","--vcf","--fasta",reference,"--pick_allele","--gencode_basic","--everything","--check_alleles","--failed","--force_overwrite","--fork",str(ncores))#"--no_stats")       
       
       vcf_reader = vcf.Reader(filename=tempvcf)
       vcf_writer = vcf.Writer(open(vcfouts[i], 'w'),vcf_reader)

      #if pick allele is not active pick only those annotation that fit the most likely allele
      #where is the allele position in the csq string ?
  
       for record in vcf_reader:
           max_pos = record.INFO['PP'].index(max(record.INFO['PP']))
           record.ID =   record.INFO["CSQ"][max_pos].split("|")[10]    
           vcf_writer.write_record(record)            
        
          
       vcf_writer.close()   

   os.unlink(tempvcf) 