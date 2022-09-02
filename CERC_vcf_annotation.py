#!/usr/bin/env python3  
# -*- coding: utf-8 -*-
#--------------------------------------------------------
# Author: Justin Pelletier
# Creation date: 08/08/2022
# Version 1.0
#--------------------------------------------------------
#Imports
import pandas as pd 
import numpy as np 
import re
import sys
import os
import os.path
import datetime
import gzip
import shutil
from optparse import OptionParser
from collections import OrderedDict
import cyvcf2
from cyvcf2 import VCF, Writer




    
def readGTF(infile):
    """
    Function from the source code of gtfparse (https://github.com/openvax/gtfparse)
    Reads a GTF file and labels the respective columns in agreement with GTF file standards:
    'seqname','source','feature','start','end','score','strand','frame','attribute'.
    :param infile: path/to/file.gtf
    :returns: a Pandas dataframe of the respective GTF
    """
    df = pd.read_table(infile, sep='\t', comment="#", header=None, dtype=str)
    df.columns = ['seqname','source','feature','start','end','score','strand','frame','attribute']
    return df


def retrieve_GTF_field(field,gtf):
    """
    Function from the source code of gtfparse (https://github.com/openvax/gtfparse)
    Returns a field of choice from the attribute column of the GTF
    :param field: field to be retrieved
    :returns: a Pandas dataframe with one columns containing the field of choice
    """
    inGTF = gtf.copy()
    def splits(x):
        l = x.split(";")
        l = [ s.split(" ") for s in l]
        res = np.nan
        for s in l:
            if field in s:
                if '"' in s[-1]:
                    res = s[-1][1:-1]
                else:
                    res = s[-1]
        return res
    inGTF[field] = inGTF['attribute'].apply(lambda x: splits(x))
    return inGTF[[field]]
	

	
def genes_in(chrom, pos, GTF):
    """
    Returns a list of overlapping gene_id from the GTF file from the variant
    :param field: chromosome, position, GTF dataframe
    :returns: ordered list of gene_id from the gencode file overlapping the variant
    """
    list_genes = list(dict.fromkeys(GTF[(GTF['end'].astype(int) >= pos) & (GTF['start'].astype(int) < pos) ]['gene_id']))
    #list_genes = list(dict.fromkeys(GTF[ (GTF['seqname'] == chrom) & (GTF['end'].astype(int) >= pos) & (GTF['start'].astype(int) < pos) ]['gene_id']))
    return list_genes
    	

def genes_around(chrom, pos, GTF):
    """
    Returns a list of surrounding gene_id from the GTF file of genes within +/- 200KB from the variant and its associated gtf file
    :param field: chromosome, position, GTF dataframe
    :returns: ordered list of gene_id from the gencode file of genes within +/- 200KB from the variant and a gtf file of genes within +/- 200KB 
    """
    new_gtf = GTF[(GTF['end'].astype(int) <= pos+200000) & (GTF['start'].astype(int) > pos-200000) ]
    list_genes = list(dict.fromkeys(new_gtf['gene_id']))
    #list_genes = list(dict.fromkeys(GTF[ (GTF['seqname'] == chrom) & (GTF['end'].astype(int) <= pos+200000) & (GTF['start'].astype(int) > pos-200000) ]['gene_id']))
    return list_genes, new_gtf
	
	
	
	
def gene_nearest(chrom, pos, GTF):
    """
    Returns the nearest gene's gene_id from the GTF from the variant
    :param field: chromosome, position, GTF dataframe
    :returns: list of one gene_id from the gencode file
    """
    # Get the start of the closest gene from the variant
    ordered_distance_start=GTF.iloc[(GTF['start'].astype(int) - pos).abs().argsort(),:]
    good_chr_start=ordered_distance_start[(ordered_distance_start['seqname'] == chrom)]
    list_start_pos=list(dict.fromkeys(good_chr_start.iloc[:1]['start']))
    list_start_id=list(dict.fromkeys(good_chr_start.iloc[:1]['gene_id']))
    # Get the end of the closest gene from the variant
    ordered_distance_end=GTF.iloc[(GTF['end'].astype(int) - pos).abs().argsort(),:]
    good_chr_end=ordered_distance_end[(ordered_distance_end['seqname'] == chrom)]
    list_end_pos=list(dict.fromkeys(good_chr_end.iloc[:1]['end']))
    list_end_id=list(dict.fromkeys(good_chr_end.iloc[:1]['gene_id']))
    
    # Compare the proximity of the closest start and end of the identified genes to get the closests
    if(abs(int(list_end_pos[0]) - pos) < abs(int(list_start_pos[0]) - pos)):
        list_genes = list_end_id
    else:
        list_genes = list_start_id
    return list_genes
 



# MAIN
def main():
	
    """
    First step: initialize the option and parse them
    """
        
    # Parser to get options from the command line
    parser = OptionParser()
    
    # Define the available options, their action, destination and help message
    parser.add_option("-f", "--vcf", action="store", type="string",  dest="vcf_file", help="Input VCF file")
    parser.add_option("-g", "--gtf", action="store", type="string",  dest="gencode_file", help="Input Gencode GTF file")
    parser.add_option("-i", "--genes-in", action="store_true", dest="genes_in", default=False, help="Generate the list of overlapping genes")
    parser.add_option("-a", "--genes-around", action="store_true", dest="genes_around", default=False, help="Generate the list of genes within +/-200 000 base pairs")
    parser.add_option("-n", "--genes-nearest", action="store_true", dest="gene_nearest", default=False, help="Get the nearest gene")
    parser.add_option("-o", "--output", action="store", type="string",  dest="output_file", default="output", help="Prefix of the outputs files")
    

    # Start timer for input files reading time
    begin_time = datetime.datetime.now()

    # Retrieve options from the parser
    (options, args) = parser.parse_args()
    option_dict = vars(options)
    
    # Parse and modify the gtf
    gencode_file = option_dict['gencode_file']
    gtf_file = readGTF(gencode_file)
    # Adds a column for gene_id in the GTF file
    gtf_file["gene_id"] = retrieve_GTF_field("gene_id",gtf_file)
    
	
    # Stores the output filename
    # If the option is called store the new output name, else, use the input VCF filename
    if(option_dict['output_file'] != "output"):
        output_file = option_dict['output_file'] + ".vcf.gz.tmp"
        output_file_print = option_dict['output_file'] + ".vcf.gz"
    else:
        output_file = option_dict['vcf_file'] + ".tmp"
        output_file_print = option_dict['vcf_file']	

    # Read the GTF file	
    gtf_file = readGTF(gencode_file)
    # Adds a column for gene_id in the GTF file
    gtf_file["gene_id"] = retrieve_GTF_field("gene_id",gtf_file)
    
    # Get the options 
    f_genes_in = option_dict['genes_in']
    f_genes_around = option_dict['genes_around']
    f_gene_nearest = option_dict['gene_nearest']

    # Reads the vcf file
    vcf_file=option_dict['vcf_file']
    vcf= VCF(vcf_file)

    # Print the input and output files
    print("Gencode file : " + gencode_file)
    print("VCF file: " + vcf_file)
    print("Output file : " + output_file_print)

    # End timer of file reading
    input_time=(datetime.datetime.now() - begin_time)
    print("Input reading time: " + str(input_time))

    # Start timer for execution time
    begin_time = datetime.datetime.now()
    
    # Print out the selected options
    print("Selected options; ")
    if(f_genes_in == True):
        print(" Overlapping enes ")
    if(f_genes_around == True):
        print(" Genes within +/- 200kb ")
    if(f_gene_nearest == True):
        print(" Nearest gene ")


    """
    Second step: Create INFO field for the options in the vcf file
                 Launch the selected options to fill the INFO fields
    """

    # Add an info field for each option selected
    vcf.add_info_to_header({'ID': 'GENES_IN', 'Description': 'overlapping gene', 'Type':'Character', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'GENES_200KB', 'Description': 'gene within 200KB', 'Type':'Character', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'GENE_NEAREST', 'Description': 'closest gene', 'Type':'Character', 'Number': '1'})

	
    # Open a Writer obejct to scan the vcf file and output the modified one
    w = Writer(output_file, vcf)
    
    # Loops through the vcf file and execute the selected options for each variant
    for variant in vcf:
        
        # Subset the GTF for the appropriate chromosome
        gtf=gtf_file[(gtf_file['seqname'] == variant.CHROM)]

        # Set the default value of genes_in_empty and genes_around_empty as true  
        genes_in_empty = True
        genes_around_empty = True
        # Launch the option to get surrounding genes within 200KB of the variant
        if(f_genes_around == True):
            genes_200kb, purged_gtf = genes_around(variant.CHROM, variant.POS, gtf)
            # Changes the value of genes_around_empty 
            if (len(genes_200kb) != 0):
                genes_around_empty = False
        
        # Launch the function to get overlaping genes
        if(f_genes_in == True):
            # If the genes_around funtion is called, give the genes_in function a reduced size gtf of 400kb genes to look at
            if(f_genes_around == True and genes_around_empty == False):
                # Launch the function only if there are genes around
                if(genes_around_empty == False):
                    genes = genes_in(variant.CHROM, variant.POS, purged_gtf)
                else: 
                    genes = []
            else:
                genes = genes_in(variant.CHROM, variant.POS, gtf)
            # Writes it in the INFO field if the list overlapping genes and set the genes_in_empty variable
            if (len(genes) != 0):
                variant.INFO["GENES_IN"] = ",".join(genes)
                genes_in_empty=False
            else:
                variant.INFO["GENES_IN"] = "."
        

        # Writes it in the INFO field if the list of genes within +/- 200KB 
        if(f_genes_around == True):
            if (len(genes_200kb) != 0):
                variant.INFO["GENES_200KB"] = ",".join(genes_200kb)
            else:
                variant.INFO["GENES_200KB"] = "."


        # Launch the option to get nearest gene if the genes has no overlapping gene
        if(f_gene_nearest == True and genes_in_empty == True):
            # When genes around isn't empty, pass the resulting gtf to the function to reduce execution time
            if(genes_around_empty == False):
                genes = gene_nearest(variant.CHROM, variant.POS, purged_gtf)
            else:
                genes = gene_nearest(variant.CHROM, variant.POS, gtf)
            # Writes it in the INFO field if the nearest gene 
            if (len(genes) != 0):
                variant.INFO["GENE_NEAREST"] = ",".join(genes)
            else:
                variant.INFO["GENE_NEAREST"] = "."
        # Launch the option to get nearest gene if the genes has overlapping gene, skip the call of the gene_nearest function 
        if(f_gene_nearest == True and genes_in_empty == False):
            variant.INFO["GENE_NEAREST"] = "."
	
        #write the variant with its new INFO field in the output file
        w.write_record(variant)
		
    
    """
    Last step: Write and close the modified VCF file then compress it
    """
    #closes the vcf file	
    w.close(); vcf.close()
    
    #Compress output file
    with open(output_file, 'rb') as f_in, gzip.open(output_file_print, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    #remove uncompressed file
    os.remove(output_file)


    # Compute the total execution time of the script
    total_time=(datetime.datetime.now() - begin_time)
    print("Execution time: " + str(total_time))
    		
    	
    	
if __name__ == "__main__":
    main()

# Run to get the required packages
# pip install -r requirements.txt

