INTRODUCTION
------------

CERC-Genomic-Medecine skill test2 (https://github.com/CERC-Genomic-Medicine/skills_test_2)

This command line tool writen in python is made to annotate a VCF file from a GENCODE GTF file. 



REQUIREMENTS
------------

Python 3.9.6

Required python packages are listed in: requirements.txt 
To install them, use : ```pip install -r requirements.txt```

List of packages used:
-pandas
-numpy
-re
-sys
-os
-datetime
-gzip
-shutil
-optparse
-collections
-cyvcf2


INSTALLATION
------------
Download the archive 

Uncompress the files : tar -zxvf CERC_vcf_annotation.tar.gz


USAGE
------------
Before launching any analysis make sure all the required packages are install
try : pip install -r requirements.txt


The program is a python3 script and can be launched in the shell.

e.g.: python3 CERC_vcf_annotation.py <parameters> <options>


OPTIONS
-----------

List of possible options;

-f / --vcf : Input VCF file

-g / --gtf :  Input Gencode GTF file

-i / --genes-in  :  Generate the list of overlapping genes

-a / --genes-around :  Generate the list of genes within +/-200000 base pairs

-n / --genes-nearest : Get the nearest gene 

-o / --output : Prefix of the outputs files. If not specified, output will be redirect as the input VCF file

-h / --help  : Show help message and exit


OUTPUTS	
-----------


InputFile.vcf.gz : Input VCF file with its INFO field modified with additional information required by the tool (GENES_IN, GENES_200KB, GENE_NEAREST)



EXAMPLES
-----------

> python3 CERC_vcf_annotation.py -f <VCF_filename>.vcf.gz -g <GENCODE_filename>.gtf.gz -o <output>

Will modify the header to add INFO field annotation but the INFO field will be untouched.
Results wil be written in a new VCF file named <output>.vcf.gz

> python3 CERC_vcf_annotation.py -f <VCF_filename>.vcf.gz -g <GENCODE_filename>.gtf.gz -i -o <output>

Will modify the header to add INFO field annotation and add the list of overlapping genes in the INFO field.
Results wil be written in a new VCF file named <output>.vcf.gz

> python3 CERC_vcf_annotation.py -f <VCF_filename>.vcf.gz -g <GENCODE_filename>.gtf.gz -a -o <output>

Will modify the header to add INFO field annotation and add the list of genes within +/- 200KB in the INFO field.
Results wil be written in a new VCF file named <output>.vcf.gz

> python3 CERC_vcf_annotation.py -f <VCF_filename>.vcf.gz -g <GENCODE_filename>.gtf.gz -n -o <output>

Will modify the header to add INFO field annotation and add the nearest genes in the INFO field.
Results wil be written in a new VCF file named <output>.vcf.gz

> python3 CERC_vcf_annotation.py -f <VCF_filename>.vcf.gz -g <GENCODE_filename>.gtf.gz -i -a -n -o <output>

Will modify the header to add INFO field annotation and add the list of overlapping genes, the list of genes within +/- 200KB and the nearest genes in the INFO field.
Results wil be written in a new VCF file named <output>.vcf.gz


> python3 CERC_vcf_annotation.py -f <VCF_filename>.vcf.gz -g <GENCODE_filename>.gtf.gz -i -a -n 

Will modify the header to add INFO field annotation and add the list of overlapping genes, the list of genes within +/- 200KB and the nearest genes in the INFO field.
Results wil be written in the input VCF file <VCF_filename>.vcf.gz


AUTHOR
-----------
PELLETIER Justin (https://mhi-omics.org/people/justin-pelletier/)
email: justin.pelletier@umontreal.ca

