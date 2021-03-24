# CAUS

Pipeline for Chromosome Assignment Using Synteny and it also fixes 3D-DNA's contig problems 

Say if you have two assemblies target.fasta and reference.fasta

You want to assign scaffolds in target.fasta to the reference assembly reference.fasta:

## How to run

        /tmp/casu/src/caus target.fasta reference.fasta assign.fasta > try.out

This gives you the assigned assembly file: assign.fasta 


### Note:                                                                         
	1. you need to give the full path of the caus  file.                     \
      	2. target.fasta and reference.fasta should be in your working directory. \
      	3. If you have a Segmentation fault (core dumped), please use try this   \

         /tmp/casu/caus -len 1000000 arget.fasta reference.fasta assign.fasta > try.out

## Download and Compile:

Requirements for compiling: gcc

	$ git clone https://github.com/wtsi-hpag/caus.git
	$ cd caus
	$ bash install.sh 

## Running with 3D-DNA 
3D-DNA is a popular HiC scaffolding pipeline. 
However, it makes so many breakpoints when examining HiC coverage profile.
After scaffolding with 3D-DNA, contigs are sometimes fragmented.
You can run caus to fix this problem. 

### Say you have contigs.fasta scaffolds-3D-DNA.fasta

	/tmp/casu/src/caus contigs.fasta scaffolds-3D-DNA.fasta scaffolds-caus3D.fasta > try.out

### Here                                                              
 	contigs.fasta is produced with long or short reads;       \ 
	scaffolds-3D-DNA.fasta is an output file from 3D-DNA;     \
	scaffolds-caus3D.fasta is a new scaffolding file.         \

More detailed information and examples can be found in the PPT file. 


Please contact Zemin Ning ( zn1@sanger.ac.uk ) for any further information. 



