
Description

CAUS is a pipeline for Chromosome Assignment Using Synteny 

Say if you have two assemblies target.fasta and reference.fasta

You want to assign scaffolds in target.fasta to the reference assembly reference.fasta:

You run

/tmp/casu/caus arget.fasta reference.fasta assign.fasta > try.out

This gives you the assigned assembly file: assign.fasta 


Note: 1. you need to give the full path of the caus  file.
      2. target.fasta and reference.fasta should be in your working directory.

Download and Compile:

Requirements for compiling: gcc

	$ git clone https://github.com/SangerHpag/scanPAV.git
	$ cd scanPAV 
	$ make 

Please contact Zemin Ning ( zn1@sanger.ac.uk ) for any further information. 



