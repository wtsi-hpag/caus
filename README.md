
Description

CAUS is a pipeline for Chromosome Assignment Using Synteny 

Say if you have two assemblies target.fasta and reference.fasta

You want to assign scaffolds in target.fasta to the reference assembly reference.fasta:

You run

/tmp/casu/caus arget.fasta reference.fasta assign.fasta > try.out

This gives you the assigned assembly file: assign.fasta 


Note: 1. you need to give the full path of the caus  file.
      2. target.fasta and reference.fasta should be in your working directory.
      3. If you have a Segmentation fault (core dumped), please use try this

         /tmp/casu/caus -len 1000000 arget.fasta reference.fasta assign.fasta > try.out

Download and Compile:

Requirements for compiling: gcc

	$ git clone https://github.com/wtsi-hpag/caus.git
	$ cd caus
	$ make 

Please contact Zemin Ning ( zn1@sanger.ac.uk ) for any further information. 



