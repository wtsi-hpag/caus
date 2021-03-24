# CAUS

Pipeline for Chromosome Assignment Using Synteny 

Say if you have two assemblies target.fasta and reference.fasta

You want to assign scaffolds in target.fasta to the reference assembly reference.fasta:

## How to run

        /tmp/casu/src/caus target.fasta reference.fasta assign.fasta > try.out

This gives you the assigned assembly file: assign.fasta 


Note: 1. you need to give the full path of the caus  file.
      2. target.fasta and reference.fasta should be in your working directory.
      3. If you have a Segmentation fault (core dumped), please use try this

         /tmp/casu/caus -len 1000000 arget.fasta reference.fasta assign.fasta > try.out

##Download and Compile:

Requirements for compiling: gcc

	$ git clone https://github.com/wtsi-hpag/caus.git
	$ cd caus
	$ bash install.sh 

## Other applications
3D-DNA is a popular HiC scaffolding pipeline. 
However, it makes so many breakpoints when examining HiC coverage profile.
After scaffolding with 3D-DNA, contigs are sometimes fragmented.
Run caus to fix this problem, say you have contigs.fasta scaffolds-3D-DNA.fasta

	/tmp/casu/src/caus contigs.fasta scaffolds-3D-DNA.fasta scaffolds-caus3D.fasta > try.out

More detailed information and examples can be found in the PPT file. 


Please contact Zemin Ning ( zn1@sanger.ac.uk ) for any further information. 



