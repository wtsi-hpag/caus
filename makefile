# Makefile for caus 
CC= gcc
CFLAGS= -O2
LFLAGS= -lm

CAUSPOUT= fast.o caus_seqout.o 
CAUSFAST= fast.o caus_fasta.o
CAUSSHED= fast.o caus_shred.o
CAUSPASA= caus_assign.o 
CAUSSMT1= caus_smalt1.o 
CAUSSMT2= caus_smalt2.o 
CAUSSMT3= caus_smalt3.o 
CAUSPIPE= caus.o 

all  : caus_seqout caus_fasta caus_shred caus_assign caus_smalt1 caus_smalt2 caus_smalt3 caus

caus_seqout: makefile $(CAUSPOUT)
	$(CC) $(CFLAGS) -o $@ $(CAUSPOUT) $(LFLAGS) 
	chmod o-r caus_seqout
	cp caus_seqout caus-bin
       
caus_fasta: makefile $(CAUSFAST)
	$(CC) $(CFLAGS) -o $@ $(CAUSFAST) $(LFLAGS) 
	chmod o-r caus_fasta 
	cp caus_fasta caus-bin 

caus_shred: makefile $(CAUSSHED)
	$(CC) $(CFLAGS) -o $@ $(CAUSSHED) $(LFLAGS) 
	chmod o-r caus_shred 
	cp caus_shred caus-bin 

caus_assign: makefile $(CAUSPASA)
	$(CC) $(CFLAGS) -o $@ $(CAUSPASA) $(LFLAGS) 
	chmod o-r caus_assign 
	cp caus_assign caus-bin

caus_smalt1: makefile $(CAUSSMT1)
	$(CC) $(CFLAGS) -o $@ $(CAUSSMT1) $(LFLAGS) 
	chmod o-r caus_smalt1 
	cp caus_smalt1 caus-bin

caus_smalt2: makefile $(CAUSSMT2)
	$(CC) $(CFLAGS) -o $@ $(CAUSSMT2) $(LFLAGS) 
	chmod o-r caus_smalt2 
	cp caus_smalt2 caus-bin

caus_smalt3: makefile $(CAUSSMT3)
	$(CC) $(CFLAGS) -o $@ $(CAUSSMT3) $(LFLAGS) 
	chmod o-r caus_smalt3 
	cp caus_smalt3 caus-bin

caus: makefile $(CAUSPIPE)
	$(CC) $(CFLAGS) -o $@ $(CAUSPIPE) $(LFLAGS) 
	chmod o-r caus
	cp caus caus-bin


