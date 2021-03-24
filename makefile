# Makefile for caus
CC= gcc
CFLAGS= -O2
LFLAGS= -lm

SOURCES= caus_seqout.c caus_fasta.c fastq2fasta.c caus_shred.c caus_assign.c caus_smalt0.c caus_smalt1.c caus_smalt2.c caus_smalt3.c caus_clean.c caus.c

OBJS = $(patsubst %.c,%.o,$(SOURCES)) fast.o
EXECS = $(patsubst %.c,%,$(SOURCES))
EXECS_BIN = $(patsubst %.c,caus-bin/%,$(SOURCES))
COMPILE = $(CC) $(CFLAGS)


all:  cleanall iprint $(OBJS) executables clean oprint

executables:
	for exe in $(EXECS);  do $(COMPILE) -o $$exe $$exe.o fast.o $(LFLAGS); cp $$exe caus-bin/.; done

%.o: %.c fasta.h
	$(CC) $(CFLAGS)  -c $< -o $@

iprint:
	@echo '+++ Compiling All ... '

oprint:
	@echo 'All Done '


clean:
	@echo '+++ Cleaning Up ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS)
	@cp caus-bin/caus .

cleanall:
	@echo '+++ Cleaning All ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS) fast.o
	@rm -f $(EXECS_BIN)
