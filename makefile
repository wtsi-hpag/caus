# Makefile for easyChain 
CC= gcc
C2C= g++
CFLAGS= -O2
CFLAG3= -O3
CFLASS= -static 
LFLAGS= -lm
PFLAGS= -pthread

SOURCES=easyChain_fasta.c easyChain_shred.c easyChain.c


OBJS = $(patsubst %.c,%.o,$(SOURCES)) fast.o
EXECS = $(patsubst %.c,%,$(SOURCES))
EXECS_BIN = $(patsubst %.c,easyChain-bin/%,$(SOURCES))
COMPILE = $(CC) $(CFLAGS) $(PFLAGS)

all:  cleanall iprint $(OBJS) executables clean oprint

executables:
	for exe in $(EXECS);  do $(COMPILE) -o $$exe $$exe.o fast.o; cp $$exe easyChain-bin/.; done

%.o: %.c fasta.h
	$(CC) $(CFLAGS) $(PFLAGS)  -c $< -o $@

iprint:
	@echo '+++ Compiling All ... '

oprint:
	@echo 'All Done '

clean: 
	@echo '+++ Cleaning Up ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS)
	@cp easyChain-bin/easyChain .

	$(C2C) $(PFLAGS) $(CFLAG3) $(CFLASS) checkError.c -o checkError
cleanall: 
	@echo '+++ Cleaning All ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS) fast.o
	@rm -f $(EXECS_BIN)
