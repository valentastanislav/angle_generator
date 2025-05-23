ROOTLIB := $(shell root-config --glibs)
ROOTINC := $(shell root-config --cflags)
ROOTVERSION := $(shell root-config --version | tr -c -d [:digit:])
ROOTCHECK := $(shell if [ $(ROOTVERSION) -gt 60000 ]; then echo true; else echo false; fi )

BIN		= ./
CC		= g++
CCFLAG	= -O -Wall -fPIC -lstdc++

dicebox_to_GEANT:		RTG_energy.o Fk.o dicebox_to_GEANT.o
ifeq ($(strip $(ROOTCHECK)),true)
	$(info have good ROOT version [${ROOTVERSION}])
else
	$(error have bad ROOT version [${ROOTVERSION}])
endif
	$(info $$ROOTLIB is [${ROOTLIB}])
	$(info $$ROOTINC is [${ROOTINC}])
	$(CC) $(CCFLAG) -o $@ $^ $(ROOTLIB)

dicebox_to_GEANT.cc:	Fk.h
dicebox_to_GEANT.o:		dicebox_to_GEANT.cc
	$(CC) $(ROOTINC) -c $^

Fk.cc:					Fk.h
Fk.o:					Fk.cc
	$(CC) -c $^

RTG_energy.cc:			RTG_energy.h
RTG_energy.o:			RTG_energy.cc
	$(CC) -c $^

clean:
	rm -f *.o dicebox_to_GEANT