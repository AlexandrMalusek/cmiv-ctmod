SHELL     = sh
COMPILER  = gcc

ifeq ($(COMPILER),gcc) 
# GNU C++ compiler settings
CPATH     =
CC        = $(CPATH)gcc
LD        = $(CPATH)g++
CXX       = $(CPATH)g++
CPPFLAGS  = -I$(ROOTSYS)/include
CXXFLAGS  = -g -Wall -fPIC -I$(ROOTSYS)/include \
            -I$(CTMODSYS)/src -DUSE_ROOT -L$(ROOTSYS)/lib
LDFLAGS   = -g -shared -rdynamic
else
ifeq ($(COMPILER),icc)
# Intel C++ compiler settings
CPATH     =
CC        = $(CPATH)icc
LD        = $(CPATH)icc
CXX       = $(CPATH)icc
CPPFLAGS  = -I$(ROOTSYS)/include
CXXFLAGS  = -g -O -fPIC -I$(ROOTSYS)/include
LDFLAGS   = -shared -g
endif
endif

PFLAGS    = 
RFLAGS    = -DUSE_ROOT

include $(CTMODSYS)/src/Makefile_var

ctmod : ctmod.o getMaterialManager.o getGeometry.o getTomograph.o TInputFile.o $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ ctmod.o getMaterialManager.o getGeometry.o \
	getTomograph.o TInputFile.o $(OBJ) obj/ctmodDict.o \
	`root-config --cflags --libs`

getMaterialManager.o : getMaterialManager.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

getGeometry.o : getGeometry.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

getTomograph.o : getTomograph.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

TInputFile.o : TInputFile.cc TInputFile.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY : clean symlinks install

clean :
	rm -f getMaterialManager.o getGeometry.o getTomograph.o ctmod.o \
	TInputFile.o ctmod

symlinks :
	ln -fs $(CTMODSYS)/src/obj obj
	ln -fs $(CTMODSYS)/material material
	ln -fs $(CTMODSYS)/detector detector
	ln -fs $(CTMODSYS)/spectra spectra

install :
	cp ctmod ../bin
