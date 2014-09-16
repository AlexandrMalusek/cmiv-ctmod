SHELL     = sh
COMPILER  = gcc
#COMPILER  = icc

ifeq ($(COMPILER),gcc) 
# GNU C++ compiler settings
#CPATH     = /usr/local/gcc/3.3.2/bin/
CPATH     =
CC        = $(CPATH)gcc
LD        = $(CPATH)g++
CXX       = $(CPATH)g++
CPPFLAGS  = -I$(ROOTSYS)/include
CXXFLAGS  = -fPIC -pipe -I$(ROOTSYS)/include -Wall -g
# CXXFLAGS  = -fPIC -pipe -I$(ROOTSYS)/include -Wall -O2
LDFLAGS   = -g -shared -rdynamic
else
ifeq ($(COMPILER),icc)
# Intel C++ compiler settings
#CPATH     = /usr/local/intel/cc_80/bin/
CPATH     =
CC        = $(CPATH)icc
LD        = $(CPATH)icc
CXX       = $(CPATH)icc
CPPFLAGS  = -I$(ROOTSYS)/include
CXXFLAGS  = -g -fPIC -I$(ROOTSYS)/include
LDFLAGS   = -shared -g
endif
endif

RFLAGS    = -DUSE_ROOT

ROOTLIBS  := $(shell $(ROOTSYS)/bin/root-config --libs)

LIB      = libRCTmod.so
LIBDIR   = ../lib

all : $(LIB)

### General suffix rules
obj/%.d: %.cc
	$(SHELL) -ec '$(CXX) -MM $(CPPFLAGS) $(RFLAGS) $< \
	  | sed '\''s?\($*\)\.o[ :]*?\1.o $@ : ?g'\'' > $@; \
	  [ -s $@ ] || rm -f $@'

obj/%.o: %.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(RFLAGS) -o $@ $<

include Makefile_var

### Dependency files
DEPR = $(OBJ:.o=.d)

### Include dependency files
include $(DEPR)


### Libraries
$(LIB) : $(OBJ) obj/ctmodDict.o
	$(LD) $(LDFLAGS) -o $(LIB) $(OBJ) obj/ctmodDict.o
	chmod 775 $(LIB)
	rlibmap -f -o libRCTmod.rootmap -l libRCTmod.so -d libHist \
          libGraf libGraf3d libGpad libTree libRint libPostscript libMatrix \
          -c ctmodLinkDef.h

### Special dependencies
obj/ctmodDict.o : ctmodDict.cc

ctmodDict.cc : $(HDRS)
	@echo "Generating Dictionary ..."
	rootcint -v -f ctmodDict.cc -c -DUSE_ROOT $(HDRS) ctmodLinkDef.h

### Miscelaneous targets
.PHONY : clean install show tags

clean:
	rm -f $(LIB)
	rm -f $(OBJ)
	rm -f $(DEPR)
	rm -f ctmodDict.cc ctmodDict.h obj/ctmodDict.o

install: $(LIB)
	cp $(LIB) libRCTmod.rootmap $(LIBDIR)

show:
	@echo -n "ROOTLIBS: "
	@echo $(ROOTLIBS)

tags:
	ls *.{cc,h} | grep -v ctmodDict > file_names
	etags --language=c++ $$(cat file_names)

ebrowse:
	ls *.{cc,h} | grep -v ctmodDict > file_names
	ebrowse $$(cat file_names)
