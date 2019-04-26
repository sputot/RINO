
PREFIX=aaf_lib
LIB=libaaf.a
LIBNAME   = $(LIB:%.a=%)
SHAREDLIB = $(LIBNAME).so
DEBUG = 

MAKE = make -j
MAKEFILE = Makefile

CCC = g++
	CCFLAGS += -DLIN64 -Wall -DUSE_AAF_EXTENSIONS -Wno-deprecated -O3 #-g #-O3
CCFLAGS += -Weffc++ $(DEBUG)

CXXFLAGS += -I.. -I$(FILIBHOME)
LDFLAGS  += -L.

AR = ar
ARFLAGS = ruv
RANLIB = ranlib

RM = rm -f 

DOXYGEN = doxygen

SOURCES.c =	aa_rounding.c 	aa_interval.c 	aa_aaftrigo.c 	aa_aafapprox.c \
		aa_aafarithm.c 	aa_aafcommon.c 	aa_exceptions.c	aa_aafspecial.c aa_aafinnerflow.c

SOURCES.h =	aa.h aa_aaf.h 	aa_interval.h 	aa_rounding.h	aa_exceptions.h 

SOURCES = $(SOURCES.h)  $(SOURCES.c)

OBJECTS = $(SOURCES.c:%.c=%.o)

OBJECTS_SO = $(SOURCES.c:%.c=%.so.o)

default: all

all:	 $(LIB) $(SHAREDLIB)

static:  $(LIB)

shared:  $(SHAREDLIB)

debug:
	@$(MAKE) -f $(MAKEFILE) cleaner
	@$(MAKE) -f $(MAKEFILE) DEBUG='-g'

$(LIB): $(SOURCES) $(OBJECTS)
	@echo "\n - Compiling $(LIB) -\n"
	@rm -f $(LIB)
	$(AR) $(ARFLAGS) $(LIB) $(OBJECTS)
	$(RANLIB) $(LIB) 

# %.o: %.c $(SOURCES.h)
# 	$(CCC) $(CCFLAGS) $(CXXFLAGS) -c $(@:%.o=%.c) -o $@ 
%.o: %.c $(SOURCES.h)
	$(CCC) $(CCFLAGS) $(CXXFLAGS) -g -O -c $(@:%.o=%.c) -o $@                                        

$(SHAREDLIB): $(OBJECTS_SO)
	g++ -shared -Wall -o $@  $(OBJECTS_SO)

%.so.o: %.c $(SOURCES.h)
	$(CCC) $(CCFLAGS) $(CXXFLAGS) -fPIC -c $(@:%.so.o=%.c) -o $@                                             

clean::
	$(RM) $(OBJECTS) $(OBJECTS_SO)

cleaner:: clean
	$(RM) $(LIB) $(SHAREDLIB) config.status
	$(RM) -r doxygen

doxygen: $(SOURCES) doxygen.*
	$(DOXYGEN) doxygen.cfg
	@echo "doxygen done."

