CXX = g++ 

MAX_ORDER = 40

LIBS = -laaf -lprim -lgsl -llapack -lblas -lstdc++ -lyaml-cpp

CURRENT_DIR = $(shell pwd)

CXXFLAGS = -ggdb -Wall -frounding-math -DMAXORDER=$(MAX_ORDER) -I. -I$(FILIBHOME) -I$(FADBADHOME) -I/usr/local/include \
	 -I$(CURRENT_DIR)/aaflib-0.1 -fpermissive -std=c++11

LDFLAGS  +=  -L/usr/local/lib -L$(CURRENT_DIR)/aaflib-0.1 -L$(FILIBHOME)/libprim/.libs/

SOURCES_utils = utils.cpp matrix.cpp network_handler.cpp inner.cpp ode_def.cpp ode_integr.cpp dde_integr.cpp discrete_system.cpp

SOURCES.h = filib_interval.h fadbad_aa.h utils.h network_handler.h matrix.h inner.h ode_def.h ode_integr.h dde_integr.h  discrete_system.h

SOURCES = $(SOURCES_utils) main.cpp

OBJECTS_utils = $(SOURCES_utils:%.cpp=%.o)

OBJECTS = $(SOURCES:%.cpp=%.o)

all: $(OBJECTS) main

main : $(SOURCES) $(OBJECTS) 
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS_utils) $@.o $(LIBS)

clean:
	-rm *.o main

.cpp.o:
	$(CXX) $(CXXFLAGS)  -c $<

%.o: %.cpp $(SOURCES.h)
	$(CXX) $(CXXFLAGS) -c $(@:%.o=%.cpp) -o $@ 

