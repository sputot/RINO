CXX = g++ -g 

MAX_ORDER = 40


#ONNX_active = TRUE
ONNX_active = FALSE

CURRENT_DIR = $(shell pwd)


LIBS = -laaf -lprim -lgsl -llapack -lblas -lstdc++ -lyaml-cpp


CXXFLAGS = -ggdb -Wall -frounding-math -DMAXORDER=$(MAX_ORDER) -DUSE_AAF_EXTENSIONS -I. -I$(FILIBHOME) -I$(FADBADHOME) -I/usr/local/include \
	 -I$(CURRENT_DIR)/aaflib-0.1 -fpermissive -std=c++11

LDFLAGS  +=  -L/usr/local/lib -L$(CURRENT_DIR)/aaflib-0.1 -L$(FILIBHOME)/libprim/.libs/

SOURCES_utils = utils.cpp matrix.cpp network_handler.cpp inner.cpp ode_def.cpp ode_integr.cpp dde_integr.cpp discrete_system.cpp

SOURCES.h = filib_interval.h fadbad_aa.h utils.h network_handler.h matrix.h inner.h ode_def.h ode_integr.h dde_integr.h  discrete_system.h

SOURCES = $(SOURCES_utils) main.cpp

OBJECTS_utils = $(SOURCES_utils:%.cpp=%.o)

ifeq ($(ONNX_active),TRUE)
	SHERLOCK_DIR = $(CURRENT_DIR)/sherlock_2_reduct/
	CXXFLAGS += -I$(SHERLOCK_DIR)/src -I$(SHERLOCK_DIR)/eigen_file
	LIBS += -lprotobuf
	SOURCES.h += sherlock.h
	OBJS_SHERLOCK = $(SHERLOCK_DIR)/src/onnx.pb.o $(SHERLOCK_DIR)/src/sherlock.o $(SHERLOCK_DIR)/src/network_computation.o \
$(SHERLOCK_DIR)/src/configuration.o $(SHERLOCK_DIR)/src/nodes.o $(SHERLOCK_DIR)/src/computation_graph.o $(SHERLOCK_DIR)/src/parsing_onnx.o
endif


OBJECTS = $(SOURCES:%.cpp=%.o)

all: $(OBJECTS) main

main : $(SOURCES) $(OBJECTS) 
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS_utils) $(OBJS_SHERLOCK) $@.o $(LIBS)
	mv ./main ./rino

clean:
	-rm *.o main

.cpp.o:
	$(CXX) $(CXXFLAGS)  -c $<

%.o: %.cpp $(SOURCES.h)
	$(CXX) $(CXXFLAGS) -c $(@:%.o=%.cpp) -o $@ 

