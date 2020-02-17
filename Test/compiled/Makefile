# location of the Boost Python include files and library
PYTHON_INCLUDE = `python3-config --includes`
PYTHON_LIB=`python3-config --libs`

BOOST_LIB=-lboost_python36 -lboost_numpy36
LINUX_INC=-I$(CONDA_PREFIX)/include/
LINUX_LIB=-L$(CONDA_PREFIX)/lib/

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix libs/,$(notdir $(CPP_FILES:.cpp=.o)))

MODULES := $(wildcard src/*.C)
MODULES_OBJ_FILES := $(addprefix ./,$(notdir $(MODULES:.C=.o)))
MODULES_SHARED_LIBS := $(addprefix ./,$(notdir $(MODULES:.C=.so)))



all: $(MODULES_SHARED_LIBS)

libs/%.o: src/%.cpp interface/%.h
	g++ -std=c++11 $(LINUX_INC) $(PYTHON_INCLUDE) -I./interface -O2 -fPIC -c -o $@ $< 

#pack helpers in lib
libpyhtt.so: $(OBJ_FILES) 
	g++ -std=c++11 -o $@ -shared -fPIC  $(LINUXADD) -fPIC  libs/*.o $(PYTHON_LIB) $(BOOST_LIB) $(LINUX_LIB)

%.o: src/%.C 
	g++ -std=c++11 -O2 -I./interface $(PYTHON_INCLUDE) $(LINUX_INC) -fPIC -c -o $(@) $<

%.so: %.o libpyhtt.so
	g++ -std=c++11 -o $(@) -shared -fPIC $(LINUXADD) $< $(PYTHON_LIB) $(BOOST_LIB) -L./ -lpyhtt $(LINUX_LIB)

clean: 
	rm -f $(OBJ_FILES) $(MODULES_SHARED_LIBS) $(MODULES_OBJ_FILES) libpyhtt.so
