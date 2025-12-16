CXXFLAGS_BUILD=-DPW_DEFINE_SCALAR_DOUBLE
PRECISION=double
CXX = icpx
CXXFLAGS = -fPIC -g -O3 -std=c++14 -qopenmp -I$(PW_SDK_PATH)/include $(CXXFLAGS_BUILD)
CXXFLAGS += -I./third_party

SOURCES =  example.cpp
SOURCES += example_module.cpp
SOURCES += example_coupling_friction_sliding.cpp

OBJ_DIR = obj
OBJECTS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

BIN_DIR      = bin
TARGET_EXEC  = $(BIN_DIR)/pw.executable.example
LIBRARIES   += -lpw.api.$(PRECISION)
LIBRARIES   += -lpw.opengl.$(PRECISION)
LIBRARIES   += -lpw.viewer.$(PRECISION)
LIBRARIES   += -lpw.postprocess.$(PRECISION)
LIBRARIES   += -lpw.dem.$(PRECISION)
LIBRARIES   += -lpw.mps.$(PRECISION)
LIBRARIES   += -lpw.core.$(PRECISION)
LIBRARIES   += -lboost_filesystem
LIBRARIES   += -lboost_program_options
LIBRARIES   += -lboost_serialization

.PHONY : all prepare

all : prepare $(TARGET_EXEC)

clean :
	rm -rf bin obj

prepare :
	mkdir -p $(BIN_DIR) $(OBJ_DIR)

$(TARGET_EXEC) : $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) -L$(PW_SDK_PATH)/../lib  $(LIBRARIES)

$(OBJ_DIR)/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
