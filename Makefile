EXAMPLE = 11_lbm_eel_3dof
OLB_ROOT := ../openlb

-include $(OLB_ROOT)/config.mk
include $(OLB_ROOT)/rules.mk

INCLUDE_DIRS := include $(OLB_ROOT)/src

ifneq ($(USE_EMBEDDED_DEPENDENCIES), OFF)
INCLUDE_DIRS += \
	$(OLB_ROOT)/external/zlib \
	$(OLB_ROOT)/external/tinyxml2

LDFLAGS := -L$(OLB_ROOT)/external/lib $(LDFLAGS)

dependencies:
	$(MAKE) -C $(OLB_ROOT)/external

clean-dependencies:
	$(MAKE) -C $(OLB_ROOT)/external clean
else
.PHONY: dependencies clean-dependencies
endif

LDFLAGS += -L$(OLB_ROOT)/build/lib

# Enable OpenMP for manual LBM loops (independent of OpenLB's PARALLEL_MODE)
CXXFLAGS += -fopenmp
LDFLAGS  += -fopenmp

SRC_FILES := \
	src/main.cpp \
	$(wildcard src/core/*.cpp) \
	$(wildcard src/physics/*.cpp) \
	$(wildcard src/solver/*.cpp) \
	$(wildcard src/io/*.cpp)

OBJ_FILES := $(SRC_FILES:.cpp=.o)
DEP_FILES := $(SRC_FILES:.cpp=.d)
INCLUDE_FLAGS := $(addprefix -I,$(INCLUDE_DIRS))

all: dependencies core $(EXAMPLE)

core:
	$(MAKE) -C $(OLB_ROOT) core

%.d: %.cpp
	@$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $< -MM -MT $(@:.d=.o) >$@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c -o $@ $<

$(EXAMPLE): $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) -o $@ -lolbcore $(LDFLAGS)

.PHONY: onlysample
onlysample: $(EXAMPLE)

.PHONY: run
run: $(EXAMPLE)
	./$(EXAMPLE)

.PHONY: clean-tmp
clean-tmp:
	rm -f tmp/*.* tmp/vtkData/*.* tmp/vtkData/data/*.* tmp/imageData/*.* tmp/imageData/data/*.* tmp/gnuplotData/*.* tmp/gnuplotData/data/*.*

.PHONY: clean-core
clean-core:
	make -C $(OLB_ROOT) clean-core

.PHONY: clean
clean: clean-tmp
	rm -f $(OBJ_FILES) $(DEP_FILES) $(EXAMPLE)

-include $(DEP_FILES)
