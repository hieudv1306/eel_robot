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

TEST_DIR := tmp/test_bins
TEST_CXXFLAGS := $(CXXFLAGS) -Iinclude
TEST_HEADERS := $(wildcard include/core/*.hpp include/physics/*.hpp)
TEST_BINS := \
	$(TEST_DIR)/test_gait \
	$(TEST_DIR)/test_material \
	$(TEST_DIR)/test_metrics \
	$(TEST_DIR)/test_soft_backbone \
	$(TEST_DIR)/test_geometry

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

$(TEST_DIR):
	mkdir -p $@

$(TEST_DIR)/test_gait: tests/test_gait.cpp src/physics/gait.cpp $(TEST_HEADERS) | $(TEST_DIR)
	$(CXX) $(TEST_CXXFLAGS) -o $@ $(filter %.cpp,$^)

$(TEST_DIR)/test_material: tests/test_material.cpp src/core/params.cpp src/physics/material.cpp src/physics/soft_rod.cpp $(TEST_HEADERS) | $(TEST_DIR)
	$(CXX) $(TEST_CXXFLAGS) -o $@ $(filter %.cpp,$^)

$(TEST_DIR)/test_metrics: tests/test_metrics.cpp src/physics/diagnostics.cpp $(TEST_HEADERS) | $(TEST_DIR)
	$(CXX) $(TEST_CXXFLAGS) -o $@ $(filter %.cpp,$^)

$(TEST_DIR)/test_soft_backbone: tests/test_soft_backbone.cpp src/core/params.cpp src/physics/gait.cpp src/physics/material.cpp src/physics/soft_rod.cpp src/physics/soft_backbone.cpp $(TEST_HEADERS) | $(TEST_DIR)
	$(CXX) $(TEST_CXXFLAGS) -o $@ $(filter %.cpp,$^)

$(TEST_DIR)/test_geometry: tests/test_geometry.cpp src/core/params.cpp src/physics/gait.cpp src/physics/geometry.cpp src/physics/markers.cpp src/physics/material.cpp src/physics/soft_rod.cpp src/physics/soft_backbone.cpp $(TEST_HEADERS) | $(TEST_DIR)
	$(CXX) $(TEST_CXXFLAGS) -o $@ $(filter %.cpp,$^)

.PHONY: test
test: $(TEST_BINS)
	@set -e; for test_bin in $(TEST_BINS); do \
		echo "Running $$test_bin"; \
		$$test_bin; \
	done

.PHONY: clean-tmp
clean-tmp:
	rm -f tmp/*.* tmp/vtkData/*.* tmp/vtkData/data/*.* tmp/imageData/*.* tmp/imageData/data/*.* tmp/gnuplotData/*.* tmp/gnuplotData/data/*.*

.PHONY: clean-core
clean-core:
	make -C $(OLB_ROOT) clean-core

.PHONY: clean
clean: clean-tmp
	rm -f $(OBJ_FILES) $(DEP_FILES) $(EXAMPLE)
	rm -rf $(TEST_DIR)

-include $(DEP_FILES)
