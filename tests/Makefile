ROOT_DIR  	    := ../
TESTS_DIR 	    := ./
SRC_DIR         := $(ROOT_DIR)src/
BUILD_DIR 	    := $(ROOT_DIR)build/
TESTS_OUT_DIR   := $(BUILD_DIR)tests/
DEBUG_OUT_DIR   := $(BUILD_DIR)debug/
TARGET    	    := $(TESTS_OUT_DIR)connectionTest

RAND_INCLUDE := $(ROOT_DIR)../random123/include/
JSON_INCLUDE := $(ROOT_DIR)../json/include/

TESTS := $(shell find $(TESTS_DIR) -name "*.cpp" | xargs -I {} basename {})
SRCS  := $(shell find $(SRC_DIR) -name "*.cpp" | xargs -I {} basename {})
EXCLUDES := main.cpp
SRCS :=$(filter-out $(EXCLUDES), $(SRCS))
TEST_OBJS  := $(TESTS:%.cpp=$(TESTS_OUT_DIR)%.o)
SRC_OBJS   += $(SRCS:%.cpp=$(DEBUG_OUT_DIR)%.o)
DEPS  := $(OBJS:.o=.d)

INC_DIRS  := $(shell find $(TESTS_DIR) -type d)
INC_DIRS  += $(SRC_DIR)
INC_DIRS  += $(RAND_INCLUDE)
INC_DIRS  += $(JSON_INCLUDE)

INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CXX 	 := g++
CPPFLAGS := -g -Wall -Wextra -std=c++11 -pedantic
CXXFLAGS := $(INC_FLAGS) -MMD -MP 

LD 		 := g++
LD_FLAGS := -g -Wall -Wextra

MKDIR := mkdir -p
RM    := rm -rf

$(TARGET): $(TEST_OBJS) $(SRC_OBJS)
	$(LD) $(LD_FLAGS) $(TEST_OBJS) $(SRC_OBJS) -o $@

$(TESTS_OUT_DIR)%.o: $(TESTS_DIR)%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(DEBUG_OUT_DIR)%.o: $(SRC_DIR)%.cpp 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(TESTS_OUT_DIR)

-include $(DEPS)

