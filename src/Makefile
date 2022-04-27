ROOT_DIR  := ../
SRC_DIR   := ./
BUILD_DIR := $(ROOT_DIR)build/
DEBUG_OUT_DIR := $(BUILD_DIR)debug/
TARGET 	  := $(DEBUG_OUT_DIR)connectivityAlg

RAND_INCLUDE := $(ROOT_DIR)/../random123/include/
JSON_INCLUDE := $(ROOT_DIR)/../json/include/
SRCS := $(shell find $(SRC_DIR) -name "*.cpp" | xargs -I {} basename {})
OBJS := $(SRCS:%.cpp=$(DEBUG_OUT_DIR)%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS  := $(shell find $(SRC_DIR) -type d)
INC_DIRS  += $(RAND_INCLUDE)
INC_DIRS  += $(JSON_INCLUDE)

INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CXX 	 := g++ 
CPPFLAGS := -g -Wall -Wextra -std=c++11 -pedantic
CXXFLAGS := $(INC_FLAGS) -MMD -MP 

LD 		 := g++
LD_FLAGS := -g -Wall -Wextra

MKDIR := mkdir -p
RM 	  := rm -rf

$(TARGET): $(OBJS)
	$(LD) $(LD_FLAGS) $(OBJS) -o $@

$(DEBUG_OUT_DIR)%.o: $(SRC_DIR)%.cpp 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(DEBUG_OUT_DIR)

-include $(DEPS)

