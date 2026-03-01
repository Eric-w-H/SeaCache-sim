###################################################
# Constants
###################################################
TARGET := scache

SRC_DIRS := ./src
BUILD_DIR := ./build
OUTPUT_DIR := ./output

# find the source files, extract the filenames, 
# stick them in the build dir as .o
SRC := $(shell find $(SRC_DIRS) -name '*.cpp')
FILENAMES := $(basename $(notdir $(SRC)))
OBJS := $(FILENAMES:%=$(BUILD_DIR)/%.o)

# Flags for g++
CPPFLAGS := -O3 -Wall -Wextra -Werror

# Phony targets (do not represent a file)
.PHONY: clean

###################################################
# Targets
###################################################
$(TARGET): $(BUILD_DIR)/$(TARGET)
	@mkdir -p $(OUTPUT_DIR)
	cp $(BUILD_DIR)/$(TARGET) ./$(TARGET)

$(BUILD_DIR)/$(TARGET): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(BUILD_DIR)/$(TARGET)

# Construct a unique compilation step for each src->object
# to minimize re-building
define OBJ_COMP_TEMPLATE =
$(BUILD_DIR)/$(basename $(notdir $(1))).o: $(1)
	@mkdir -p $(BUILD_DIR)
	g++ -c $(CPPFLAGS) $(1) -o $$@
endef

# Instantiate & evaluate each of the object build steps
$(foreach f,$(SRC),$(eval $(call OBJ_COMP_TEMPLATE,$(f))))

clean:
	-rm -r $(BUILD_DIR) $(OUTPUT_DIR) $(TARGET)
