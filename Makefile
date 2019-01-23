#===============================================================
# Makefile for building C Code in command line environement
# using the GCC Open Source Compiler
# Edited by Enzo IGLESIS
#===============================================================

######################################

RM						:= rm -rvf
MKDIR					= mkdir -p -- $@

DONE						= [\e[32mDONE\e[0;m]:
NOTE					= [\e[36mNOTE\e[0;m]:

BOLD					=\e[1m
UNDERLINE				=\e[4m
DEFAULT					=\e[0;m

######################################

ifndef TARGET
    $(info Target unspecified, default: TARGET=main.out)
    TARGET 			:= main.out
endif

INC_DIR				:= libraries

SOURCES 			:= $(patsubst src/%.c, %.c, $(wildcard src/*.c))
SOURCES             += $(patsubst $(INC_DIR)/src/%.c, %.c, $(wildcard $(INC_DIR)/src/*.c))

OBJ_DIR 			:= output


######################################

CC              		:= gcc
GDB			    		:= gdb
INCLUDES 				:= -I $(INC_DIR) -lm -lpthread
COMMONFLAGS 			:= -g -gstrict-dwarf -Wall -pthread
CFLAGS 					:= $(COMMONFLAGS)
LDFLAGS 				:= $(COMMONFLAGS)

######################################

VPATH := src:$(INC_DIR)/src/

OBJECTS := $(addsuffix .o,$(addprefix $(OBJ_DIR)/, $(basename $(SOURCES))))

######################################

all: $(OBJ_DIR)/$(TARGET)

$(OBJECTS): | $(OBJ_DIR)

$(OBJ_DIR):
	@echo "created directory '$(BOLD)$(UNDERLINE)$(OBJ_DIR)/$(DEFAULT)'"
	@$(MKDIR) $(OBJ_DIR)

$(OBJ_DIR)/%.o: %.c
	@echo ======================================================================
	@echo "$(NOTE) Generating $(BOLD)$(UNDERLINE)$@$(DEFAULT)"
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/$(TARGET): $(OBJECTS)
	@echo ======================================================================
	@echo "$(NOTE) Linking objects and generating output binary $(BOLD)$(UNDERLINE)$@$(DEFAULT)"
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(INCLUDES)
	@echo "$(DONE) output binary: $(BOLD)$(UNDERLINE)$@$(DEFAULT)"

debug: all
	$(GDB) $(OBJ_DIR)/$(TARGET)

clean:
	@$(RM) $(OBJ_DIR)

rebuild: clean all
