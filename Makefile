CXX		  := g++
CXX_FLAGS := -g -W -std=c++11 -fPIC 
# CXX_FLAGS += -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer
LDPATH    :=
LD_FLAGS  := -lm -lpthread

BIN		:= bin
BUILD   := build
INCLUDE	:= include
LIB		:= lib
SRC		:= src

SRCEXT    := cpp
STATICEXT := a
SHAREDEXT := so

EXECUTE := $(BIN)/main
SOURCES := $(shell find $(SRC) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRC)/%, $(BUILD)/%, $(SOURCES:.$(SRCEXT)=.o))
DEPENDS := $(OBJECTS:%.o=%.d)

LIBNAME := 
LIBRARY := $(BIN)/lib$(LIBNAME).$(STATICEXT)
LIBSRCS :=
LIBSRCS += $(SRC)/$(LIBNAME).$(SRCEXT)
LIBOBJS := $(patsubst $(SRC)/%, $(BUILD)/%, $(LIBSRCS:.$(SRCEXT)=.o))

$(EXECUTE): $(OBJECTS)
	$(CXX) -o $(EXECUTE) $^ -L$(LDPATH) $(LD_FLAGS)

$(LIBRARY): $(OBJECTS)
	@echo $(LIBNAME)
	ar -rc $(LIBRARY) $(LIBOBJS)
	$(CXX) -shared $(CXX_FLAGS) -o $(LIBRARY:.$(STATICEXT)=.$(SHAREDEXT)) $(LIBOBJS)

$(BUILD)/%.o: $(SRC)/%.$(SRCEXT)
	@$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -MM -MT $@ -MF $(patsubst %.o, %.d, $@) $<
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -o $@ -c $< 

-include $(DEPENDS)

.PHONY:clean 
clean:
	-rm -r $(BIN)/* $(BUILD)/*

.PHONY:all 
all: $(EXECUTE) $(LIBRARY)

.PHONY:lib 
lib: $(LIBRARY)

run:$(EXECUTE)
	clear
	./$(EXECUTE)