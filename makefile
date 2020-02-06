# Compiler options
CC = g++
OPTCFLAGS = -Ofast
CFLAGS = -Wall -std=c++11 $(OPTCFLAGS)
LDFLAGS = -pthread -lstdc++ -lz -llzma -lbz2

# Directory organisation
BASEDIR = .
BIN = $(BASEDIR)/bin
SRC = $(BASEDIR)/src
BUILD = $(BASEDIR)/build
INCLUDE = $(BASEDIR)/include
CPP = $(wildcard $(SRC)/*.cpp)

# Target
TARGET = psass

# Variables
OBJS = $(addprefix $(BUILD)/, $(notdir $(CPP:.cpp=.o)))

# Rules

all: htslib init print-OBJS $(TARGET)

print-%  : ; @echo $* = $($*)

init:
	mkdir -p $(BUILD) $(BUILD)
	mkdir -p $(BIN) $(BIN)

htslib:
	$(MAKE) -C include/htslib

clean-htslib:
	$(MAKE) -C include/htslib clean

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/$(TARGET) $^ $(INCLUDE)/htslib/libhts.a $(LDFLAGS)

$(BUILD)/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -c -o $@ $^

clean:
	rm -rf $(BUILD)/*.o
	rm -rf $(BIN)/$(TARGET)

clean-all: clean clean-htslib

rebuild: clean $(TARGET)

rebuild-all: clean-all all
