# Compiler options
CC = g++
OPTCFLAGS = -Ofast
CFLAGS = -Wall -std=c++11 $(OPTCFLAGS)
LDFLAGS = -lstdc++ -lm

# Directory organisation
BASEDIR = .
BIN = $(BASEDIR)/bin
SRC = $(BASEDIR)/src
INCLUDE = $(BASEDIR)/include
BUILD = $(BASEDIR)/build
CPP = $(wildcard $(SRC)/*.cpp)

# Target
TARGET = psass

# Variables
OBJS = $(addprefix $(BUILD)/, $(notdir $(CPP:.cpp=.o)))

# Rules

all: init print-OBJS $(TARGET)

print-%  : ; @echo $* = $($*)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/$(TARGET) $^ $(LDFLAGS)

$(BUILD)/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -c -o $@ $^

clean:
	rm -rf $(BUILD)/*.o
	rm -rf $(BIN)/$(TARGET)

init:
	mkdir -p $(BUILD) $(BUILD)
	mkdir -p $(BIN) $(BIN)

rebuild: clean $(TARGET)
: clean $(TARGET)
