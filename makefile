# Compiler options
CC = g++
OPTCFLAGS = -Ofast
CFLAGS = -Wall -std=c++11 $(OPTCFLAGS)
LDFLAGS = -pthread -static-libstdc++ -lz

# Directory organisation
BASEDIR = .
BIN = $(BASEDIR)/bin
SRC = $(BASEDIR)/src
BUILD = $(BASEDIR)/build
CPP = $(wildcard $(SRC)/*.cpp)

# Target
TARGET = poolsex

# Variables
OBJS = $(addprefix $(BUILD)/, $(notdir $(CPP:.cpp=.o)))

# Rules

all: init print-OBJS $(TARGET)

print-%  : ; @echo $* = $($*)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$(TARGET) $^ $(LDFLAGS)

$(BUILD)/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

clean:
	rm -rf $(BUILD)/*.o
	rm -rf $(BIN)/$(TARGET)

init:
	mkdir -p $(BUILD) $(BUILD)
	mkdir -p $(BIN) $(BIN)

rebuild: clean $(TARGET)
