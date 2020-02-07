# Compiler options
CC = g++
OPTCFLAGS = -Ofast
CFLAGS = -Wall -std=c++11 $(OPTCFLAGS)
LDFLAGS_PSASS = -pthread -lstdc++ -lz -llzma -lbz2
LDFLAGS_KPOOL = -pthread -lstdc++ -lz

# Directory organisation
BASEDIR = .
BIN = $(BASEDIR)/bin
SRC = $(BASEDIR)/src
BUILD = $(BASEDIR)/build
INCLUDE = $(BASEDIR)/include
CPP = $(wildcard $(SRC)/*.cpp)

# Target
TARGETS = psass kpool

all: htslib init $(TARGETS)

init:
	mkdir -p $(BUILD) $(BUILD)
	mkdir -p $(BIN) $(BIN)

htslib:
	$(MAKE) -C include/htslib

clean-htslib:
	$(MAKE) -C include/htslib clean

psass: $(BUILD)/analyze.o  $(BUILD)/gff_file.o  $(BUILD)/output_handler.o  $(BUILD)/pair_data.o  $(BUILD)/pileup_converter.o  $(BUILD)/pileup.o  $(BUILD)/pool_data.o  $(BUILD)/psass.o
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/psass $^ $(INCLUDE)/htslib/libhts.a $(LDFLAGS_PSASS)

kpool: $(BUILD)/kpool.o $(BUILD)/kpool_merge.o
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/kpool $^ $(LDFLAGS_KPOOL)

$(BUILD)/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -c -o $@ $^

clean:
	rm -rf $(BUILD)/*.o
	rm -rf $(BIN)/*

clean-all: clean clean-htslib

rebuild: clean $(TARGETS)

rebuild-all: clean-all all
