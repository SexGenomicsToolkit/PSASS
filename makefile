ifndef HTSLIB_CORES
	HTSLIB_CORES = 1
endif

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
TARGETS = $(BIN)/psass $(BIN)/kpool

.PHONY: all
all: include/htslib $(BIN) $(BUILD) $(TARGETS)

$(BUILD):
	mkdir -p $(BUILD)

$(BIN):
	mkdir -p $(BIN)

include/htslib_configured:
	cd include/htslib && ./configure
	touch $@

include/htslib: include/htslib_configured
	$(MAKE) -C include/htslib -j $(HTSLIB_CORES)

clean-htslib:
	rm include/htslib_configured
	$(MAKE) -C include/htslib clean

$(BIN)/psass: $(BUILD)/analyze.o  $(BUILD)/gff_file.o  $(BUILD)/output_handler.o  $(BUILD)/pair_data.o  $(BUILD)/pileup_converter.o  $(BUILD)/pileup.o  $(BUILD)/pool_data.o  $(BUILD)/psass.o
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/psass $^ $(INCLUDE)/htslib/libhts.a $(LDFLAGS_PSASS)

$(BIN)/kpool: $(BUILD)/kpool.o $(BUILD)/kpool_merge.o $(BUILD)/kpool_filter.o
	$(CC) $(CFLAGS) -I $(INCLUDE) -o $(BIN)/psass $^ $(LDFLAGS_KPOOL)

$(BUILD)/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -c -o $@ $^

.PHONY: clean
clean:
	rm -rf $(BUILD)/*.o
	rm -rf $(BIN)/*

clean-all: clean clean-htslib

rebuild: clean $(TARGETS)

rebuild-all: clean-all all
