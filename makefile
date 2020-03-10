# Compiler options
CC = gcc
CXX = g++
OPTCXXFLAGS = -Ofast
CXXFLAGS = -Wall -std=c++11 $(OPTCFLAGS)
LDFLAGS_PSASS = -pthread -lstdc++ -lz -llzma -lbz2

# Directory organisation
BASEDIR = .
BIN = $(BASEDIR)/bin
BUILD = $(BASEDIR)/build
INCLUDE = $(BASEDIR)/include
SRC = $(BASEDIR)/src

# Get number of parallel jobs
MAKE_PID := $(shell echo $$PPID)
JOBS := $(shell ps T | sed -n 's/.*$(MAKE_PID).*$(MAKE).* \(-j\|--jobs\) *\([0-9][0-9]*\).*/\2/p')
ifeq ($(JOBS),)
	JOBS := 1
endif

# Dummy flag to indicate that htslib was configured
HTSLIB_CONF_FLAG = $(INCLUDE)/.htslib_configured

# Target
TARGETS = $(BIN)/psass

# Declare phony targets (i.e. targets which are not files)
.PHONY: all clean clean-all rebuild rebuild-all

# Main rule
all: $(BIN) $(BUILD) $(TARGETS)

# Build directory
$(BUILD):
	if [ ! -e $@ ]; then \
		mkdir $@;        \
	fi;

# Bin directory
$(BIN):
	if [ ! -e $@ ]; then \
		mkdir $@;        \
	fi;

# Special rule to configure htslib with autoconf only on first build and full rebuild (using a dummy flag file)
$(HTSLIB_CONF_FLAG):
	cd $(INCLUDE)/htslib && ./configure CC=$(CC) --disable-libcurl
	touch $@

# Build htslib
$(INCLUDE)/htslib/libhts.a: $(HTSLIB_CONF_FLAG)
	$(MAKE) CC=$(CC) -C include/htslib -j $(JOBS)

# Clean htslib (run make clean and remove configure dummy flag)
clean-htslib:
	rm -f $(HTSLIB_CONF_FLAG)
	$(MAKE) CC=$(CC) -C include/htslib clean

# Linking for psass
$(BIN)/psass: $(BUILD)/analyze.o  $(BUILD)/gff_file.o  $(BUILD)/output_handler.o  $(BUILD)/pair_data.o  $(BUILD)/pileup_converter.o  $(BUILD)/pileup.o  $(BUILD)/pool_data.o  $(BUILD)/psass.o $(INCLUDE)/htslib/libhts.a
	$(CXX) $(CXXFLAGS) -I $(INCLUDE) -o $(BIN)/psass $^ $(LDFLAGS_PSASS)

# Build a single object file. Added htslib as dependency so that it is build before object files
$(BUILD)/%.o: $(SRC)/%.cpp $(INCLUDE)/htslib/libhts.a
	$(CXX) $(CXXFLAGS) -I $(INCLUDE) -c -o $@ $<

# Clean PSASS files
clean:
	rm -rf $(BUILD)/*.o
	rm -rf $(BIN)/*

# Clean all files
clean-all: clean clean-htslib

# Rebuild PSASS only
rebuild:
	$(MAKE) clean
	$(MAKE) -j $(JOBS)

# Rebuild PSASS and dependencies
rebuild-all:
	$(MAKE) clean-all
	$(MAKE) -j $(JOBS)
