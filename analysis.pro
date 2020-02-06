TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += include/

SOURCES += \
    pileup.cpp \
    src/main.cpp \
    src/gff_file.cpp \
    src/pileup.cpp \
    src/psass.cpp \
    src/pool_data.cpp \
    src/pair_data.cpp \
    src/output_handler.cpp \
    src/pileup_converter.cpp

HEADERS += \
    include/CLI11/CLI11.hpp \
    pileup.h \
    src/arg_parser.h \
    src/parameters.h \
    src/pileup.h \
    src/utils.h \
    src/gff_file.h \
    src/psass.h \
    src/pool_data.h \
    src/pair_data.h \
    src/input_data.h \
    src/output_handler.h \
    src/pileup_converter.h

DISTFILES += \
    include/CLI11/LICENSE
