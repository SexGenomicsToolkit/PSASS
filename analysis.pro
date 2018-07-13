TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/main.cpp \
    src/arg_parser.cpp \
    src/analysis.cpp \
    src/gff_file.cpp

HEADERS += \
    src/arg_parser.h \
    src/parameters.h \
    src/analysis.h \
    src/utils.h \
    src/gff_file.h
