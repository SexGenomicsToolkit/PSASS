TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/main.cpp \
    src/output.cpp \
    src/utils.cpp \
    src/vfc_parsing.cpp \
    src/stats.cpp

HEADERS += \
    src/output.h \
    src/utils.h \
    src/variant_data.h \
    src/vcf_parsing.h \
    src/stats.h
