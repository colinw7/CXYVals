TEMPLATE = app

TARGET = CQXYVals

DEPENDPATH += .

QT += widgets

QMAKE_CXXFLAGS += -std=c++11

#CONFIG += debug

# Input
SOURCES += \
CQXYVals.cpp \

HEADERS += \
CQXYVals.h \

DESTDIR     = ../bin
OBJECTS_DIR = ../obj
LIB_DIR     = ../lib

INCLUDEPATH += \
../include \
.

unix:LIBS += \
-L$$LIB_DIR \
-lCXYVals \
