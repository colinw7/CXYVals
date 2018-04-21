TEMPLATE = app

TARGET = CQXYVals

DEPENDPATH += .

MOC_DIR = .moc

QT += widgets

QMAKE_CXXFLAGS += -std=c++14

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
