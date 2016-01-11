QT = core gui widgets network printsupport svg

TARGET = qtgnuplotlib-example
TEMPLATE = app
LIBS += -lQtGnuplot

CONFIG += debug

SOURCES = qtgnuplotlib-example.cpp \
    visualization.cpp

FORMS += \
    visualization.ui

HEADERS += \
    visualization.h

