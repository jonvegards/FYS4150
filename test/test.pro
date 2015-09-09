TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

INCLUDEPATH += /usr/local/include
#LIBS += -L/usr/local/lib -llapack -lbase -larmadillo #-framework Accelerate

LIBS += -L/usr/local/lib -llapack -lblas -larmadillo

