#TEMPLATE = app
#QT += qml quick
#CONFIG += c++11
#CONFIG -= app_bundle

TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include
#LIBS += -L/usr/local/lib -framework Accelerate
LIBS += -L/usr/local/lib -llapack -lblas -larmadillo

SOURCES += main.cpp
    #lu_decomposition.cpp

RESOURCES += qml.qrc

# Default rules for deployment.
include(deployment.pri)

#DISTFILES += \
#    plotting_results.py

HEADERS += \
    lu_decomposition.h
