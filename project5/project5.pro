TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    arma_solve.cpp \
    num_solve.cpp \
    lib.cpp

HEADERS += \
    arma_solve.h \
    num_solve.h \
    lib.h
QMAKE_MAC_SDK = macosx10.11 # Something wrong with SDK... (idk lol)
INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -larmadillo
