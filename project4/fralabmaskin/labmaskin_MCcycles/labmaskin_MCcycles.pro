TEMPLATE = app
CONFIG += console c++11 -O3
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    lib.cpp \

HEADERS += \
    lib.h \

# MPI Settings
QMAKE_CXX = mpic++
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

# Settings for my own computer
QMAKE_CFLAGS += $$system(/usr/local/bin/mpic++ --showme:compile)
QMAKE_LFLAGS += $$system(/usr/local/bin/mpic++ --showme:link)
QMAKE_CXXFLAGS += $$system(/usr/local/bin/mpic++ --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(/usr/local/bin/mpic++ --showme:compile) -DMPICH_IGNORE_CXX_SEEK

# MPI Settings for terminal machines @ UiO
#QMAKE_CFLAGS += $$system(~/build/bin/mpic++ --showme:compile)
#QMAKE_LFLAGS += $$system(~/build/bin/mpic++ --showme:link)
#QMAKE_CXXFLAGS += $$system(~/build/bin/mpic++ --showme:compile) -DMPICH_IGNORE_CXX_SEEK
#QMAKE_CXXFLAGS_RELEASE += $$system(~/build/bin/mpic++ --showme:compile) -DMPICH_IGNORE_CXX_SEEK

