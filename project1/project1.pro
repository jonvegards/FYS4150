TEMPLATE = app

QT += qml quick
CONFIG += c++11

CONFIG -= app_bundle

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -framework Accelerate


SOURCES += main.cpp \
    #lu_decomposition.cpp
    lu_decomposition.cpp

RESOURCES += qml.qrc

# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =

# Default rules for deployment.
include(deployment.pri)

DISTFILES += \
    plotting_results.py

HEADERS += \
    lu_decomposition.h
