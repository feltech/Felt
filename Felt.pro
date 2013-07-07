#-------------------------------------------------
#
# Project created by QtCreator 2013-04-28T03:22:56
#
#-------------------------------------------------

QT      += core

TARGET = Felt
TEMPLATE = lib

DEFINES += FELT_LIBRARY

SOURCES += felt.cpp

HEADERS += felt.h\
        felt_global.h \
    felt/Surface.hpp \
    felt/Grid.hpp

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}

unix|win32: LIBS += -lgomp
