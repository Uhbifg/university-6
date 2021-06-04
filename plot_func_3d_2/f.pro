TEMPLATE = app
TARGET = a
INCLUDEPATH += . include src

DEFINES += QT_DEPRECATED_WARNINGS

QT+= opengl core widgets gui 

#DESTDIR = bin/
OBJECTS_DIR = obj/
INCLUDE_DIR = include/
SOURCE_DIR = src/
MOC_DIR = moc/
CONFIG += warn_on qt
CONFIG -= debug_and_release

# QMAKE_LIBS += -fopenmp #-lgomp
# QMAKE_CXXFLAGS += -fopenmp
# QMAKE_LFLAGS += -fopenmp

# Input
HEADERS += \
	$${INCLUDE_DIR}/window.h \
	$${INCLUDE_DIR}/functions.h 
SOURCES += \
	$${SOURCE_DIR}/main.cpp \
	$${SOURCE_DIR}/window.cpp \
	$${SOURCE_DIR}/functions.cpp 

#                  ./a -1 1 -2 2 6 6 3
#                  ./a 0 1 2 5 6 6 3
#        ./a 0 1 0 1 2 2 5