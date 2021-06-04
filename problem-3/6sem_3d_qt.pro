QT += opengl
win32:LIBS += -lOpenGL32

HEADERS += glwidget.h \
           funcs.h
SOURCES += main.cpp \
	   glwidget.cpp \
	   glwidget_tools.cpp \
	   funcs.cpp \
	   approx_tools.cpp
