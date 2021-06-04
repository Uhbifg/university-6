
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>

#include "window.h"

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenuBar *toolbar = new QMenuBar(window);
    Window *graph_area = new Window(window);
    QAction *action;

    if (graph_area->parse_command_line(argc, argv)) {
        QMessageBox::warning(0, "Wrong input arguments!",
                             "Wrong input arguments!");
        return -1;
    }

    action = toolbar->addAction("&Change function", graph_area, SLOT(change_func()));
    action->setShortcut(QString("0"));

    action = toolbar->addAction("E&xit", window, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));

    action = toolbar->addAction("&Change mode", graph_area, SLOT(change_mode()));
    action->setShortcut(QString("1"));

    action = toolbar->addAction("&Upsacle", graph_area, SLOT(upscale()));
    action->setShortcut(QString("2"));

    action = toolbar->addAction("&Downscale", graph_area, SLOT(downscale()));
    action->setShortcut(QString("3"));

    action = toolbar->addAction("&n * 2", graph_area, SLOT(more_n()));
    action->setShortcut(QString("4"));

    action = toolbar->addAction("&n / 2", graph_area, SLOT(less_n()));
    action->setShortcut(QString("5"));

    action = toolbar->addAction("&p++", graph_area, SLOT(more_p()));
    action->setShortcut(QString("6"));

    action = toolbar->addAction("&p--", graph_area, SLOT(less_p()));
    action->setShortcut(QString("7"));


    toolbar->setMaximumHeight(30);
    window->setMenuBar(toolbar);
    window->setCentralWidget(graph_area);
    window->setWindowTitle("Graph");

    window->show();
    app.exec();
    delete window;
    return 0;
}
