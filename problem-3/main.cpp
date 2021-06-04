#include <qapplication.h>
#include <qgl.h>
#include <QtWidgets/QMessageBox>

#include "glwidget.h"
#include "approx_tools.h"
int main(int argc, char** argv)
{
	QApplication app(argc, argv);
	myGLWidget myglw;


    if (myglw.parse_command_line(argc, argv)) {
        QMessageBox::warning(0, "Wrong input arguments!",
                             "Hint: min_x, max_x, min_y, max_y, nx, ny, func_id");
        return -1;
    }

	app.setActiveWindow(&myglw);
	myglw.show();
	return app.exec();
}
