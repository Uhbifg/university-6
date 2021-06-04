#include "glwidget.h"
#include <qpainter.h>

#include <stdio.h>
#include <math.h>

#include "approx_tools.cpp"



#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MAX(a, b) ((a) < (b) ? (b) : (a))

#define X_DRAW_ACC 1
#define Y_DRAW_ACC 1

void myGLWidget::initializeGL()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);

	setDefaultCamera();
}

void myGLWidget::paintGL()
{
    QPainter painter(this);
    painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
	setProjectionMatrix();

	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glBegin(GL_QUADS);

	int	x_n = X_DRAW_ACC * nx, y_n = Y_DRAW_ACC * ny;
	double x, y, z;
	glColor3d(1.0,0.0,0.0);
    switch (mode) {
        case 0:
            for(int i = 0; i < x_n - 1; i++){
                for(int j = 0; j < y_n - 1; j++){
                    glColor3d(1.0 * (x_n - i) / x_n, 1.0 * j / y_n, 0.0);
                    x = (bx - ax) * i / (x_n - 1) + ax;
                    y = (by - ay) * j / (y_n - 1) + ay;
                    z = calc_f(x, y);
                    glVertex3d(x, y, z);
                    x = (bx - ax) * (i + 1) / (x_n - 1) + ax;
                    y = (by - ay) * j / (y_n - 1) + ay;
                    z = calc_f(x, y);
                    glVertex3d(x, y, z);
                    x = (bx - ax) * (i + 1) / (x_n - 1) + ax;
                    y = (by - ay) * (j + 1) / (y_n - 1) + ay;
                    z = calc_f(x, y);
                    glVertex3d(x, y, z);
                    x = (bx - ax) * i / (x_n - 1) + ax;
                    y = (by - ay) * (j + 1) / (y_n - 1) + ay;
                    z = calc_f(x, y);
                    glVertex3d(x, y, z);
                }
            }
            break;
        case 1:
            auto F = new double[nx * ny][16];
            printf("1");
            method_1_prod(f_vals, nx, ny, x_vals, y_vals, func_id, F);
            for(int i = 0; i < x_n - 1; i++){
                for(int j = 0; j < y_n - 1; j++){
                    glColor3d(1.0 * (x_n - i) / x_n, 1.0 * j / y_n, 0.0);
                    x = (bx - ax) * i / (x_n - 1) + ax;
                    y = (by - ay) * j / (y_n - 1) + ay;
                    z = method_compute(x, y, x_vals, y_vals, nx, ny, F);
                    glVertex3d(x, y, z);
                    x = (bx - ax) * (i + 1) / (x_n - 1) + ax;
                    y = (by - ay) * j / (y_n - 1) + ay;
                    z = method_compute(x, y, x_vals, y_vals, nx, ny, F);
                    glVertex3d(x, y, z);
                    x = (bx - ax) * (i + 1) / (x_n - 1) + ax;
                    y = (by - ay) * (j + 1) / (y_n - 1) + ay;
                    z = method_compute(x, y, x_vals, y_vals, nx, ny, F);
                    glVertex3d(x, y, z);
                    x = (bx - ax) * i / (x_n - 1) + ax;
                    y = (by - ay) * (j + 1) / (y_n - 1) + ay;
                    z = method_compute(x, y, x_vals, y_vals, nx, ny, F);
                    glVertex3d(x, y, z);
                }
            }
            break;
    }

	glEnd();

    painter.beginNativePainting();
    painter.setPen("blue");
    painter.drawText(0, 20, f_name);
    painter.drawText(0, 40, mode_name);
    painter.drawText(0, 60, scale_name);
    painter.drawText(0, 80, n_name);
    painter.drawText(0, 100, p_name);
    painter.drawText(0, 120, angle_name);
    painter.drawText(0, 140, max_f_name);
    painter.endNativePainting();
    painter.end();

	glDisable(GL_DEPTH_TEST);
}

void myGLWidget::resizeGL(int nWidth, int nHeight)
{
	glViewport(0, 0, nWidth, nHeight);
	aspect = 1.0 * nWidth / nHeight;
	update();
}

#define ANGLE_DIFF	(5)
#define POSITION_DIFF	(0.1)



void myGLWidget::setProjectionMatrix()
{
	GLfloat view[16] = {0}, projection[16] = {0}, tmp[16] = {0};

	GLfloat  vnear = 5, vtop = 2, vbottom = -vtop, vright = vtop * aspect, vleft = -vright;

	projection[0] = 2 * vnear / (vright - vleft);
	projection[8] = (vright + vleft) / (vright - vleft);
	projection[5] = 2 * vnear / (vtop - vbottom);
	projection[9] = (vtop + vbottom) / (vtop - vbottom);
	projection[10] = - 1;
	projection[14] = - 2 * vnear;
	projection[11] = -1;

	GLfloat	cam_x, cam_y, cam_z;
	cam_x = 0;
	cam_y = 0;
	cam_z = camera_p;

	view[0] = 1;
	view[6] = -1;
	view[9] = 1;
	view[15] = 1;

	view[12] = -cam_x;
	view[13] = -cam_y;
	view[14] = -cam_z;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glRotated(angle_h, 0, 0, 1);
	glGetFloatv(GL_PROJECTION_MATRIX, tmp);

	glLoadMatrixf(projection);
	glMultMatrixf(view);

	glRotated(angle_h, 0, 0, 1);
	glRotated(angle_v, tmp[0], tmp[4], tmp[8]);
}

void myGLWidget::setDefaultCamera()
{
	camera_p = 7;
	angle_h = 45;
	angle_v = 20;
	aspect = 1.0 * width() / height();
}
