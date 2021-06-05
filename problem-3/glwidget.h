#ifndef _my_widget
#define _my_widget

#include <QOpenGLWidget>

#include <QKeyEvent>

#include <qnamespace.h>
#include "approx_tools.h"

class myGLWidget : public QOpenGLWidget {
    Q_OBJECT
public:
    int parse_command_line(int argc, char *argv[]);

protected:
    virtual void paintGL();

    virtual void initializeGL();

    virtual void resizeGL(int nWidth, int nHeight);

    virtual void keyPressEvent(QKeyEvent *e);

    void setProjectionMatrix();
    double calc_der_y(double x, double y);
    double calc_der_x(double x, double y);
    void set_func_name();
    void set_mode_name();
    double set_max_f();
    void update_names();
    void init_func();
    double calc_f(double x, double y);
    void setDefaultCamera();
    double ax, bx, ay, by, ax_args, bx_args, ay_args, by_args;
    int nx, ny, k, func_id;
    QString scale_name, p_name, n_name, mode_name, f_name, angle_name, max_f_name;
    float angle_h, angle_v;
    float camera_p;
    float aspect;
    int mode, scale, p;
    double max_f;
    double *f_vals;
    double *x_vals;
    double *y_vals;
    double *F;
    bool method_changed = true;
    int empty = 1;
};

#endif
