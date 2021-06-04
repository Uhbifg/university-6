#include "glwidget.h"
#include "math.h"
#include "funcs.h"

#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MAX(a, b) ((a) < (b) ? (b) : (a))

#define MAX_ACC 1000000
int myGLWidget::parse_command_line(int argc, char *argv[]) {
    if (argc != 8) {
        return -1;
    }

    if (sscanf(argv[1], "%lf", &ax) != 1
        || sscanf(argv[2], "%lf", &bx) != 1
        || sscanf(argv[3], "%lf", &ay) != 1
        || sscanf(argv[4], "%lf", &by) != 1
        || fabs(ax - bx) < 1.e-6
        || fabs(ay - by) < 1.e-6
        || sscanf(argv[5], "%d", &nx) != 1
        || sscanf(argv[6], "%d", &ny) != 1
        || nx < 5 || ny < 5
        || sscanf(argv[7], "%d", &func_id) != 1
        || (func_id < 0 || (func_id > 7))) {
        return -2;
    }

    mode = 0;
    scale = 0;
    p = 0;
    angle_h = 45;
    ax_args = ax;
    bx_args = bx;
    ay_args = ay;
    by_args = by;
    max_f = 0;
    init_func();
    update_names();
    return 0;
}

void myGLWidget::init_func(){
    max_f = 0;
    if(!f_vals){
        delete[] f_vals;
        delete[] x_vals;
        delete[] y_vals;
    }

    f_vals = new double[nx * ny];
    x_vals = new double[nx];
    y_vals = new double[ny];
    // init points
    double point = 0;
    for(int i = 0; i < nx; i++){
        point = (bx - ax) / (nx - 1) * i + ax;
        x_vals[i] = point;
    }
    for(int i = 0; i < ny; i++){
        point = (by - ay) / (ny - 1) * i + ay;
        y_vals[i] = point;
    }

    // init function values
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            f_vals[i + j * nx] = calc_f(x_vals[i], y_vals[j]);
        }
    }
    max_f = set_max_f();
    for(int i = 0; i < nx; i++){
        point = (bx - ax) / (nx - 1) * i + ax;
        x_vals[i] = point;
    }
    for(int i = 0; i < ny; i++){
        point = (by - ay) / (ny - 1) * i + ay;
        y_vals[i] = point;
    }

    // init function values
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            f_vals[i + j * nx] = calc_f(x_vals[i], y_vals[j]);
        }
    }
    update();
}


void myGLWidget::keyPressEvent(QKeyEvent* e){
    bool f_change = false;
    switch (e->key()) {
        case Qt::Key_C:
            setDefaultCamera();
            break;
        case Qt::Key_Up:
            angle_v = MIN(angle_v + 5.0, 80);
            break;
        case Qt::Key_Down:
            angle_v = MAX(angle_v - 5.0, -80);
            break;
        case Qt::Key_Left:
            angle_h -= 5.0;
            break;
        case Qt::Key_Right:
            angle_h += 5.0;
            break;
        case Qt::Key_Plus:
            camera_p = MAX(camera_p - 0.1, 7);
            break;
        case Qt::Key_Minus:
            camera_p += 0.1;
            break;
        case Qt::Key_0:
            f_change = true;
            func_id = (func_id + 1) % 8;
            break;
        case Qt::Key_1:
            mode = (mode + 1) % 3;
            break;
        case Qt::Key_2:
            scale += 1;
            ax = ax_args * pow(2, scale);
            bx = bx_args * pow(2, scale);
            ay = ay_args * pow(2, scale);
            by = by_args * pow(2, scale);
            break;
        case Qt::Key_3:
            scale -= 1;
            ax = ax_args * pow(2, scale);
            bx = bx_args * pow(2, scale);
            ay = ay_args * pow(2, scale);
            by = by_args * pow(2, scale);
            break;
        case Qt::Key_4:
            if(4 * nx * ny < MAX_ACC){
                f_change = true;
                nx *= 2;
                ny *= 2;
            }
            break;
        case Qt::Key_5:
            if(nx / 2 >= 5){
                nx /= 2;
                f_change = true;
            }
            if(ny / 2 >= 5){
                ny /= 2;
                f_change = true;
            }
            break;
        case Qt::Key_6:
            f_change = true;
            p += 1;
            break;
        case Qt::Key_7:
            f_change = true;
            p -= 1;
            break;
        case Qt::Key_8:
            angle_h = angle_h + 15.0;
            while(angle_h > 360){
                angle_h -= 360;
            }
            break;
        case Qt::Key_9:
            angle_h = angle_h - 15.0;
            while(angle_h < -360){
                angle_h += 360;
            }
            break;
    }
    if(f_change){
        init_func();
    }
    update_names();
    update();
}

void myGLWidget::update_names(){
    scale_name = QStringLiteral("scale = %1").arg(scale);
    p_name = QStringLiteral("p = %1").arg(p);
    n_name = QStringLiteral("nx = %1, ny = %2").arg(nx).arg(ny);
    angle_name = QStringLiteral("angle = %1 deg").arg(angle_h);
    max_f_name = QStringLiteral("func absolute maximum = %1").arg(set_max_f());
    set_func_name();
    set_mode_name();
}

void myGLWidget::set_mode_name() {
    switch (mode) {
        case 0:
            mode_name = QStringLiteral("mode 0 (function)");
            break;
        case 1:
            mode_name = QStringLiteral("mode 1 (first approximation)");
            break;
        case 2:
            mode_name = QStringLiteral("mode 2 (second approximation)");
            break;
    }
}

void myGLWidget::set_func_name(){
    switch (func_id) {
        case 0:
            f_name = QStringLiteral("f(x, y) = 1");
            break;
        case 1:
            f_name = QStringLiteral("f(x, y) = x");
            break;
        case 2:
            f_name = QStringLiteral("f(x, y) = y");
            break;
        case 3:
            f_name = QStringLiteral("f(x, y) = x + y");
            break;
        case 4:
            f_name = QStringLiteral("f(x, y) = sqrt(x ** 2 + y ** 2)");
            break;
        case 5:
            f_name = QStringLiteral("f(x, y) = x ** 2 + y ** 2");
            break;
        case 6:
            f_name = QStringLiteral("f(x, y) = exp(x ** 2 - y ** 2)");
            break;
        case 7:
            f_name = QStringLiteral("f(x, y) = 1 / (25 * (x ** 2 + y ** 2) + 1)");
            break;
    }
}


double myGLWidget::set_max_f(){
    double val1 = fabs(calc_f(x_vals[0], y_vals[0]));
    double val2 = fabs(calc_f(x_vals[0], y_vals[ny - 1]));
    double val3 = fabs(calc_f(x_vals[nx - 1], y_vals[0]));
    double val4 = fabs(calc_f(x_vals[nx - 1], y_vals[ny - 1]));
    double val5 = fabs(calc_f(x_vals[nx / 2], y_vals[nx / 2]));
    return MAX(MAX(MAX(MAX(val1, val2), val3), val4), val5);
}


double myGLWidget::calc_f(double x, double y) {
    double f;
    switch (func_id) {
        case 0:
            f = f0(x, y);
            break;
        case 1:
            f = f1(x, y);
            break;
        case 2:
            f = f2(x, y);
            break;
        case 3:
            f = f3(x, y);
            break;
        case 4:
            f = f4(x, y);
            break;
        case 5:
            f = f5(x, y);
            break;
        case 6:
            f = f6(x, y);
            break;
        case 7:
            f = f7(x, y);
            break;
    }
    if (p != 0 && fabs(x - x_vals[nx / 2]) < 1e-6 && fabs(y - y_vals[ny / 2]) < 1e-6){
        f += p * 0.1 * max_f;
    }
    return f;
}











