
#include <QPainter>
#include <stdio.h>
#include <math.h>
#include "window.h"
#include <QTextStream>
#include "approximate_tools.h"

#include "funcs.cpp"

#define DEFAULT_GRAPG_ACC 100;





Window::Window(QWidget *parent)
        : QWidget(parent) {
    a = -10;
    b = 10;
    n = 10;
    p = 0;
    mode = 0;
    func_id = 0;
    scale = 0;
    mode_name = "mode = 1";
    scale_name = QStringLiteral("scale = %1").arg(scale);
    p_name = QStringLiteral("p = %1").arg(p);
    der = new double[2];
    add_nodes = new double[4];
    x = new double[n];
    d = new double[n];
    f = new double[n];
    a_ans = new double[4*n];
}

QSize Window::minimumSizeHint() const {
    return QSize(100, 100);
}

QSize Window::sizeHint() const {
    return QSize(1000, 1000);
}

int Window::parse_command_line(int argc, char *argv[]) {
    if (argc != 5) {
        return -1;
    }

    if (sscanf(argv[1], "%lf", &a_args) != 1
        || sscanf(argv[2], "%lf", &b_args) != 1
        || fabs(b_args - a_args) < 1.e-6
        || sscanf(argv[3], "%d", &n) != 1
        || n <= 3 || sscanf(argv[4], "%d", &func_id) != 1
        || (func_id < 0 || (func_id > 6))) {
        return -2;
    }
    a = a_args;
    b = b_args;
    n_name = QStringLiteral("n = %1").arg(n);
    update_func();
    return 0;
}

//eval_appr(double x, double a_border, double b_border, int n, double *a, double *nodes)

void Window::draw_approx(double (*func)(double, double, double, int, double*, double*), QPainter &painter){
    double x1, x2, y1, y2;
    double max_y, min_y;
    double graph_acc = DEFAULT_GRAPG_ACC;
    double delta_y, delta_x = (b - a) / graph_acc / 10;
    QPen pen_black(Qt::black, 0, Qt::SolidLine);
    QPen pen_red(Qt::red, 0, Qt::SolidLine);
    QPen pen_green(Qt::green, 0, Qt::SolidLine);
    painter.setPen(pen_green);

    // calculate min and max for current function
    max_y = min_y = 0;
    for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
        y1 = (*func)(x1, a, b, n, a_ans, x);
        if (y1 < min_y)
            min_y = y1;
        if (y1 > max_y)
            max_y = y1;
    }

    delta_y = 0.01 * (max_y - min_y);
    min_y -= delta_y;
    max_y += delta_y;

    // save current Coordinate System
    painter.save();

    // make Coordinate Transformations
    painter.translate(0.5 * width(), 0.5 * height());
    painter.scale(width() / (b - a), -height() / (max_y - min_y));
    painter.translate(-0.5 * (a + b), -0.5 * (min_y + max_y));
    // draw approximated line for graph
    x1 = a;
    y1 = (*func)(x1, a, b, n, a_ans, x);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
        y2 = (*func)(x2, a, b, n, a_ans, x);
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        x1 = x2, y1 = y2;
    }
    x2 = b;
    y2 = (*func)(x2, a, b, n, a_ans, x);
    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

    // restore previously saved Coordinate System
    painter.restore();
}

void Window::draw_function(double (*func)(double, double, double, int, double*, double*, double, double, int), QPainter &painter){
    QTextStream out(stdout);
    double x1, x2, y1, y2;
    double max_y, min_y;
    double graph_acc = DEFAULT_GRAPG_ACC;
    double delta_y, delta_x = (b - a) / graph_acc / 10;
    QPen pen_black(Qt::black, 0, Qt::SolidLine);
    QPen pen_red(Qt::red, 0, Qt::SolidLine);
    QPen pen_blue(Qt::blue, 0, Qt::SolidLine);
    painter.setPen(pen_black);

    // calculate min and max for current function
    max_y = min_y = 0;
    for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
        y1 = (*func)(x1, a, b, n, a_ans, x, x[n/2], 0, p);
        if (y1 < min_y)
            min_y = y1;
        if (y1 > max_y)
            max_y = y1;
    }

    double max_f = max(fabs(max_y), fabs(min_y));
    for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
        y1 = (*func)(x1, a, b, n, a_ans, x, x[n/2], max_f, p);
        if (y1 < min_y)
            min_y = y1;
        if (y1 > max_y)
            max_y = y1;
    }
    delta_y = 0.01 * (max_y - min_y);
    min_y -= delta_y;
    max_y += delta_y;

    // save current Coordinate System
    painter.save();

    // make Coordinate Transformations
    painter.translate(0.5 * width(), 0.5 * height());
    painter.scale(width() / (b - a), -height() / (max_y - min_y));
    painter.translate(-0.5 * (a + b), -0.5 * (min_y + max_y));
    // draw approximated line for graph
    x1 = a;
    y1 = (*func)(x1, a, b, n, a_ans, x, x[n/2], max_f, p);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
        y2 = (*func)(x2, a, b, n, a_ans, x, x[n/2], max_f, p);
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        x1 = x2, y1 = y2;
    }
    x2 = b;
    y2 = (*func)(x2, a, b, n, a_ans, x, x[n/2], max_f, p);
    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

    // draw axis
    painter.setPen(pen_red);
    painter.drawLine(a, 0, b, 0);
    painter.drawLine(0, max_y, 0, min_y);
    out << QString("max(|F_max|, |F_min|) = %1").arg(max(fabs(max_y - delta_y), fabs(min_y + delta_y))) << endl;

    // restore previously saved Coordinate System
    painter.restore();
    painter.setPen(pen_blue);
    painter.drawText(0, 120, QString("max(|F_max|, |F_min|) = %1").arg(max(fabs(max_y - delta_y), fabs(min_y + delta_y))));
}

double Window::draw_residual(double (*func)(double, double, double, int, double*, double*, double, double, int), QPainter &painter, double max_y_f){
    QTextStream out(stdout);
    double x1, x2, y1, y2;
    double max_y, min_y;
    double graph_acc = 100;
    double delta_y, delta_x = (b - a) / graph_acc / 10;
    QPen pen_black(Qt::black, 0, Qt::SolidLine);
    QPen pen_red(Qt::red, 0, Qt::SolidLine);
    QPen pen_blue(Qt::blue, 0, Qt::SolidLine);
    painter.setPen(pen_black);

    // calculate min and max for current function
    max_y = min_y = 0;
    for (x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
        y1 = (*func)(x1, a, b, n, a_ans, x, x[n/2], max_y_f, p);
        if (y1 < min_y)
            min_y = y1;
        if (y1 > max_y)
            max_y = y1;
    }

    double max_f = max(fabs(max_y), fabs(min_y));

    delta_y = 0.01 * (max_y - min_y);
    min_y -= delta_y;
    max_y += delta_y;

    // save current Coordinate System
    painter.save();

    // make Coordinate Transformations
    painter.translate(0.5 * width(), 0.5 * height());
    painter.scale(width() / (b - a), -height() / (max_y - min_y));
    painter.translate(-0.5 * (a + b), -0.5 * (min_y + max_y));
    // draw approximated line for graph
    x1 = a;
    y1 = (*func)(x1, a, b, n, a_ans, x, x[n/2], max_y_f, p);
    for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) {
        y2 = (*func)(x2, a, b, n, a_ans, x, x[n/2], max_f, p);
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        x1 = x2, y1 = y2;
    }
    x2 = b;
    y2 = (*func)(x2, a, b, n, a_ans, x, x[n/2], max_y_f, p);
    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));

    // draw axis
    painter.setPen(pen_red);
    painter.drawLine(a, 0, b, 0);
    painter.drawLine(0, max_y, 0, min_y);
    out << QString("max(|F_max|, |F_min|) = %1").arg(max(fabs(max_y - delta_y), fabs(min_y + delta_y))) << endl;

    // restore previously saved Coordinate System
    painter.restore();
    painter.setPen(pen_blue);
    return max(fabs(max_y - delta_y), fabs(min_y + delta_y));
}

double Window::render_residual_1(QPainter &painter, double max_y){
    double m1 = 0;
    switch (func_id) {
        case 0:
            m1 = draw_residual(residual_1_0, painter, max_y);
            break;
        case 1:
            m1 = draw_residual(residual_1_1, painter, max_y);
            break;
        case 2:
            m1 = draw_residual(residual_1_2, painter, max_y);
            break;
        case 3:
            m1 = draw_residual(residual_1_3, painter, max_y);
            break;
        case 4:
            m1 = draw_residual(residual_1_4, painter, max_y);
            break;
        case 5:
            m1 = draw_residual(residual_1_5, painter, max_y);
            break;
        case 6:
            m1 = draw_residual(residual_1_6, painter, max_y);
            break;
    }
    return m1;
}

void Window::render_function(QPainter &painter){
    switch (func_id) {
        case 0:
            draw_function(f_0, painter);
            break;
        case 1:
            draw_function(f_1, painter);
            break;
        case 2:
            draw_function(f_2, painter);
            break;
        case 3:
            draw_function(f_3, painter);
            break;
        case 4:
            draw_function(f_4, painter);
            break;
        case 5:
            draw_function(f_5, painter);
            break;
        case 6:
            draw_function(f_6, painter);
            break;
    }
}



/// render graph
void Window::paintEvent(QPaintEvent * /* event */) {
    QTextStream out(stdout);
    QPainter painter(this);
    delete[] x;
    delete[] a_ans;
    delete[] f;
    delete[] d;
    x = new double[n];
    d = new double[n];
    f = new double[n];
    a_ans = new double[4*n];

    double max_y = 0, m1, m2;
    if(p != 0){
        double graph_acc = DEFAULT_GRAPG_ACC;
        double delta_x = (b - a) / graph_acc / 10;
        double y1 = 0;
        for (double x1 = a; x1 - b < 1.e-6; x1 += delta_x) {
            y1 = calc_f(func_id, x1, x, x[n/2], 0, p);
            if (fabs(y1) > max_y) {
                max_y = fabs(y1);
            }
        }
    }

    for (int i = 0; i < n; i++) {
        x[i] = a + (b - a) * i / (n - 1);
        f[i] = calc_f(func_id, x[i], x, x[n/2], max_y, p);
    }


    /* border info calc*/
    der[0] = calc_der(func_id, x[0]);
    der[1] = calc_der(func_id, x[n - 1]);
    add_nodes[0] = a - (b - a) / (n - 1);
    add_nodes[2] = a + (b - a) * n / (n - 1);
    add_nodes[1] = calc_f(func_id, add_nodes[0], x, 0, 0, 0);
    add_nodes[3] = calc_f(func_id, add_nodes[2], x, 0, 0, 0);

    switch(mode){
        case 0:
            render_function(painter);
            method_init_1(n, x, f, a_ans, d, der);
            draw_approx(eval_appr, painter);
            break;
        case 1:
            render_function(painter);
            method_init_2(n, x, f, a_ans, d, add_nodes);
            draw_approx(eval_appr, painter);
            break;
        case 2:
            render_function(painter);
            method_init_1(n, x, f, a_ans, d, der);
            draw_approx(eval_appr, painter);
            method_init_2(n, x, f, a_ans, d, add_nodes);
            draw_approx(eval_appr, painter);
            break;
        case 3:
            method_init_1(n, x, f, a_ans, d, der);
            m1 = render_residual_1(painter, max_y);
            method_init_2(n, x, f, a_ans, d, der);
            m2 = render_residual_1(painter, max_y);
            painter.drawText(0, 120, QString("max(|F_max|, |F_min|) = %1").arg(max(m1, m2)));
            break;
    }



    painter.setPen("blue");
    painter.drawText(0, 20, f_name);
    painter.drawText(0, 40, mode_name);
    painter.drawText(0, 60, scale_name);
    painter.drawText(0, 80, n_name);
    painter.drawText(0, 100, p_name);
}
