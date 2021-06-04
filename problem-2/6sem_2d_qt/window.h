
#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

class Window : public QWidget
{
	Q_OBJECT

private:
	int func_id;
	QString f_name;
    QString mode_name;
    QString scale_name;
    QString n_name;
    QString p_name;
	double a;
	double b;
	double a_args;
	double b_args;
    int n;
	double *f;
	double *x;
    double *a_ans; // 4*n
    double *d; // n
    double *der; // 2
    double *add_nodes; // 4
	int mode;
    int scale;
    int p;
public:
	Window(QWidget *parent);

	QSize minimumSizeHint() const;
	QSize sizeHint() const;

	int parse_command_line(int argc, char *argv[]);
	void update_func();
	void update_mode();
    void update_scale();
    void update_n();
    void update_p();
    void draw_function(double (*func)(double, double, double, int, double*, double*, double, double, int), QPainter &painter);
    void render_function(QPainter &painter);
    void draw_approx(double (*func)(double, double, double, int, double*, double*), QPainter &painter);
    double render_residual_1(QPainter &painter, double max_y);
    double draw_residual(double (*func)(double, double, double, int, double*, double*, double, double, int), QPainter &painter, double max_y_f);
public slots:
    void more_n();
    void less_n();
    void more_p();
    void less_p();
    void upscale();
    void downscale();
	void change_func();
    void change_mode();
protected:
	void paintEvent(QPaintEvent *event);
};

#endif
