#include "window.h"

void Window::update_p() {
    p_name = QStringLiteral("p = %1").arg(p);
    update();
}


void Window::update_n(){
    n_name = QStringLiteral("n = %1").arg(n);
    update();
}


void Window::update_mode(){
    mode_name = QStringLiteral("mode = %1").arg(mode);
    update();
}



/// update function for drawing (or calc)
void Window::update_func() {
    switch (func_id) {
        case 0:
            f_name = "k = 0, f (x) = 1";
            break;
        case 1:
            f_name = "k = 1, f (x) = x";
            break;
        case 2:
            f_name = "k = 2, f (x) = x * x";
            break;
        case 3:
            f_name = "k = 3, f (x) = x * x * x";
            break;
        case 4:
            f_name = "k = 4, f (x) = x * x * x * x";
            break;
        case 5:
            f_name = "k = 5, f (x) = exp(x)";
            break;
        case 6:
            f_name = "k = 6, f (x) = 1 / (25 x * x + 1)";
            break;
    }
    update();
}


void Window::change_func() {
    func_id = (func_id + 1) % 7;
    update_func();
}

void Window::change_mode(){
    mode = (mode + 1) % 4;
    update_mode();
}

void Window::upscale(){
    scale += 1;
    update_scale();
}


void Window::downscale(){
    scale -= 1;
    update_scale();
}


void Window::update_scale(){
    scale_name = QStringLiteral("scale = %1").arg(scale);
    a = a_args * pow(2, scale);
    b = b_args * pow(2, scale);
    update();
}


void Window::more_n(){
    n *= 2;
    update_n();
}


void Window::less_n(){
    if (n/2 >= 5){
        n /= 2;
        update_n();
    }else{
        QTextStream out(stdout);
        out << QString("too small n to decrease") << endl;
    }
}

void Window::less_p() {
    p -= 1;
    update_p();
}

void Window::more_p() {
    p += 1;
    update_p();
}