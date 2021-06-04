#ifndef _my_widget
#define _my_widget

#include"functions.h"

#include <qgl.h>
#include <qnamespace.h>
#include <QVector2D>
#include <QVector3D>
#include <QtOpenGL>

#include<QGridLayout>

#include<QLineEdit>
#include<QTextEdit>
#include<QLabel>
#include<QPushButton>
#include<QSpinBox>
#include<QDoubleSpinBox>
#include<QComboBox>
#include<QCheckBox>

#include<QVector>
#include<QImage>
#include<Qt>
#include<QApplication>
#include<QMouseEvent> 
#include<QColor>

#include<QCoreApplication>

#include<QKeyEvent>
#include<QMainWindow>
#include<QWidget>
#include<QPainter>
#include<QTransform>


class GLWindow: public QGLWidget
{
  Q_OBJECT
	public:
		double AX, AY, BX, BY;
		int NX, NY, curFunction, functionsCount, curMethod, methodsCount;
		int np, span, p, curRot, maxRot;

		double FMin, FMax, FM;
		
		double Ax[16], Ay[16];
		double *X, *Y, *F, *DxF, *DyF, *DxDyF;
		double *precoeffs1, *coeffs1;//, *coeffs2;

		double (*funcarray[8])(double, double);
		double (*Dxfuncarray[8])(double, double);
		double (*Dyfuncarray[8])(double, double);
		double (*DxDyfuncarray[8])(double, double);


		GLWindow(double, double, double, double, int, int, int);
		
		~GLWindow()
		{
			if(X != NULL)
				free(X);
			if(Y != NULL)
				free(Y);
			if(F != NULL)
				free(F);
			if(DxF != NULL)
				free(DxF);
			if(DyF != NULL)
				free(DyF);
			if(DxDyF != NULL)
				free(DxDyF);
			if(precoeffs1 != NULL)
				free(precoeffs1);
			if(coeffs1 != NULL)
				free(coeffs1);
		}

		void findMinMax();

		void initXF();

		void setProjectionMatrix();
		void setDefaultCamera();

		float	angle_v;
		float	camera_p;

	    QMatrix4x4 projection;
	    QMatrix4x4 view;
    	QMatrix4x4 model;

	protected:
		virtual void paintGL();
		virtual void initializeGL();
		virtual void resizeGL(int nWidth, int nHeight);

		void drawFunction();
		void drawMethod1();
};


class Window:public QMainWindow
{
  Q_OBJECT
	public:

		QPoint pos = QPoint(50,50);
		QSize size = QSize(700,700);

		int toolbar_width = 250;

	// 	//____________________________________

		Window(double, double, double, double, int, int, int);

		~Window(){}

		int functionsCount, methodsCount, maxRot;

		QString methods[4], functions[8];
		
		GLWindow * glwindow;
		QWidget * toolbar;
		QGridLayout * grid;

		QLabel * NXLabel, * NYLabel, * KLabel, * SpanLabel, * PLabel, * RotLabel;
		QLabel  * AXLabel, * AYLabel, * BXLabel, * BYLabel;
		QLabel * MethodLabel, * FunctionLabel, * MinLabel, * MaxLabel;

		QSpinBox * NXBox, * NYBox, * KBox, * SpanBox, * PBox, * RotBox;
		QDoubleSpinBox * AXBox, * AYBox, * BXBox, * BYBox;

		void createSBox(QLabel*, QSpinBox*, int, int, int, int);
		void createDSBox(QLabel*, QDoubleSpinBox*, double, int);

		void updateMethodLabel()
		{
			MethodLabel->setText(methods[glwindow->curMethod]);
		}

		void updateFunctionLabel()
		{
			FunctionLabel->setText(functions[glwindow->curFunction]);
		}

		void updateMinMaxLabel()
		{
			MinLabel->setText(QString("FMin = %1").arg(glwindow->FMin));
			MaxLabel->setText(QString("FMax = %1").arg(glwindow->FMax));
		}

		void initGUI();

	public slots:
		void upd() 
		{ 
			this->update(); 
			this->updateMinMaxLabel(); 
			this->glwindow->updateGL(); 
		}


	protected:
		virtual void keyPressEvent(QKeyEvent *);
		virtual void resizeEvent(QResizeEvent *);
};

#endif
