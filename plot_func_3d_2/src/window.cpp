#include "window.h"

#include <stdio.h>
#include <math.h>


#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MAX(a, b) ((a) < (b) ? (b) : (a))
//==============================================
void GLWindow::findMinMax()
{
	double a11, a12, a21, a22;

	a11 = funcarray[curFunction](AX,AY);
	a12 = funcarray[curFunction](AX,BY);
	a21 = funcarray[curFunction](BX,AY);
	a22 = funcarray[curFunction](BX,BY);

	if (3 >= curFunction && curFunction >= 1)
	{
		FMin = a11;
		FMax = a22;
	}

	switch(curFunction)
	{
		case 0:
			FMin = 0;
			FMax = 2;
			break;
		case 4:
			if(AX <= 0 && BX >= 0 && AY <= 0 && BY >= 0)
				FMin = 0;
			else if (AX <= 0 && BX >= 0)
				FMin = min(f4(0, AY), f4(0, BY));
			else if (AY <= 0 && BY >= 0)
				FMin = min(f4(AX, 0), f4(BX, 0));
			else
				FMin = min(min(a11, a12), min(a21, a22));
			FMax = max(max(a11, a12), max(a21, a22));
			break;
		case 5:
			if(AX <= 0 && BX >= 0 && AY <= 0 && BY >= 0)
				FMin = 0;
			else if (AX <= 0 && BX >= 0)
				FMin = min(f5(0, AY), f5(0, BY));
			else if (AY <= 0 && BY >= 0)
				FMin = min(f5(AX, 0), f5(BX, 0));
			else
				FMin = min(min(a11, a12), min(a21, a22));
			FMax = max(max(a11, a12), max(a21, a22));
			break;
		case 6:
			if(AX <= 0 && BX >= 0 && AY <= 0 && BY >= 0){
				FMin = min(f6(0, AY), f6(0, BY));
				FMax = max(f6(AX, 0), f6(BX, 0));
			}
			else if (AX <= 0 && BX >= 0){
				FMin = min(f6(0, AY), f6(0, BY));
				FMax = max(max(a11, a12), max(a21, a22));
			}
			else if (AY <= 0 && BY >= 0){
				FMin = min(min(a11, a12), min(a21, a22));
				FMax = max(f6(AX, 0), f6(BX, 0));
			}
			else{
				FMin = min(min(a11, a12), min(a21, a22));
				FMax = max(max(a11, a12), max(a21, a22));
			}
			break;
		case 7:
			if(AX <= 0 && BX >= 0 && AY <= 0 && BY >= 0)
				FMax = 1;
			else if (AX <= 0 && BX >= 0)
				FMax = max(f7(0, AY), f7(0, BY));
			else if (AY <= 0 && BY >= 0)
				FMax = max(f7(AX, 0), f7(BX, 0));
			else
				FMax = max(max(a11, a12), max(a21, a22));
			FMin = min(min(a11, a12), min(a21, a22));
			break;
	}

	FM = max(fabs(FMin), fabs(FMax));
}

Window::Window(const double ax, const double ay, const double bx, const double by,
					int nx, int ny, int k)
{
	functionsCount = 8;
	// curMethod = 0;
	methodsCount = 4;

	// p = 0;
	// span = 0;
	// curRot = 0;
	maxRot = 72;

	functions[0] = QString("f(x) = 1");
	functions[1] = QString("f(x) = x");
	functions[2] = QString("f(x) = y");
	functions[3] = QString("f(x) = x + y");
	functions[4] = QString("f(x) = sqrt(x^2 + y^2)");
	functions[5] = QString("f(x) = x^2 + y^2");
	functions[6] = QString("f(x) = exp(x^2 - y^2)");
	functions[7] = QString("f(x) = 1/(25(x^2+y^2)+1)");

	methods[0] = QString("method 1");
	methods[1] = QString("method 2");
	methods[2] = QString("methods 1 and 2");
	methods[3] = QString("function");

	glwindow = new GLWindow(ax, ay, bx, by, nx, ny, k);
	glwindow->resize(size);
	glwindow->move(QPoint(toolbar_width,0));
	glwindow->show();	

	initGUI();
}

GLWindow::GLWindow(const double ax, const double ay, const double bx, const double by,
					int nx, int ny, int k):
	AX(ax), AY(ay), BX(bx), BY(by),
	NX(nx), NY(ny),
	curFunction(k)
{
	functionsCount = 8;
	curMethod = 0;
	methodsCount = 4;

	p = 0;
	span = 0;
	curRot = 9;
	maxRot = 72;

	funcarray[0] = f0;
	funcarray[1] = f1;
	funcarray[2] = f2;
	funcarray[3] = f3;
	funcarray[4] = f4;
	funcarray[5] = f5;
	funcarray[6] = f6;
	funcarray[7] = f7;

	Dxfuncarray[0] = Dxf0;
	Dxfuncarray[1] = Dxf1;
	Dxfuncarray[2] = Dxf2;
	Dxfuncarray[3] = Dxf3;
	Dxfuncarray[4] = Dxf4;
	Dxfuncarray[5] = Dxf5;
	Dxfuncarray[6] = Dxf6;
	Dxfuncarray[7] = Dxf7;

	Dyfuncarray[0] = Dyf0;
	Dyfuncarray[1] = Dyf1;
	Dyfuncarray[2] = Dyf2;
	Dyfuncarray[3] = Dyf3;
	Dyfuncarray[4] = Dyf4;
	Dyfuncarray[5] = Dyf5;
	Dyfuncarray[6] = Dyf6;
	Dyfuncarray[7] = Dyf7;

	DxDyfuncarray[0] = DxDyf0;
	DxDyfuncarray[1] = DxDyf1;
	DxDyfuncarray[2] = DxDyf2;
	DxDyfuncarray[3] = DxDyf3;
	DxDyfuncarray[4] = DxDyf4;
	DxDyfuncarray[5] = DxDyf5;
	DxDyfuncarray[6] = DxDyf6;
	DxDyfuncarray[7] = DxDyf7;

	X = NULL; Y = NULL;
	F = NULL; DxF = NULL; DyF = NULL; DxDyF = NULL;
	precoeffs1 = NULL; coeffs1 = NULL;

	initXF();

	glClearColor(1.0, 1.0, 1.0, 1.0);

	setDefaultCamera();
}
//==============================================
void GLWindow::initXF()
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

	X = (double*)malloc(NX * sizeof(double));
	Y = (double*)malloc(NY * sizeof(double));
	F = (double*)malloc(NX * NY * sizeof(double));
	DxF = (double*)malloc(NX * NY * sizeof(double));
	DyF = (double*)malloc(NX * NY * sizeof(double));
	DxDyF = (double*)malloc(NX * NY * sizeof(double));

	precoeffs1 = (double*)malloc(NX * NY * 16 * sizeof(double));
	coeffs1 = (double*)malloc(NX * NY * 16 * sizeof(double));

	double dx = (BX - AX)/(NX - 1);
	double dy = (BY - AY)/(NY - 1);

	for (int i = 0; i < NX; ++i)
		X[i] = i * dx + AX;

	for (int i = 0; i < NY; ++i)
		Y[i] = i * dy + AY;

	for (int i = 0; i < NX; ++i)
		for (int j = 0; j < NY; ++j)
			{
				F[j*NX + i] = funcarray[curFunction](X[i], Y[j]);
				DxF[j*NX + i] = Dxfuncarray[curFunction](X[i], Y[j]);
				DyF[j*NX + i] = Dyfuncarray[curFunction](X[i], Y[j]);
				DxDyF[j*NX + i] = DxDyfuncarray[curFunction](X[i], Y[j]);
			}


	double hx = 1.0/dx, hx2 = 1.0/(dx*dx);
	double hy = 1.0/dy, hy2 = 1.0/(dy*dy);

	Ax[0] = 1; Ax[1] = 0; Ax[2] = 0; Ax[3] = 0;
	Ax[4] = 0; Ax[5] = 1; Ax[6] = 0; Ax[7] = 0;
	Ax[8] = -3*hx2; Ax[9] = -2*hx; Ax[10] = 3*hx2; Ax[11] = -hx;
	Ax[12] = 2*hx2; Ax[13] = hx; Ax[14] = -2*hx2; Ax[15] = hx;

	Ay[0] = 1; Ay[1] = 0; Ay[2] = 0; Ay[3] = 0;
	Ay[4] = 0; Ay[5] = 1; Ay[6] = 0; Ay[7] = 0;
	Ay[8] = -3*hy2; Ay[9] = -2*hy; Ay[10] = 3*hy2; Ay[11] = -hy;
	Ay[12] = 2*hy2; Ay[13] = hy; Ay[14] = -2*hy2; Ay[15] = hy;


	findMinMax();
	printf("%lf \n", FM);

	F[(NY/2) * NX + NX/2] += p * FM / 10;


	for (int i = 0; i < NX; ++i){
		for (int j = 0; j < NY; ++j)
		{
			precoeffs1[(j*NX + i)*16] = F[j*NX + i];
			precoeffs1[(j*NX + i)*16 + 2] = F[(j+1)*NX + i];
			precoeffs1[(j*NX + i)*16 + 8] = F[j*NX + i+1];
			precoeffs1[(j*NX + i)*16 + 10] = F[(j+1)*NX + i+1];

			precoeffs1[(j*NX + i)*16 + 4] = DxF[j*NX + i];
			precoeffs1[(j*NX + i)*16 + 6] = DxF[(j+1)*NX + i];
			precoeffs1[(j*NX + i)*16 + 12] = DxF[j*NX + i+1];
			precoeffs1[(j*NX + i)*16 + 14] = DxF[(j+1)*NX + i+1];

			precoeffs1[(j*NX + i)*16 + 1] = DyF[j*NX + i];
			precoeffs1[(j*NX + i)*16 + 3] = DyF[(j+1)*NX + i];
			precoeffs1[(j*NX + i)*16 + 9] = DyF[j*NX + i+1];
			precoeffs1[(j*NX + i)*16 + 11] = DyF[(j+1)*NX + i+1];

			precoeffs1[(j*NX + i)*16 + 5] = DxDyF[j*NX + i];
			precoeffs1[(j*NX + i)*16 + 7] = DxDyF[(j+1)*NX + i];
			precoeffs1[(j*NX + i)*16 + 13] = DxDyF[j*NX + i+1];
			precoeffs1[(j*NX + i)*16 + 15] = DxDyF[(j+1)*NX + i+1];
		}
	}

	calculateCoeffs1(NX, NY, Ax, Ay, precoeffs1, coeffs1);
	// calculateCoeffs1(NX, NY, precoeffs1, coeffs1, Ax, Ay);

	// printf("%lf\n", coeffs1[0]);

	// double dx2 = 1.0/dx, dx3 = 1.0/(dx*dx);
	// double dy2 = 1.0/dy, dy3 = 1.0/(dy*dy);

	// // double dx2 = 1.0/(dx*dx), dx3 = 1.0/(dx*dx*dx);
	// // double dy2 = 1.0/(dy*dy), dy3 = 1.0/(dy*dy*dy);

	// Ax[0]=1; Ax[1]=0; Ax[2]=0; Ax[3]=0;
	// Ax[4]=0; Ax[5]=1; Ax[6]=0; Ax[7]=0;
	// Ax[8]=-3*dx3; Ax[9]=-2*dx2; Ax[10]=3*dx3; Ax[11]=-dx2;
	// Ax[12]=2*dx3; Ax[13]=dx2; Ax[14]=-2*dx3; Ax[15]=dx2;

	// Ay[0]=1; Ay[1]=0; Ay[2]=0; Ay[3]=0;
	// Ay[4]=0; Ay[5]=1; Ay[6]=0; Ay[7]=0;
	// Ay[8]=-3*dy3; Ay[9]=-2*dy2; Ay[10]=3*dy3; Ay[11]=-dy2;
	// Ay[12]=2*dy3; Ay[13]=dy2; Ay[14]=-2*dy3; Ay[15]=dy2;

// 	// calculateCoeffs2(N, X, F, coeffs2, DF2);
}


//=================================================
//=================================================

void Window::initGUI()
{	
	resize(toolbar_width, size.height());
	move(pos);

	toolbar = new QWidget(this);
	toolbar->resize(toolbar_width, size.height());
	// toolbar->move(QPoint(size.width() - toolbar_width,0));
	toolbar->move(QPoint(0,0));
	toolbar->setStyleSheet(".QWidget {background-color: #123;}");	

	grid = new QGridLayout(toolbar);
	// grid = new QGridLayout(this);
	grid->setSpacing(10);
	toolbar->setLayout(grid);

	//_________________________________

	int h = 0;

	NXLabel = new QLabel("nx = ");
	NXBox = new QSpinBox(this);
	createSBox(NXLabel, NXBox, glwindow->NX, 2, 500, h); // set min=5
	QObject::connect(NXBox, QOverload<int>::of(&QSpinBox::valueChanged),
	[=](int i){ 
		glwindow->NX = i;
		glwindow->initXF();
	});
	++h;

	NYLabel = new QLabel("ny = ");
	NYBox = new QSpinBox(this);
	createSBox(NYLabel, NYBox, glwindow->NY, 2, 500, h);
	QObject::connect(NYBox, QOverload<int>::of(&QSpinBox::valueChanged),
	[=](int i){ 
		glwindow->NY = i;
		glwindow->initXF();
	});
	++h;


	KLabel = new QLabel("k = ");
	KBox = new QSpinBox(this);
	createSBox(KLabel, KBox, glwindow->curFunction, 0, functionsCount - 1, h);
	QObject::connect(KBox, QOverload<int>::of(&QSpinBox::valueChanged),
	[=](int i){ 
		glwindow->curFunction = i;
		glwindow->initXF();
	});
	++h;

	//_______________________________

	AXLabel = new QLabel("ax = ");
	AXBox = new QDoubleSpinBox(this);
	createDSBox(AXLabel, AXBox, glwindow->AX, h);
	QObject::connect(AXBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	[=](double d){ 
		glwindow->AX = d;
		glwindow->initXF();
	});
	++h;

	AYLabel = new QLabel("ay = ");
	AYBox = new QDoubleSpinBox(this);
	createDSBox(AYLabel, AYBox, glwindow->AY, h);
	QObject::connect(AYBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	[=](double d){ 
		glwindow->AY = d;
		glwindow->initXF();
	});
	++h;

	BXLabel = new QLabel("bx = ");
	BXBox = new QDoubleSpinBox(this);
	createDSBox(BXLabel, BXBox, glwindow->BX, h);
	QObject::connect(BXBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	[=](double d){ 
		glwindow->BX = d;
		glwindow->initXF();
	});
	++h;

	BYLabel = new QLabel("by = ");
	BYBox = new QDoubleSpinBox(this);
	createDSBox(BYLabel, BYBox, glwindow->BY, h);
	QObject::connect(BYBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
	[=](double d){ 
		glwindow->BY = d;
		glwindow->initXF();
	});
	++h;

	//__________________________________

	SpanLabel = new QLabel("span = ");
	SpanBox = new QSpinBox(this);
	createSBox(SpanLabel, SpanBox, 0, -10, 10, h);
	QObject::connect(SpanBox, QOverload<int>::of(&QSpinBox::valueChanged),
	[=](int i){ glwindow->span = i; });
	++h;

	PLabel = new QLabel("p = ");
	PBox = new QSpinBox(this);
	createSBox(PLabel, PBox, 0, -10, 10, h);
	QObject::connect(PBox, QOverload<int>::of(&QSpinBox::valueChanged),
	[=](int i){ 
		glwindow->p = i;
		glwindow->initXF();
	});
	++h;

	RotLabel = new QLabel("rot = ");
	RotBox = new QSpinBox(this);
	createSBox(RotLabel, RotBox, 9, 0, maxRot - 1, h);
	QObject::connect(RotBox, QOverload<int>::of(&QSpinBox::valueChanged),
	[=](int i){ 
		glwindow->curRot = i;
		// glwindow->initXF();
	});
	++h;

	//__________________________________

	QPushButton * draw = new QPushButton("draw", this);
	QObject::connect(draw, SIGNAL(clicked()), this, SLOT(upd()));
	grid->addWidget(draw, h, 0, 1, 2);
	++h;

	MethodLabel = new QLabel(methods[glwindow->curMethod]);
	MethodLabel->setStyleSheet("color: rgba(255,255,255)");
	grid->addWidget(MethodLabel, h, 0, 1, 3);
	++h;

	FunctionLabel = new QLabel(functions[glwindow->curFunction]);
	FunctionLabel->setStyleSheet("color: rgba(255,255,255)");
	grid->addWidget(FunctionLabel, h, 0, 1, 3);
	++h;

	//__________________________________

	MinLabel = new QLabel(QString("FMin = %1").arg(glwindow->FMin));
	MinLabel->setStyleSheet("color: rgba(255,255,255)");
	grid->addWidget(MinLabel, h, 0, 1, 3);
	++h;

	MaxLabel = new QLabel(QString("FMax = %1").arg(glwindow->FMax));
	MaxLabel->setStyleSheet("color: rgba(255,255,255)");
	grid->addWidget(MaxLabel, h, 0, 1, 3);
	++h;

	// glwindow->glClearColor(1.0, 1.0, 1.0, 1.0);

	// glwindow->setDefaultCamera();
}


//=================================================
//=================================================

void Window::createSBox(QLabel *label, QSpinBox *sbox, int val, int a, int b, int h)
{
	label->setStyleSheet("color: rgba(255,255,255)");
	grid->addWidget(label, h, 0, 1, 1);

	sbox->setRange(a, b);
	sbox->setValue(val);
	sbox->setSingleStep(1);
	grid->addWidget(sbox, h, 1, 1, 2);
}

void Window::createDSBox(QLabel *label, QDoubleSpinBox *dsbox, double val, int h)
{
	label->setStyleSheet("color: rgba(255,255,255)");
	grid->addWidget(label, h, 0, 1, 1);

	dsbox->setRange(-100, 100);
	dsbox->setDecimals(5);
	dsbox->setSingleStep(0.01);
	dsbox->setValue(val);
	grid->addWidget(dsbox, h, 1, 1, 2);
}

//=================================================
//=================================================

void Window::keyPressEvent(QKeyEvent *event) 
{
	switch (event->key()) {
		case Qt::Key_Escape:
			qApp->quit();
			break;
		case Qt::Key_0:
			glwindow->curFunction = (glwindow->curFunction + 1) % functionsCount;
			KBox->setValue(glwindow->curFunction);
			this->updateFunctionLabel();
			break;
		case Qt::Key_1:
			glwindow->curMethod = (glwindow->curMethod + 1) % methodsCount;
			this->updateMethodLabel();
			break;
		case Qt::Key_2:
			SpanBox->setValue(glwindow->span + 1);
			break;
		case Qt::Key_3:
			SpanBox->setValue(glwindow->span - 1);
			break;
		case Qt::Key_4:
			NXBox->setValue(glwindow->NX*2);
			NYBox->setValue(glwindow->NY*2);
			break;
		case Qt::Key_5:
			NXBox->setValue(glwindow->NX/2);
			NYBox->setValue(glwindow->NY/2);
			break;
		case Qt::Key_6:
			PBox->setValue(glwindow->p + 1);
			break;
		case Qt::Key_7:
			PBox->setValue(glwindow->p - 1);
			break;
		case Qt::Key_C:
			glwindow->setDefaultCamera();
			break;
		case Qt::Key_W:
			glwindow->angle_v = MIN(glwindow->angle_v + 5.0, 80);
			break;
		case Qt::Key_S:
			glwindow->angle_v = MAX(glwindow->angle_v - 5.0, -80);
			break;
		case Qt::Key_A:
			glwindow->curRot = (glwindow->curRot + 1) % maxRot;
			RotBox->setValue(glwindow->curRot);
			break;
		case Qt::Key_D:
			glwindow->curRot = (glwindow->curRot + glwindow->maxRot - 1) % maxRot;
			RotBox->setValue(glwindow->curRot);
			break;
		case Qt::Key_Plus:
			glwindow->camera_p = MAX(glwindow->camera_p - 0.1, -0.5);
			break;
		case Qt::Key_Minus:
			glwindow->camera_p = MIN(glwindow->camera_p + 0.1, 0.5);
			break;
	} 

	this->updateMinMaxLabel();
	glwindow->updateGL();
}

void Window::resizeEvent(QResizeEvent *event)
{
	// size = QSize(event->size().width(), event->size().height());

	// toolbar->resize(toolbar_width, size.height());
	// toolbar->move(QPoint(size.width() - toolbar_width,0));

	// glwindow->resize(size.width() - toolbar_width, size.height());
	// glwindow->move(QPoint(0,0));
	size = QSize(event->size().width() + 700, event->size().height());

	toolbar->resize(toolbar_width, event->size().height());
	toolbar->move(QPoint(0,0));

	glwindow->resize(700, size.height());
	glwindow->move(QPoint(toolbar_width,0));
}
//==============================================
//==============================================
//==============================================

void GLWindow::initializeGL()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);

	setDefaultCamera();
}

void GLWindow::paintGL()
{
	setProjectionMatrix();

	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);

	glBegin(GL_QUADS);

	if(curMethod == 0 || curMethod == 2)
		drawMethod1();

	if(curMethod == 3 || curMethod == 2)
		drawFunction();

	glEnd();

	glDisable(GL_DEPTH_TEST);
}

void GLWindow::resizeGL(int nWidth, int nHeight)
{
	glViewport(0, 0, nWidth, nHeight);
	updateGL();
}

// #define ANGLE_DIFF	(5)
// #define POSITION_DIFF	(0.1)

void GLWindow::setProjectionMatrix()
{
	glViewport(0, 0, width(), height());

    projection.setToIdentity();
    projection.perspective(45.0, width()/height(), 2, 100);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf((const GLfloat*)projection.data());
    // glLoadIdentity();

    view.setToIdentity();
    view.translate(0,0,-(1 + camera_p)*(FMax-FMin + BX-AX + BY-AY));
    view.rotate(270.0, QVector3D(1, 0, 0));

    view.rotate(angle_v, QVector3D(1, 0, 0));
    view.rotate(curRot * (360.0 / maxRot), QVector3D(0, 0, 1));

    view.translate(-(BX+AX)/2, -(BY+AY)/2, -(FMax+FMin)/2);

	double scale = pow(2, span);

    model.setToIdentity();
    model.scale(scale, scale, 1.0f);
    model =  view * model;

    // glMatrixMode(GL_MODELVIEW);
    // glLoadMatrixf((const GLfloat*)view.data());
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf((const GLfloat*)model.data());
}

void GLWindow::setDefaultCamera()
{
	camera_p = 0;
	angle_v = 0;
}

//==============================================

void GLWindow::drawFunction()
{
	glColor3d(1.0,0.0,0.0);
	
	double	x, y, z;

	int npX = max(width()/NX, 2);
	int npY = max(width()/NY, 2);

	double dx = (BX - AX) / ((NX - 1)*(npX + 1)),
		dy = (BY - AY) / ((NY - 1)*(npY + 1));

	for (int i = 0; i < (NX - 1)*(npX + 1); i++)
		for (int j = 0; j < (NY - 1)*(npY + 1); j++) {

			glColor3d(1.0 * ((NX - 1)*(npX + 1) - i) / ((NX - 1)*(npX + 1)), 1.0 * j / ((NY - 1)*(npY + 1)), 0.0);
			x = dx * i + AX;
			y = dy * j + AY;
			z = funcarray[curFunction](x, y);
			glVertex3d(x, y, z);
			x = dx * (i + 1) + AX;
			y = dy * j + AY;
			z = funcarray[curFunction](x, y);
			glVertex3d(x, y, z);
			x = dx * (i + 1) + AX;
			y = dy * (j + 1) + AY;
			z = funcarray[curFunction](x, y);
			glVertex3d(x, y, z);
			x = dx * i + AX;
			y = dy * (j + 1) + AY;
			z = funcarray[curFunction](x, y);
			glVertex3d(x, y, z);
		}
}

void GLWindow::drawMethod1()
{
	glColor3d(1.0,0.0,0.0);
	
	double	x, y, z;

	int npX = max(width()/NX/10, 2);
	int npY = max(width()/NY/10, 2);

	double dx = (BX - AX) / ((NX - 1)*(npX - 1)),
		dy = (BY - AY) / ((NY - 1)*(npY - 1));

	for (int i = 0; i < NX - 1; i++)
		for (int j = 0; j < NY - 1; j++)
			for (int ik = 0; ik < npX - 1; ik++)
				for (int jk = 0; jk < npY - 1; jk++) {
					glColor3d(0.0, 1.0 * (NX - i) * (npX - ik) / NX / npX, 1.0 * j/ NY);
					x = X[i] + dx * ik;
					y = Y[j] + dy * jk;
					z = Pf(i, j, NX, x, y, X, Y, coeffs1);
					glVertex3d(x, y, z);
					x = X[i] + dx * (ik + 1);
					y = Y[j] + dy * jk;
					z = Pf(i, j, NX, x, y, X, Y, coeffs1);
					glVertex3d(x, y, z);
					x = X[i] + dx * (ik + 1);
					y = Y[j] + dy * (jk + 1);
					z = Pf(i, j, NX, x, y, X, Y, coeffs1);
					glVertex3d(x, y, z);
					x = X[i] + dx * ik;
					y = Y[j] + dy * (jk + 1);
					z = Pf(i, j, NX, x, y, X, Y, coeffs1);
					glVertex3d(x, y, z);
				}
}
