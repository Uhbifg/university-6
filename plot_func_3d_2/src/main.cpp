#include"window.h"

#include <qapplication.h>
#include <qgl.h>

int main(int argc, char *argv[])
{
	double ax, ay, bx, by;
	int nx, ny, k;

	if (argc != 8)
	{
		printf("%d\n", argc);
		printf("wrong number of arguments\n");
		return -1;
	}
	printf("@\n");

	ax = atof(argv[1]);
	bx = atof(argv[2]);
	ay = atof(argv[3]);
	by = atof(argv[4]);
	nx = atoi(argv[5]);
	ny = atoi(argv[6]);
	k = atoi(argv[7]);

	if(ax >= bx || ay >= by || k < 0 || nx < 0 || ny < 0 || k > 7)
	{
		printf("wrong arguments\n");
		return -1;
	}

	QApplication::setApplicationName("Function3D");
	QApplication app(argc, argv);

	Window my_window(ax, ay, bx, by, nx, ny, k);
	
	// app.setActiveWindow(&my_window);
	my_window.show();

	app.exec();

	return 0;
}



