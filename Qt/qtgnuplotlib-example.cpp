//#include <QApplication>
#include "visualization.h"

int main(int argc, char* argv[])
{
	QApplication app(argc, argv);
    QCoreApplication::setAttribute(Qt::AA_DontUseNativeMenuBar);
    Visualization v;
    v.show();

    return app.exec();
}
