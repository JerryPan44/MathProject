#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <QMainWindow>
#include <QtGnuplot/QtGnuplotWidget.h>
#include <QtGnuplot/QtGnuplotInstance.h>

namespace Ui {
class Visualization;
}

class Visualization : public QMainWindow
{
    Q_OBJECT

public:
    explicit Visualization(QWidget *parent = 0);
    ~Visualization();

private:
    QtGnuplotWidget * widgetPoints1;
    QtGnuplotWidget * widgetPoints2;
    QtGnuplotWidget * finalWidget;
    QtGnuplotInstance *instance;
    QtGnuplotInstance *instance2;
    QtGnuplotInstance *finalInstance;
    QString plotSolutions(QStringList & solutions);
    void findMinAndMaxXandY(QStringList & solutions,
                            double & xMin, double & xMax,
                            double & yMin, double & yMax);
protected:
    bool eventFilter(QObject *obj, QEvent *event);


private slots:
    void on_actionExit_triggered();

    void on_actionOpen_triggered();

    void on_readFromFile_clicked();

    void on_selectFile_clicked();

    void on_fromPoints_clicked();

    void on_solve_clicked();

    void on_InsertPoints_clicked();

    void on_InsertPoints2_clicked();

private:
    Ui::Visualization *ui;
};

#endif // VISUALIZATION_H
