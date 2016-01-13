#include <QFileDialog>
#include <QFile>
#include <QMessageBox>
#include <QTextStream>
#include <QString>
#include "visualization.h"
#include "ui_visualization.h"

using namespace std;

Visualization::Visualization(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Visualization)
{
    ui->setupUi(this);
    finalWidget = new QtGnuplotWidget();

    this->widgetPoints1 = new QtGnuplotWidget();        //different widgets because different output
    this->widgetPoints1->installEventFilter(this);
    this->widgetPoints1->setStatusLabelActive(true);

    this->widgetPoints2 = new QtGnuplotWidget();
    this->widgetPoints2->installEventFilter(this);
    this->widgetPoints2->setStatusLabelActive(true);

    instance = new QtGnuplotInstance();
    instance->setWidget(widgetPoints1);
    ui->selectFile->setEnabled(false);
    ui->InsertPoints->setEnabled(false);
    ui->solve->setEnabled(false);
    ui->InsertPoints2->setEnabled(false);
}

Visualization::~Visualization()
{
    delete ui;
    delete instance;
    delete finalWidget;
    delete widgetPoints1;
    delete widgetPoints2;
}

void Visualization::on_actionExit_triggered()
{
    qApp->quit();
}

void Visualization::on_actionOpen_triggered()
{
}

bool Visualization::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::MouseButtonPress)
    {
        if (obj == this->widgetPoints1) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent *>(event);
            if (mouseEvent->button() == Qt::LeftButton) {
                ui->pointsTxt->append(this->widgetPoints1->getStatusLabel()->text().replace(QRegExp("(, |,  )"), " ").replace(QRegExp("\n "), "\n").replace(0, 1, ""));
            }
        }
        if (obj == this->widgetPoints2) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent *>(event);
            if (mouseEvent->button() == Qt::LeftButton) {
                ui->points2Txt->append(this->widgetPoints2->getStatusLabel()->text().replace(QRegExp("(, |,  )"), " ").replace(QRegExp("\n "), "\n").replace(0, 1, ""));
            }
        }
    }
    return QObject::eventFilter(obj, event);
}






void Visualization::on_readFromFile_clicked()
{
    ui->selectFile->setEnabled(true);
    ui->InsertPoints->setEnabled(false);
    ui->InsertPoints2->setEnabled(false);
}

void Visualization::on_selectFile_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QString(),
                tr("All Files (*.*)"));

    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly)) {
            QMessageBox::critical(this, tr("Error"), tr("Could not open file"));
            return;
        }
        QTextStream in(&file);
        QString str = in.readAll();
        ui->equationsTxt->setText(str);
        file.close();
        this->ui->solve->setEnabled(true);
    }
}

void Visualization::on_fromPoints_clicked()
{
        ui->selectFile->setEnabled(false);
        ui->InsertPoints->setEnabled(true);
        ui->InsertPoints2->setEnabled(true);
}

void Visualization::on_solve_clicked()
{
    if(ui->readFromFile->isChecked())
    {
        QString program = "../exe";
        QStringList arguments;
        arguments << "-read" <<"-d1"<<"3"<<"-d2"<<"3"<<"-solve"<<"8";
        QProcess * executable = new QProcess();
        QString myQString = this->ui->equationsTxt->toPlainText();
        executable->start(program, arguments);
        QByteArray toChar = myQString.toLatin1();
        const char *str = toChar.data();
        int size = toChar.length();
        executable->write(str, size);
        executable->closeWriteChannel();
        if(!executable->waitForFinished()) // beware the timeout default parameter
            qDebug() << "executing program failed with exit code" << executable->exitCode();
        else
            ui->outputTxt->setText(QString(executable->readAllStandardOutput()));
        executable->close();
    }
    if(ui->fromPoints->isChecked())
    {
        QString program = "../exe";
        QStringList arguments;
        int d1 = ui->d1Box->text().toInt();
        int d2 = ui->d2Box->text().toInt();
        if(d1 < 0 || d2 < 0)
        {
            QMessageBox messageBox;
            messageBox.critical(0,"Error","d1 and d2 must be positive!");
            messageBox.setFixedSize(500,200);
        }
        arguments << "-points" <<"-d1"<<ui->d1Box->text()<<"-d2"<<ui->d2Box->text()<<"-solve"<<"8";
        QProcess * executable = new QProcess();
        QString points = this->ui->pointsTxt->toPlainText();
        points.append("\n-d \n");
        points.append(this->ui->points2Txt->toPlainText());
        points.append("\n");
        QByteArray toChar = points.toLatin1();
        const char *str = toChar.data();
        int size = toChar.length();

//        QByteArray toChar2 = points2.toLatin1();
//        const char *str2 = toChar2.data();
//        int size2 = toChar2.length();

        executable->start(program, arguments);

        executable->waitForStarted();

        if(executable->write(str, size) < size)
        {
            QMessageBox::critical(this, tr("Error"), tr("error in writing to process " + executable->readAllStandardOutput()));
            executable->terminate();
            return;
        }
        executable->closeWriteChannel();
        if(!executable->waitForFinished()) // beware the timeout default parameter
        {
            QMessageBox::critical(this, tr("Error"), tr("executing program failed with exit code" + executable->exitCode() + executable->readAllStandardOutput()));
            return;
        }
        else
        {
            ui->outputTxt->setText(QString(executable->readAllStandardOutput()));
            executable->close();
            QFile file("InterpolationEquations.txt");
            if (!file.open(QIODevice::ReadOnly)) {
                QMessageBox::critical(this, tr("Error"), tr("Could not open file"));
                return;
            }
            QTextStream in(&file);
            QString str = in.readAll();
            if(str == "")
            {
                QMessageBox::critical(this, tr("Error"), tr("Interpolation did not find the equations"));
                return;
            }
            ui->equationsTxt->setText(str);
        }
    }
    instance->setWidget(finalWidget);
    finalWidget->show();
    finalWidget->resize(QSize(800,600));
    QStringList functions = this->ui->equationsTxt->toPlainText().split("\n");
    QString f1 = functions.at(0);
    f1 = f1.replace(QString("^"),QString("**"));
    QString f2 = functions.at(1);
    f2 = f2.replace(QString("^"),QString("**"));

    QFile file("solutions.txt");
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::critical(this, tr("Error"), tr("Could not open file"));
        return;
    }

    QTextStream in(&file);
    QString str = in.readAll();
    if(str == "")
    {
        QMessageBox::critical(this, tr("Error"), tr("The system has no solution"));
        return;
    }

    QStringList solutions = str.split("\n");
    double xMin, xMax, yMin, yMax;
    findMinAndMaxXandY(solutions, xMin, xMax, yMin, yMax);
    *instance <<\
    "set yrange ["+QString::number(yMin)+":"+QString::number(yMax)+"]\nset xrange["+QString::number(xMin)+":"+QString::number(xMax)+"]\nset isosamples 500,500\nf(x,y)="+f1+
                 "\nf2(x,y)= "+f2+"\nset contour\nset cntrparam levels discrete 0\nset view 0,0\nunset ztics\nunset surface\nsplot f(x,y), f2(x,y)\n";
}

void Visualization::findMinAndMaxXandY(QStringList & solutions,
                        double & xMin, double & xMax,
                        double & yMin, double & yMax)
{
    xMin = solutions.at(0).split(" ").at(0).toDouble();
    xMax = xMin;
    yMin = solutions.at(0).split(" ").at(1).toDouble();
    yMax = yMin;
    double xCurr, yCurr;
    for(int i = 1; i < solutions.length() - 1 ; i++)        //length - 1 because last one doesnt count
    {
        xCurr = solutions.at(i).split(" ").at(0).toDouble();
        yCurr = solutions.at(i).split(" ").at(1).toDouble();

        if(xCurr > xMax)
            xMax = xCurr;
        if(xCurr < xMin)
            xMin = xCurr;

        if(yCurr < yMin)
            yMin = yCurr;
        if(yCurr > yMax)
            yMax = yCurr;

    }
}

void Visualization::on_InsertPoints_clicked()
{
    instance->setWidget(widgetPoints1);
    this->widgetPoints1->show();
    widgetPoints1->resize(QSize(800,600));
    ui->pointsTxt->setText("");
    *instance <<\
    "set yrange [-30:30]\nset xrange [-30:30]\nset isosamples 500,500\nset contour\nset cntrparam levels discrete 0\nset view 0,0\nunset ztics\nunset surface\nplot 1/0\n";
    ui->solve->setEnabled(true);
}


void Visualization::on_InsertPoints2_clicked()
{
    instance->setWidget(widgetPoints2);
    this->widgetPoints2->show();
    widgetPoints2->resize(QSize(800,600));
    ui->points2Txt->setText("");
    *instance <<\
    "set yrange [-30:30]\nset xrange [-30:30]\nset isosamples 500,500\nset contour\nset cntrparam levels discrete 0\nset view 0,0\nunset ztics\nunset surface\nplot 1/0\n";
    ui->solve->setEnabled(true);
}
