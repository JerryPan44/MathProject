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

    this->widgetPoints1 = new QtGnuplotWidget();
    this->widgetPoints1->installEventFilter(this);
    this->widgetPoints1->setStatusLabelActive(true);

    this->widgetPoints2 = new QtGnuplotWidget();
    this->widgetPoints2->installEventFilter(this);
    this->widgetPoints2->setStatusLabelActive(true);

    instance = new QtGnuplotInstance();
    instance2 = new QtGnuplotInstance();
    finalInstance = new QtGnuplotInstance();
    instance->setWidget(widgetPoints1);
    instance2->setWidget(widgetPoints2);
    finalInstance->setWidget(finalWidget);
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
                ui->pointsTxt->append(this->widgetPoints1->getStatusLabel()->text());
            }
        }
        if (obj == this->widgetPoints2) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent *>(event);
            if (mouseEvent->button() == Qt::LeftButton) {
                ui->points2Txt->append(this->widgetPoints2->getStatusLabel()->text());
            }
        }
    }
    return QObject::eventFilter(obj, event);
}



void Visualization::on_pushButton_clicked()
{
    finalWidget->show();
    finalWidget->resize(QSize(800,600));
    QStringList functions = this->ui->equationsTxt->toPlainText().split("\n");
    QString f1 = functions.at(0);
    f1 = f1.replace(QString("^"),QString("**"));
    QString f2 = functions.at(1);
    f2 = f2.replace(QString("^"),QString("**"));
    *finalInstance <<\
    "set yrange [-1.5:1.5]\nset xrange [-1.5:1.5]\nset isosamples 500,500\nf(x,y)="+f1+
                 "\nf2(x,y)= "+f2+"\nset contour\nset cntrparam levels discrete 0\nset view 0,0\nunset ztics\nunset surface\nsplot f(x,y), f2(x,y)\n";
    //*instance << "set tics scale 0.75\nset xtics 1\nset ytics 1\nset yrange [-10:10]\nset xlabel 'x'\nset ylabel 'y'\nset zeroaxis\nplot \"<echo '1 2'\" notitle\n";
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
        arguments << "-read" <<"-d1"<<"2"<<"-d2"<<"3"<<"-solve"<<"8";
        QProcess * executable = new QProcess();
        QString myQString = this->ui->equationsTxt->toPlainText();
        executable->start(program, arguments);
        QByteArray toChar = myQString.toLatin1();
        const char *str = toChar.data();
        int size = toChar.length();
        executable->write(str, size);
        if(!executable->waitForFinished()) // beware the timeout default parameter
            qDebug() << "executing program failed with exit code" << executable->exitCode();
        else
            ui->outputTxt->setText(QString(executable->readAllStandardOutput()));
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
        points.append("\n");
        QString points2 = this->ui->points2Txt->toPlainText();

        QByteArray toChar = points.toLatin1();
        const char *str = toChar.data();
        int size = toChar.length();

        QByteArray toChar2 = points2.toLatin1();
        const char *str2 = toChar2.data();
        int size2 = toChar2.length();

        executable->start(program, arguments);
        executable->write(str, size);
        executable->write(str2, size2);

        if(!executable->waitForFinished()) // beware the timeout default parameter
            qDebug() << "executing program failed with exit code" << executable->exitCode();
        else
            ui->outputTxt->setText(QString(executable->readAllStandardOutput()));
    }
}

void Visualization::on_InsertPoints_clicked()
{
    this->widgetPoints1->show();
    widgetPoints1->resize(QSize(800,600));

    *instance <<\
    "set yrange [-30:30]\nset xrange [-30:30]\nset isosamples 500,500\nset contour\nset cntrparam levels discrete 0\nset view 0,0\nunset ztics\nunset surface\nplot 1/0\n";
    ui->solve->setEnabled(true);
}


void Visualization::on_InsertPoints2_clicked()
{
    this->widgetPoints2->show();
    widgetPoints2->resize(QSize(800,600));

    *instance2 <<\
    "set yrange [-30:30]\nset xrange [-30:30]\nset isosamples 500,500\nset contour\nset cntrparam levels discrete 0\nset view 0,0\nunset ztics\nunset surface\nplot 1/0\n";
    ui->solve->setEnabled(true);
}
