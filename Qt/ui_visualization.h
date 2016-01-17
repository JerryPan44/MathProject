/********************************************************************************
** Form generated from reading UI file 'visualization.ui'
**
** Created by: Qt User Interface Compiler version 5.2.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_VISUALIZATION_H
#define UI_VISUALIZATION_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Visualization
{
public:
    QAction *actionOpen;
    QAction *actionExit;
    QWidget *centralwidget;
    QTextBrowser *equationsTxt;
    QLabel *equationsLbl;
    QTextBrowser *outputTxt;
    QLabel *outputLbl;
    QPushButton *solve;
    QRadioButton *readFromFile;
    QRadioButton *fromPoints;
    QPushButton *selectFile;
    QPushButton *InsertPoints;
    QTextBrowser *pointsTxt;
    QLabel *label;
    QSpinBox *d1Box;
    QSpinBox *d2Box;
    QLabel *label_2;
    QLabel *label_3;
    QTextBrowser *points2Txt;
    QPushButton *InsertPoints2;
    QLabel *label_4;
    QMenuBar *menubar;
    QMenu *menuFile;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *Visualization)
    {
        if (Visualization->objectName().isEmpty())
            Visualization->setObjectName(QStringLiteral("Visualization"));
        Visualization->resize(1024, 768);
        actionOpen = new QAction(Visualization);
        actionOpen->setObjectName(QStringLiteral("actionOpen"));
        actionOpen->setVisible(true);
        actionExit = new QAction(Visualization);
        actionExit->setObjectName(QStringLiteral("actionExit"));
        centralwidget = new QWidget(Visualization);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        equationsTxt = new QTextBrowser(centralwidget);
        equationsTxt->setObjectName(QStringLiteral("equationsTxt"));
        equationsTxt->setGeometry(QRect(10, 130, 481, 151));
        equationsLbl = new QLabel(centralwidget);
        equationsLbl->setObjectName(QStringLiteral("equationsLbl"));
        equationsLbl->setGeometry(QRect(10, 90, 71, 17));
        outputTxt = new QTextBrowser(centralwidget);
        outputTxt->setObjectName(QStringLiteral("outputTxt"));
        outputTxt->setGeometry(QRect(10, 480, 971, 261));
        outputLbl = new QLabel(centralwidget);
        outputLbl->setObjectName(QStringLiteral("outputLbl"));
        outputLbl->setGeometry(QRect(20, 290, 141, 17));
        solve = new QPushButton(centralwidget);
        solve->setObjectName(QStringLiteral("solve"));
        solve->setGeometry(QRect(380, 420, 211, 23));
        readFromFile = new QRadioButton(centralwidget);
        readFromFile->setObjectName(QStringLiteral("readFromFile"));
        readFromFile->setGeometry(QRect(30, 10, 241, 21));
        fromPoints = new QRadioButton(centralwidget);
        fromPoints->setObjectName(QStringLiteral("fromPoints"));
        fromPoints->setGeometry(QRect(640, 10, 301, 21));
        selectFile = new QPushButton(centralwidget);
        selectFile->setObjectName(QStringLiteral("selectFile"));
        selectFile->setGeometry(QRect(20, 50, 141, 23));
        InsertPoints = new QPushButton(centralwidget);
        InsertPoints->setObjectName(QStringLiteral("InsertPoints"));
        InsertPoints->setGeometry(QRect(580, 80, 121, 23));
        pointsTxt = new QTextBrowser(centralwidget);
        pointsTxt->setObjectName(QStringLiteral("pointsTxt"));
        pointsTxt->setGeometry(QRect(550, 150, 191, 111));
        label = new QLabel(centralwidget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(560, 130, 57, 15));
        d1Box = new QSpinBox(centralwidget);
        d1Box->setObjectName(QStringLiteral("d1Box"));
        d1Box->setGeometry(QRect(330, 350, 91, 24));
        d2Box = new QSpinBox(centralwidget);
        d2Box->setObjectName(QStringLiteral("d2Box"));
        d2Box->setGeometry(QRect(520, 350, 91, 24));
        label_2 = new QLabel(centralwidget);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(330, 320, 81, 16));
        label_3 = new QLabel(centralwidget);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(520, 320, 71, 16));
        points2Txt = new QTextBrowser(centralwidget);
        points2Txt->setObjectName(QStringLiteral("points2Txt"));
        points2Txt->setGeometry(QRect(770, 150, 191, 111));
        InsertPoints2 = new QPushButton(centralwidget);
        InsertPoints2->setObjectName(QStringLiteral("InsertPoints2"));
        InsertPoints2->setGeometry(QRect(790, 90, 131, 23));
        label_4 = new QLabel(centralwidget);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(770, 130, 57, 15));
        Visualization->setCentralWidget(centralwidget);
        menubar = new QMenuBar(Visualization);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 1024, 20));
        menubar->setDefaultUp(true);
        menubar->setNativeMenuBar(false);
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        Visualization->setMenuBar(menubar);
        statusbar = new QStatusBar(Visualization);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        Visualization->setStatusBar(statusbar);

        menubar->addAction(menuFile->menuAction());
        menuFile->addSeparator();
        menuFile->addSeparator();
        menuFile->addAction(actionOpen);
        menuFile->addAction(actionExit);

        retranslateUi(Visualization);

        QMetaObject::connectSlotsByName(Visualization);
    } // setupUi

    void retranslateUi(QMainWindow *Visualization)
    {
        Visualization->setWindowTitle(QApplication::translate("Visualization", "Visualization App", 0));
        actionOpen->setText(QApplication::translate("Visualization", "Open", 0));
        actionOpen->setShortcut(QApplication::translate("Visualization", "Alt+F", 0));
        actionExit->setText(QApplication::translate("Visualization", "Exit", 0));
        actionExit->setShortcut(QApplication::translate("Visualization", "Alt+X", 0));
        equationsLbl->setText(QApplication::translate("Visualization", "Equations", 0));
        outputLbl->setText(QApplication::translate("Visualization", "Program Output", 0));
        solve->setText(QApplication::translate("Visualization", "Solve", 0));
        readFromFile->setText(QApplication::translate("Visualization", "Read system of equations from file", 0));
        fromPoints->setText(QApplication::translate("Visualization", "Construct system of equations from points", 0));
        selectFile->setText(QApplication::translate("Visualization", "select file", 0));
        InsertPoints->setText(QApplication::translate("Visualization", "Insert Points 1", 0));
        label->setText(QApplication::translate("Visualization", "Points", 0));
        label_2->setText(QApplication::translate("Visualization", "degree 1", 0));
        label_3->setText(QApplication::translate("Visualization", "degree 2", 0));
        InsertPoints2->setText(QApplication::translate("Visualization", "Insert Points 2", 0));
        label_4->setText(QApplication::translate("Visualization", "Points 2", 0));
        menuFile->setTitle(QApplication::translate("Visualization", "File", 0));
    } // retranslateUi

};

namespace Ui {
    class Visualization: public Ui_Visualization {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_VISUALIZATION_H
