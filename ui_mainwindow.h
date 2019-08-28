/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Thu 22. Aug 11:23:50 2019
**      by: Qt User Interface Compiler version 4.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>
#include <glrail3d.h>
#include <glsection.h>
#include <qcustomplot.h>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionOpen_Rail;
    QAction *actionClose;
    QAction *actionClose_Rail;
    QAction *actionShow_Curves;
    QAction *actionAbout;
    QAction *actionCreate_Rail;
    QWidget *centralWidget;
    QGroupBox *sectionGB;
    GLSection *railSec;
    QGroupBox *rail3dGB;
    QGroupBox *modifySecGB;
    QSpinBox *secNumbSp;
    QDoubleSpinBox *spacingSp;
    QLabel *label;
    QLabel *label_2;
    QPushButton *zoomIn;
    GLRail3d *rail3D;
    QPushButton *zoomOut;
    QPushButton *animatePB;
    QGroupBox *curvesGB;
    QCustomPlot *dispCurvs;
    QLabel *curvesInit;
    QGroupBox *meshRefGB;
    QPushButton *pushButton;
    QGroupBox *matPropGB;
    QLineEdit *densityText;
    QLabel *label_3;
    QLineEdit *longVelocText;
    QLabel *label_4;
    QLineEdit *shearVelText;
    QLabel *label_5;
    QLabel *matNameLab;
    QComboBox *matNameCBox;
    QGroupBox *graphPropGB;
    QLabel *label_6;
    QSpinBox *maxWaveSp;
    QLabel *label_8;
    QSpinBox *curvNumbSp;
    QPushButton *computeButt;
    QGroupBox *solTypeGB;
    QRadioButton *rbLinear;
    QRadioButton *rbQuadratic;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuDisplay;
    QMenu *menuHelp;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(871, 671);
        MainWindow->setMaximumSize(QSize(871, 671));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/images.png"), QSize(), QIcon::Normal, QIcon::Off);
        MainWindow->setWindowIcon(icon);
        actionOpen_Rail = new QAction(MainWindow);
        actionOpen_Rail->setObjectName(QString::fromUtf8("actionOpen_Rail"));
        actionClose = new QAction(MainWindow);
        actionClose->setObjectName(QString::fromUtf8("actionClose"));
        actionClose_Rail = new QAction(MainWindow);
        actionClose_Rail->setObjectName(QString::fromUtf8("actionClose_Rail"));
        actionShow_Curves = new QAction(MainWindow);
        actionShow_Curves->setObjectName(QString::fromUtf8("actionShow_Curves"));
        actionShow_Curves->setCheckable(true);
        actionShow_Curves->setChecked(true);
        actionShow_Curves->setEnabled(false);
        actionAbout = new QAction(MainWindow);
        actionAbout->setObjectName(QString::fromUtf8("actionAbout"));
        actionCreate_Rail = new QAction(MainWindow);
        actionCreate_Rail->setObjectName(QString::fromUtf8("actionCreate_Rail"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        sectionGB = new QGroupBox(centralWidget);
        sectionGB->setObjectName(QString::fromUtf8("sectionGB"));
        sectionGB->setGeometry(QRect(170, 320, 241, 291));
        QFont font;
        font.setPointSize(10);
        font.setBold(true);
        font.setWeight(75);
        sectionGB->setFont(font);
        railSec = new GLSection(sectionGB);
        railSec->setObjectName(QString::fromUtf8("railSec"));
        railSec->setGeometry(QRect(10, 20, 221, 261));
        rail3dGB = new QGroupBox(centralWidget);
        rail3dGB->setObjectName(QString::fromUtf8("rail3dGB"));
        rail3dGB->setGeometry(QRect(420, 320, 441, 291));
        rail3dGB->setFont(font);
        modifySecGB = new QGroupBox(rail3dGB);
        modifySecGB->setObjectName(QString::fromUtf8("modifySecGB"));
        modifySecGB->setGeometry(QRect(327, 20, 111, 141));
        QFont font1;
        font1.setPointSize(9);
        modifySecGB->setFont(font1);
        secNumbSp = new QSpinBox(modifySecGB);
        secNumbSp->setObjectName(QString::fromUtf8("secNumbSp"));
        secNumbSp->setGeometry(QRect(10, 40, 91, 22));
        secNumbSp->setMinimum(1);
        secNumbSp->setValue(10);
        spacingSp = new QDoubleSpinBox(modifySecGB);
        spacingSp->setObjectName(QString::fromUtf8("spacingSp"));
        spacingSp->setGeometry(QRect(10, 90, 91, 22));
        spacingSp->setDecimals(3);
        spacingSp->setMaximum(9.999);
        spacingSp->setSingleStep(0.001);
        spacingSp->setValue(0.02);
        label = new QLabel(modifySecGB);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(10, 20, 101, 16));
        QFont font2;
        font2.setPointSize(8);
        label->setFont(font2);
        label_2 = new QLabel(modifySecGB);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(10, 70, 51, 16));
        label_2->setFont(font2);
        zoomIn = new QPushButton(rail3dGB);
        zoomIn->setObjectName(QString::fromUtf8("zoomIn"));
        zoomIn->setGeometry(QRect(340, 170, 91, 23));
        rail3D = new GLRail3d(rail3dGB);
        rail3D->setObjectName(QString::fromUtf8("rail3D"));
        rail3D->setGeometry(QRect(10, 20, 311, 261));
        zoomOut = new QPushButton(rail3dGB);
        zoomOut->setObjectName(QString::fromUtf8("zoomOut"));
        zoomOut->setGeometry(QRect(340, 210, 91, 23));
        animatePB = new QPushButton(rail3dGB);
        animatePB->setObjectName(QString::fromUtf8("animatePB"));
        animatePB->setEnabled(false);
        animatePB->setGeometry(QRect(325, 247, 111, 31));
        curvesGB = new QGroupBox(centralWidget);
        curvesGB->setObjectName(QString::fromUtf8("curvesGB"));
        curvesGB->setGeometry(QRect(170, 0, 581, 311));
        curvesGB->setFont(font);
        dispCurvs = new QCustomPlot(curvesGB);
        dispCurvs->setObjectName(QString::fromUtf8("dispCurvs"));
        dispCurvs->setGeometry(QRect(10, 20, 561, 281));
        curvesInit = new QLabel(dispCurvs);
        curvesInit->setObjectName(QString::fromUtf8("curvesInit"));
        curvesInit->setGeometry(QRect(0, 0, 561, 281));
        meshRefGB = new QGroupBox(centralWidget);
        meshRefGB->setObjectName(QString::fromUtf8("meshRefGB"));
        meshRefGB->setGeometry(QRect(10, 320, 151, 81));
        meshRefGB->setFont(font);
        pushButton = new QPushButton(meshRefGB);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));
        pushButton->setGeometry(QRect(24, 30, 81, 31));
        matPropGB = new QGroupBox(centralWidget);
        matPropGB->setObjectName(QString::fromUtf8("matPropGB"));
        matPropGB->setGeometry(QRect(10, 430, 151, 181));
        matPropGB->setFont(font);
        densityText = new QLineEdit(matPropGB);
        densityText->setObjectName(QString::fromUtf8("densityText"));
        densityText->setEnabled(false);
        densityText->setGeometry(QRect(101, 59, 41, 20));
        densityText->setFont(font2);
        label_3 = new QLabel(matPropGB);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(7, 60, 81, 20));
        QFont font3;
        font3.setPointSize(8);
        font3.setBold(false);
        font3.setWeight(50);
        label_3->setFont(font3);
        longVelocText = new QLineEdit(matPropGB);
        longVelocText->setObjectName(QString::fromUtf8("longVelocText"));
        longVelocText->setEnabled(false);
        longVelocText->setGeometry(QRect(101, 101, 41, 20));
        longVelocText->setFont(font2);
        label_4 = new QLabel(matPropGB);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(6, 102, 91, 20));
        label_4->setFont(font3);
        shearVelText = new QLineEdit(matPropGB);
        shearVelText->setObjectName(QString::fromUtf8("shearVelText"));
        shearVelText->setEnabled(false);
        shearVelText->setGeometry(QRect(100, 145, 41, 20));
        shearVelText->setFont(font2);
        label_5 = new QLabel(matPropGB);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(5, 145, 91, 20));
        label_5->setFont(font3);
        matNameLab = new QLabel(matPropGB);
        matNameLab->setObjectName(QString::fromUtf8("matNameLab"));
        matNameLab->setGeometry(QRect(8, 26, 71, 16));
        matNameLab->setFont(font3);
        matNameCBox = new QComboBox(matPropGB);
        matNameCBox->setObjectName(QString::fromUtf8("matNameCBox"));
        matNameCBox->setEnabled(false);
        matNameCBox->setGeometry(QRect(88, 23, 61, 22));
        matNameCBox->setEditable(false);
        graphPropGB = new QGroupBox(centralWidget);
        graphPropGB->setObjectName(QString::fromUtf8("graphPropGB"));
        graphPropGB->setGeometry(QRect(10, 7, 151, 301));
        QFont font4;
        font4.setPointSize(9);
        font4.setBold(true);
        font4.setWeight(75);
        graphPropGB->setFont(font4);
        label_6 = new QLabel(graphPropGB);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(22, 23, 121, 16));
        label_6->setFont(font4);
        maxWaveSp = new QSpinBox(graphPropGB);
        maxWaveSp->setObjectName(QString::fromUtf8("maxWaveSp"));
        maxWaveSp->setEnabled(true);
        maxWaveSp->setGeometry(QRect(25, 46, 81, 22));
        QFont font5;
        font5.setPointSize(8);
        font5.setBold(true);
        font5.setWeight(75);
        maxWaveSp->setFont(font5);
        maxWaveSp->setMinimum(1);
        maxWaveSp->setMaximum(1000);
        maxWaveSp->setValue(200);
        label_8 = new QLabel(graphPropGB);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(23, 100, 101, 16));
        label_8->setFont(font4);
        curvNumbSp = new QSpinBox(graphPropGB);
        curvNumbSp->setObjectName(QString::fromUtf8("curvNumbSp"));
        curvNumbSp->setGeometry(QRect(27, 124, 81, 22));
        curvNumbSp->setFont(font2);
        curvNumbSp->setMinimum(5);
        curvNumbSp->setMaximum(30);
        curvNumbSp->setValue(20);
        computeButt = new QPushButton(graphPropGB);
        computeButt->setObjectName(QString::fromUtf8("computeButt"));
        computeButt->setGeometry(QRect(10, 264, 131, 31));
        solTypeGB = new QGroupBox(graphPropGB);
        solTypeGB->setObjectName(QString::fromUtf8("solTypeGB"));
        solTypeGB->setEnabled(true);
        solTypeGB->setGeometry(QRect(10, 178, 131, 71));
        rbLinear = new QRadioButton(solTypeGB);
        rbLinear->setObjectName(QString::fromUtf8("rbLinear"));
        rbLinear->setGeometry(QRect(10, 20, 82, 17));
        rbLinear->setFont(font2);
        rbLinear->setChecked(true);
        rbQuadratic = new QRadioButton(solTypeGB);
        rbQuadratic->setObjectName(QString::fromUtf8("rbQuadratic"));
        rbQuadratic->setGeometry(QRect(10, 50, 82, 17));
        rbQuadratic->setFont(font2);
        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 871, 20));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuDisplay = new QMenu(menuBar);
        menuDisplay->setObjectName(QString::fromUtf8("menuDisplay"));
        menuDisplay->setEnabled(false);
        menuHelp = new QMenu(menuBar);
        menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuDisplay->menuAction());
        menuBar->addAction(menuHelp->menuAction());
        menuFile->addAction(actionOpen_Rail);
        menuFile->addAction(actionCreate_Rail);
        menuFile->addAction(actionClose_Rail);
        menuFile->addAction(actionClose);
        menuDisplay->addAction(actionShow_Curves);
        menuHelp->addAction(actionAbout);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "RailWave", 0, QApplication::UnicodeUTF8));
        actionOpen_Rail->setText(QApplication::translate("MainWindow", "Open Rail", 0, QApplication::UnicodeUTF8));
        actionClose->setText(QApplication::translate("MainWindow", "Close", 0, QApplication::UnicodeUTF8));
        actionClose_Rail->setText(QApplication::translate("MainWindow", "Close Rail", 0, QApplication::UnicodeUTF8));
        actionShow_Curves->setText(QApplication::translate("MainWindow", "Show Curves", 0, QApplication::UnicodeUTF8));
        actionAbout->setText(QApplication::translate("MainWindow", "About", 0, QApplication::UnicodeUTF8));
        actionCreate_Rail->setText(QApplication::translate("MainWindow", "Create Rail", 0, QApplication::UnicodeUTF8));
        sectionGB->setTitle(QApplication::translate("MainWindow", "Section", 0, QApplication::UnicodeUTF8));
        rail3dGB->setTitle(QApplication::translate("MainWindow", "Rail 3D", 0, QApplication::UnicodeUTF8));
        modifySecGB->setTitle(QApplication::translate("MainWindow", "Modify Section", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Sections Number", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "Spacing", 0, QApplication::UnicodeUTF8));
        zoomIn->setText(QApplication::translate("MainWindow", "Zoom In", 0, QApplication::UnicodeUTF8));
        zoomOut->setText(QApplication::translate("MainWindow", "Zoom Out", 0, QApplication::UnicodeUTF8));
        animatePB->setText(QApplication::translate("MainWindow", "Animate", 0, QApplication::UnicodeUTF8));
        curvesGB->setTitle(QApplication::translate("MainWindow", "Dispersion Curves", 0, QApplication::UnicodeUTF8));
        curvesInit->setText(QString());
        meshRefGB->setTitle(QApplication::translate("MainWindow", "Mesh Refinement", 0, QApplication::UnicodeUTF8));
        pushButton->setText(QApplication::translate("MainWindow", "Refine", 0, QApplication::UnicodeUTF8));
        matPropGB->setTitle(QApplication::translate("MainWindow", "Material Properties", 0, QApplication::UnicodeUTF8));
        densityText->setText(QApplication::translate("MainWindow", "7850", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "Density[Kg/m3]", 0, QApplication::UnicodeUTF8));
        longVelocText->setText(QApplication::translate("MainWindow", "5996", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("MainWindow", "Long Velocity[m/s]", 0, QApplication::UnicodeUTF8));
        shearVelText->setText(QApplication::translate("MainWindow", "3260", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("MainWindow", "Shear Velocity[m/s]", 0, QApplication::UnicodeUTF8));
        matNameLab->setText(QApplication::translate("MainWindow", "Material Name", 0, QApplication::UnicodeUTF8));
        matNameCBox->clear();
        matNameCBox->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Steel", 0, QApplication::UnicodeUTF8)
        );
        graphPropGB->setTitle(QApplication::translate("MainWindow", "Graph Properties", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("MainWindow", "Max Wave Number", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("MainWindow", "Curves Number", 0, QApplication::UnicodeUTF8));
        computeButt->setText(QApplication::translate("MainWindow", "Compute", 0, QApplication::UnicodeUTF8));
        solTypeGB->setTitle(QApplication::translate("MainWindow", "Solution Type", 0, QApplication::UnicodeUTF8));
        rbLinear->setText(QApplication::translate("MainWindow", "Linear", 0, QApplication::UnicodeUTF8));
        rbQuadratic->setText(QApplication::translate("MainWindow", "Quadratic", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuDisplay->setTitle(QApplication::translate("MainWindow", "Display", 0, QApplication::UnicodeUTF8));
        menuHelp->setTitle(QApplication::translate("MainWindow", "Help", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
