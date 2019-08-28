/********************************************************************************
** Form generated from reading UI file 'designrail.ui'
**
** Created: Thu 22. Aug 11:23:50 2019
**      by: Qt User Interface Compiler version 4.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DESIGNRAIL_H
#define UI_DESIGNRAIL_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>
#include <gldesign.h>

QT_BEGIN_NAMESPACE

class Ui_DesignRail
{
public:
    QAction *actionSave_Rail;
    QAction *actionClose;
    QWidget *centralwidget;
    QGroupBox *parametersGB;
    QLineEdit *HDlineEdit;
    QLabel *HDlabel;
    QLineEdit *FDlineEdit;
    QLabel *FDlabel;
    QLineEdit *BDlineEdit;
    QLabel *BDlabel;
    QLineEdit *HWlineEdit;
    QLabel *HWlabel;
    QLineEdit *WlineEdit;
    QLabel *Wlabel;
    QLineEdit *BWlineEdit;
    QLabel *BWlabel;
    QLineEdit *R1HlineEdit;
    QLabel *R1Hlabel;
    QLineEdit *R2HlineEdit;
    QLabel *R2Hlabel;
    QLineEdit *RHBlineEdit;
    QLabel *RHBlabel;
    QLineEdit *alphalineEdit;
    QLabel *alphalabel;
    QLineEdit *betalineEdit;
    QLabel *betalabel;
    QLineEdit *thetalineEdit;
    QLabel *thetalabel;
    QGroupBox *viewGB;
    GLDesign *widget;
    QPushButton *drawPB;
    QPushButton *drawmesh;
    QPushButton *refineB;
    QPushButton *pushButton;
    QMenuBar *menubar;
    QMenu *menuFile;
    QStatusBar *statusbar;
    QToolBar *toolBar;

    void setupUi(QMainWindow *DesignRail)
    {
        if (DesignRail->objectName().isEmpty())
            DesignRail->setObjectName(QString::fromUtf8("DesignRail"));
        DesignRail->resize(545, 349);
        actionSave_Rail = new QAction(DesignRail);
        actionSave_Rail->setObjectName(QString::fromUtf8("actionSave_Rail"));
        actionClose = new QAction(DesignRail);
        actionClose->setObjectName(QString::fromUtf8("actionClose"));
        centralwidget = new QWidget(DesignRail);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        parametersGB = new QGroupBox(centralwidget);
        parametersGB->setObjectName(QString::fromUtf8("parametersGB"));
        parametersGB->setGeometry(QRect(0, 0, 261, 251));
        parametersGB->setMaximumSize(QSize(261, 251));
        QFont font;
        font.setBold(true);
        font.setUnderline(false);
        font.setWeight(75);
        font.setStrikeOut(false);
        font.setKerning(true);
        parametersGB->setFont(font);
        HDlineEdit = new QLineEdit(parametersGB);
        HDlineEdit->setObjectName(QString::fromUtf8("HDlineEdit"));
        HDlineEdit->setGeometry(QRect(40, 26, 70, 20));
        HDlabel = new QLabel(parametersGB);
        HDlabel->setObjectName(QString::fromUtf8("HDlabel"));
        HDlabel->setGeometry(QRect(11, 28, 21, 16));
        FDlineEdit = new QLineEdit(parametersGB);
        FDlineEdit->setObjectName(QString::fromUtf8("FDlineEdit"));
        FDlineEdit->setGeometry(QRect(176, 26, 70, 20));
        FDlineEdit->setInputMask(QString::fromUtf8("00+(00/00); "));
        FDlabel = new QLabel(parametersGB);
        FDlabel->setObjectName(QString::fromUtf8("FDlabel"));
        FDlabel->setGeometry(QRect(147, 28, 21, 16));
        BDlineEdit = new QLineEdit(parametersGB);
        BDlineEdit->setObjectName(QString::fromUtf8("BDlineEdit"));
        BDlineEdit->setGeometry(QRect(39, 64, 70, 20));
        BDlineEdit->setInputMask(QString::fromUtf8("00+(00/00); "));
        BDlabel = new QLabel(parametersGB);
        BDlabel->setObjectName(QString::fromUtf8("BDlabel"));
        BDlabel->setGeometry(QRect(10, 66, 21, 16));
        HWlineEdit = new QLineEdit(parametersGB);
        HWlineEdit->setObjectName(QString::fromUtf8("HWlineEdit"));
        HWlineEdit->setGeometry(QRect(176, 64, 70, 20));
        HWlineEdit->setInputMask(QString::fromUtf8("00+(00/00); "));
        HWlabel = new QLabel(parametersGB);
        HWlabel->setObjectName(QString::fromUtf8("HWlabel"));
        HWlabel->setGeometry(QRect(147, 66, 21, 16));
        WlineEdit = new QLineEdit(parametersGB);
        WlineEdit->setObjectName(QString::fromUtf8("WlineEdit"));
        WlineEdit->setGeometry(QRect(39, 104, 70, 20));
        WlineEdit->setInputMask(QString::fromUtf8("00+(00/00); "));
        Wlabel = new QLabel(parametersGB);
        Wlabel->setObjectName(QString::fromUtf8("Wlabel"));
        Wlabel->setGeometry(QRect(10, 106, 21, 16));
        BWlineEdit = new QLineEdit(parametersGB);
        BWlineEdit->setObjectName(QString::fromUtf8("BWlineEdit"));
        BWlineEdit->setGeometry(QRect(176, 100, 70, 20));
        BWlineEdit->setInputMask(QString::fromUtf8("00+(00/00); "));
        BWlabel = new QLabel(parametersGB);
        BWlabel->setObjectName(QString::fromUtf8("BWlabel"));
        BWlabel->setGeometry(QRect(147, 105, 21, 16));
        R1HlineEdit = new QLineEdit(parametersGB);
        R1HlineEdit->setObjectName(QString::fromUtf8("R1HlineEdit"));
        R1HlineEdit->setGeometry(QRect(39, 142, 70, 20));
        R1HlineEdit->setInputMask(QString::fromUtf8("00+(00/00); "));
        R1Hlabel = new QLabel(parametersGB);
        R1Hlabel->setObjectName(QString::fromUtf8("R1Hlabel"));
        R1Hlabel->setGeometry(QRect(10, 144, 31, 16));
        R2HlineEdit = new QLineEdit(parametersGB);
        R2HlineEdit->setObjectName(QString::fromUtf8("R2HlineEdit"));
        R2HlineEdit->setGeometry(QRect(175, 142, 70, 20));
        R2HlineEdit->setInputMask(QString::fromUtf8("00+(00/00); "));
        R2Hlabel = new QLabel(parametersGB);
        R2Hlabel->setObjectName(QString::fromUtf8("R2Hlabel"));
        R2Hlabel->setGeometry(QRect(146, 144, 31, 16));
        RHBlineEdit = new QLineEdit(parametersGB);
        RHBlineEdit->setObjectName(QString::fromUtf8("RHBlineEdit"));
        RHBlineEdit->setGeometry(QRect(39, 182, 70, 20));
        RHBlineEdit->setInputMask(QString::fromUtf8("00+(00/00); "));
        RHBlabel = new QLabel(parametersGB);
        RHBlabel->setObjectName(QString::fromUtf8("RHBlabel"));
        RHBlabel->setGeometry(QRect(10, 184, 31, 16));
        alphalineEdit = new QLineEdit(parametersGB);
        alphalineEdit->setObjectName(QString::fromUtf8("alphalineEdit"));
        alphalineEdit->setGeometry(QRect(174, 182, 61, 20));
        alphalineEdit->setInputMask(QString::fromUtf8("00\313\23200'; "));
        alphalabel = new QLabel(parametersGB);
        alphalabel->setObjectName(QString::fromUtf8("alphalabel"));
        alphalabel->setGeometry(QRect(132, 184, 41, 16));
        betalineEdit = new QLineEdit(parametersGB);
        betalineEdit->setObjectName(QString::fromUtf8("betalineEdit"));
        betalineEdit->setGeometry(QRect(39, 220, 61, 20));
        betalineEdit->setInputMask(QString::fromUtf8("00\313\23200'; "));
        betalabel = new QLabel(parametersGB);
        betalabel->setObjectName(QString::fromUtf8("betalabel"));
        betalabel->setGeometry(QRect(4, 222, 31, 16));
        thetalineEdit = new QLineEdit(parametersGB);
        thetalineEdit->setObjectName(QString::fromUtf8("thetalineEdit"));
        thetalineEdit->setGeometry(QRect(175, 218, 61, 20));
        thetalineEdit->setInputMask(QString::fromUtf8("00\313\23200'; "));
        thetalabel = new QLabel(parametersGB);
        thetalabel->setObjectName(QString::fromUtf8("thetalabel"));
        thetalabel->setGeometry(QRect(133, 220, 41, 16));
        viewGB = new QGroupBox(centralwidget);
        viewGB->setObjectName(QString::fromUtf8("viewGB"));
        viewGB->setGeometry(QRect(281, 5, 251, 251));
        QFont font1;
        font1.setBold(true);
        font1.setWeight(75);
        viewGB->setFont(font1);
        widget = new GLDesign(viewGB);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setGeometry(QRect(9, 19, 231, 221));
        drawPB = new QPushButton(centralwidget);
        drawPB->setObjectName(QString::fromUtf8("drawPB"));
        drawPB->setGeometry(QRect(8, 255, 91, 41));
        drawmesh = new QPushButton(centralwidget);
        drawmesh->setObjectName(QString::fromUtf8("drawmesh"));
        drawmesh->setGeometry(QRect(102, 256, 91, 41));
        refineB = new QPushButton(centralwidget);
        refineB->setObjectName(QString::fromUtf8("refineB"));
        refineB->setGeometry(QRect(200, 256, 75, 41));
        pushButton = new QPushButton(centralwidget);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));
        pushButton->setGeometry(QRect(340, 256, 151, 41));
        DesignRail->setCentralWidget(centralwidget);
        menubar = new QMenuBar(DesignRail);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 545, 20));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        DesignRail->setMenuBar(menubar);
        statusbar = new QStatusBar(DesignRail);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        DesignRail->setStatusBar(statusbar);
        toolBar = new QToolBar(DesignRail);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
        DesignRail->addToolBar(Qt::TopToolBarArea, toolBar);

        menubar->addAction(menuFile->menuAction());
        menuFile->addAction(actionSave_Rail);
        menuFile->addAction(actionClose);

        retranslateUi(DesignRail);

        QMetaObject::connectSlotsByName(DesignRail);
    } // setupUi

    void retranslateUi(QMainWindow *DesignRail)
    {
        DesignRail->setWindowTitle(QApplication::translate("DesignRail", "DesignRail", 0, QApplication::UnicodeUTF8));
        actionSave_Rail->setText(QApplication::translate("DesignRail", "Save Rail", 0, QApplication::UnicodeUTF8));
        actionClose->setText(QApplication::translate("DesignRail", "Close", 0, QApplication::UnicodeUTF8));
        parametersGB->setTitle(QApplication::translate("DesignRail", "Geometry Parameters", 0, QApplication::UnicodeUTF8));
        HDlineEdit->setInputMask(QApplication::translate("DesignRail", "00+(00/00); ", 0, QApplication::UnicodeUTF8));
        HDlineEdit->setText(QApplication::translate("DesignRail", "1+(7/8)", 0, QApplication::UnicodeUTF8));
        HDlabel->setText(QApplication::translate("DesignRail", "HD:", 0, QApplication::UnicodeUTF8));
        FDlineEdit->setText(QApplication::translate("DesignRail", "3+(13/16)", 0, QApplication::UnicodeUTF8));
        FDlabel->setText(QApplication::translate("DesignRail", "FD:", 0, QApplication::UnicodeUTF8));
        BDlineEdit->setText(QApplication::translate("DesignRail", "1+(1/8)", 0, QApplication::UnicodeUTF8));
        BDlabel->setText(QApplication::translate("DesignRail", "BD:", 0, QApplication::UnicodeUTF8));
        HWlineEdit->setText(QApplication::translate("DesignRail", "2+(21/32)", 0, QApplication::UnicodeUTF8));
        HWlabel->setText(QApplication::translate("DesignRail", "HW:", 0, QApplication::UnicodeUTF8));
        WlineEdit->setText(QApplication::translate("DesignRail", "0+(5/8)", 0, QApplication::UnicodeUTF8));
        Wlabel->setText(QApplication::translate("DesignRail", "W:", 0, QApplication::UnicodeUTF8));
        BWlineEdit->setText(QApplication::translate("DesignRail", "5+(1/2)", 0, QApplication::UnicodeUTF8));
        BWlabel->setText(QApplication::translate("DesignRail", "BW:", 0, QApplication::UnicodeUTF8));
        R1HlineEdit->setText(QApplication::translate("DesignRail", "0+(1/2)", 0, QApplication::UnicodeUTF8));
        R1Hlabel->setText(QApplication::translate("DesignRail", "R1H:", 0, QApplication::UnicodeUTF8));
        R2HlineEdit->setText(QApplication::translate("DesignRail", "0+(1/16)", 0, QApplication::UnicodeUTF8));
        R2Hlabel->setText(QApplication::translate("DesignRail", "R2H:", 0, QApplication::UnicodeUTF8));
        RHBlineEdit->setText(QApplication::translate("DesignRail", "0+(3/4)", 0, QApplication::UnicodeUTF8));
        RHBlabel->setText(QApplication::translate("DesignRail", "RHB:", 0, QApplication::UnicodeUTF8));
        alphalineEdit->setText(QApplication::translate("DesignRail", "14\313\23202'", 0, QApplication::UnicodeUTF8));
        alphalabel->setText(QApplication::translate("DesignRail", "alpha:", 0, QApplication::UnicodeUTF8));
        betalineEdit->setText(QApplication::translate("DesignRail", "1\313\23226'", 0, QApplication::UnicodeUTF8));
        betalabel->setText(QApplication::translate("DesignRail", "beta:", 0, QApplication::UnicodeUTF8));
        thetalineEdit->setText(QApplication::translate("DesignRail", "14\313\23202'", 0, QApplication::UnicodeUTF8));
        thetalabel->setText(QApplication::translate("DesignRail", "theta:", 0, QApplication::UnicodeUTF8));
        viewGB->setTitle(QApplication::translate("DesignRail", "View", 0, QApplication::UnicodeUTF8));
        drawPB->setText(QApplication::translate("DesignRail", "draw \n"
"geometry", 0, QApplication::UnicodeUTF8));
        drawmesh->setText(QApplication::translate("DesignRail", "draw \n"
"mesh", 0, QApplication::UnicodeUTF8));
        refineB->setText(QApplication::translate("DesignRail", "refine", 0, QApplication::UnicodeUTF8));
        pushButton->setText(QApplication::translate("DesignRail", "Add to Workspace", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("DesignRail", "File", 0, QApplication::UnicodeUTF8));
        toolBar->setWindowTitle(QApplication::translate("DesignRail", "toolBar", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DesignRail: public Ui_DesignRail {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DESIGNRAIL_H
