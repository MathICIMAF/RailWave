/********************************************************************************
** Form generated from reading UI file 'datawindow.ui'
**
** Created: Thu 22. Aug 11:23:50 2019
**      by: Qt User Interface Compiler version 4.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DATAWINDOW_H
#define UI_DATAWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_datawindow
{
public:
    QWidget *centralwidget;
    QGroupBox *groupBox;
    QPushButton *cancelBtn;
    QLabel *infoLb;
    QProgressBar *progressBar;

    void setupUi(QMainWindow *datawindow)
    {
        if (datawindow->objectName().isEmpty())
            datawindow->setObjectName(QString::fromUtf8("datawindow"));
        datawindow->resize(390, 120);
        datawindow->setMinimumSize(QSize(390, 120));
        datawindow->setMaximumSize(QSize(390, 120));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/images.png"), QSize(), QIcon::Normal, QIcon::Off);
        datawindow->setWindowIcon(icon);
        centralwidget = new QWidget(datawindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        groupBox = new QGroupBox(centralwidget);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setGeometry(QRect(0, 0, 390, 120));
        cancelBtn = new QPushButton(groupBox);
        cancelBtn->setObjectName(QString::fromUtf8("cancelBtn"));
        cancelBtn->setGeometry(QRect(304, 85, 75, 23));
        infoLb = new QLabel(groupBox);
        infoLb->setObjectName(QString::fromUtf8("infoLb"));
        infoLb->setGeometry(QRect(25, 20, 351, 21));
        QFont font;
        font.setFamily(QString::fromUtf8("MS Sans Serif"));
        font.setPointSize(13);
        font.setBold(true);
        font.setItalic(true);
        font.setUnderline(false);
        font.setWeight(75);
        infoLb->setFont(font);
        progressBar = new QProgressBar(groupBox);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setGeometry(QRect(22, 50, 361, 23));
        progressBar->setValue(0);
        progressBar->setTextVisible(true);
        progressBar->setOrientation(Qt::Horizontal);
        progressBar->setInvertedAppearance(false);
        progressBar->setTextDirection(QProgressBar::TopToBottom);
        datawindow->setCentralWidget(centralwidget);

        retranslateUi(datawindow);

        QMetaObject::connectSlotsByName(datawindow);
    } // setupUi

    void retranslateUi(QMainWindow *datawindow)
    {
        datawindow->setWindowTitle(QApplication::translate("datawindow", "RailWave", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QString());
        cancelBtn->setText(QApplication::translate("datawindow", "Cancel", 0, QApplication::UnicodeUTF8));
        infoLb->setText(QApplication::translate("datawindow", "This process could take a few minutes...", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class datawindow: public Ui_datawindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DATAWINDOW_H
