/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.15.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <meshviewerwidget.h>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout;
    QWidget *widget_2;
    QVBoxLayout *verticalLayout;
    QPushButton *pushButton_chargement;
    QSpacerItem *verticalSpacer;
    QPushButton *pushButton_voirBruit;
    QPushButton *pushButton_suppBruit;
    QPushButton *pushButton_suppBruitMaillage;
    QSpacerItem *verticalSpacer1;
    QPushButton *pushButton_identify_overlaping;
    QPushButton *pushButton_fix_overlaping;
    QSpacerItem *verticalSpacer2;
    QPushButton *pushButton_detectFis;
    QPushButton *pushButton_repFis;
    QSpacerItem *verticalSpacer3;
    QPushButton *pushButton_identify_holes;
    QPushButton *pushButton_fix_holes;
    QSpacerItem *verticalSpacer4;
    QPushButton *pushButton_fixall;
    QSpacerItem *verticalSpacer5;
    MeshViewerWidget *displayWidget;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1200, 800);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        widget_2 = new QWidget(centralWidget);
        widget_2->setObjectName(QString::fromUtf8("widget_2"));
        QSizePolicy sizePolicy(QSizePolicy::Maximum, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(150);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(widget_2->sizePolicy().hasHeightForWidth());
        widget_2->setSizePolicy(sizePolicy);
        widget_2->setMinimumSize(QSize(150, 0));
        verticalLayout = new QVBoxLayout(widget_2);
        verticalLayout->setSpacing(4);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(3, 3, 3, 3);
        pushButton_chargement = new QPushButton(widget_2);
        pushButton_chargement->setObjectName(QString::fromUtf8("pushButton_chargement"));
        pushButton_chargement->setMinimumSize(QSize(200, 0));

        verticalLayout->addWidget(pushButton_chargement);

        verticalSpacer = new QSpacerItem(20, 80, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        pushButton_voirBruit = new QPushButton(widget_2);
        pushButton_voirBruit->setObjectName(QString::fromUtf8("pushButton_voirBruit"));

        verticalLayout->addWidget(pushButton_voirBruit);

        pushButton_suppBruit = new QPushButton(widget_2);
        pushButton_suppBruit->setObjectName(QString::fromUtf8("pushButton_suppBruit"));

        verticalLayout->addWidget(pushButton_suppBruit);

        pushButton_suppBruitMaillage = new QPushButton(widget_2);
        pushButton_suppBruitMaillage->setObjectName(QString::fromUtf8("pushButton_suppBruitMaillage"));

        verticalLayout->addWidget(pushButton_suppBruitMaillage);

        verticalSpacer1 = new QSpacerItem(20, 80, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer1);

        pushButton_identify_overlaping = new QPushButton(widget_2);
        pushButton_identify_overlaping->setObjectName(QString::fromUtf8("pushButton_identify_overlaping"));

        verticalLayout->addWidget(pushButton_identify_overlaping);

        pushButton_fix_overlaping = new QPushButton(widget_2);
        pushButton_fix_overlaping->setObjectName(QString::fromUtf8("pushButton_fix_overlaping"));

        verticalLayout->addWidget(pushButton_fix_overlaping);

        verticalSpacer2 = new QSpacerItem(20, 80, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer2);

        pushButton_detectFis = new QPushButton(widget_2);
        pushButton_detectFis->setObjectName(QString::fromUtf8("pushButton_detectFis"));

        verticalLayout->addWidget(pushButton_detectFis);

        pushButton_repFis = new QPushButton(widget_2);
        pushButton_repFis->setObjectName(QString::fromUtf8("pushButton_repFis"));

        verticalLayout->addWidget(pushButton_repFis);

        verticalSpacer3 = new QSpacerItem(20, 80, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer3);

        pushButton_identify_holes = new QPushButton(widget_2);
        pushButton_identify_holes->setObjectName(QString::fromUtf8("pushButton_identify_holes"));

        verticalLayout->addWidget(pushButton_identify_holes);

        pushButton_fix_holes = new QPushButton(widget_2);
        pushButton_fix_holes->setObjectName(QString::fromUtf8("pushButton_fix_holes"));

        verticalLayout->addWidget(pushButton_fix_holes);

        verticalSpacer4 = new QSpacerItem(20, 80, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer4);

        pushButton_fixall = new QPushButton(widget_2);
        pushButton_fixall->setObjectName(QString::fromUtf8("pushButton_fixall"));

        verticalLayout->addWidget(pushButton_fixall);

        verticalSpacer5 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer5);


        horizontalLayout->addWidget(widget_2);

        displayWidget = new MeshViewerWidget(centralWidget);
        displayWidget->setObjectName(QString::fromUtf8("displayWidget"));

        horizontalLayout->addWidget(displayWidget);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 632, 21));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));
        pushButton_chargement->setText(QCoreApplication::translate("MainWindow", "Charger OBJ", nullptr));
        pushButton_voirBruit->setText(QCoreApplication::translate("MainWindow", "Voir le bruit", nullptr));
        pushButton_suppBruit->setText(QCoreApplication::translate("MainWindow", "Supprimer le bruit", nullptr));
        pushButton_suppBruitMaillage->setText(QCoreApplication::translate("MainWindow", "Supprimer le bruit du maillage", nullptr));
        pushButton_identify_overlaping->setText(QCoreApplication::translate("MainWindow", "Identifier Overlaping", nullptr));
        pushButton_fix_overlaping->setText(QCoreApplication::translate("MainWindow", "R\303\251parer Overlaping", nullptr));
        pushButton_detectFis->setText(QCoreApplication::translate("MainWindow", "D\303\251tecter fissure", nullptr));
        pushButton_repFis->setText(QCoreApplication::translate("MainWindow", "R\303\251parer fissure", nullptr));
        pushButton_identify_holes->setText(QCoreApplication::translate("MainWindow", "Identifier les trous", nullptr));
        pushButton_fix_holes->setText(QCoreApplication::translate("MainWindow", "R\303\251parer les trous", nullptr));
        pushButton_fixall->setText(QCoreApplication::translate("MainWindow", "R\303\251paration compl\303\250te", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
