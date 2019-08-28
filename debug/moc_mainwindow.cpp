/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Mon 26. Aug 10:42:30 2019
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      18,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x08,
      36,   11,   11,   11, 0x08,
      69,   11,   11,   11, 0x08,
      96,   11,   11,   11, 0x08,
     119,   11,   11,   11, 0x08,
     152,   11,   11,   11, 0x08,
     186,   11,   11,   11, 0x08,
     217,   11,   11,   11, 0x08,
     242,   11,   11,   11, 0x08,
     274,   11,   11,   11, 0x08,
     295,   11,   11,   11, 0x08,
     315,   11,   11,   11, 0x08,
     342,   11,   11,   11, 0x08,
     373,   11,   11,   11, 0x08,
     398,  396,   11,   11, 0x08,
     462,  432,   11,   11, 0x08,
     542,  519,   11,   11, 0x28,
     612,  594,   11,   11, 0x28,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0on_pushButton_clicked()\0"
    "on_actionCreate_Rail_triggered()\0"
    "on_actionAbout_triggered()\0"
    "on_animatePB_clicked()\0"
    "on_actionShow_Curves_triggered()\0"
    "on_spacingSp_valueChanged(double)\0"
    "on_secNumbSp_valueChanged(int)\0"
    "on_computeButt_clicked()\0"
    "on_actionClose_Rail_triggered()\0"
    "on_zoomOut_clicked()\0on_zoomIn_clicked()\0"
    "on_actionClose_triggered()\0"
    "on_actionOpen_Rail_triggered()\0"
    "Cancelled(datawindow*)\0,\0"
    "ShowDisplayMenu(bool,datawindow*)\0"
    "front,mesh,radius,rail,refine\0"
    "ShowRail(QList<QPointF>,Mesh*,double,RailDesigned*,bool)\0"
    "front,mesh,radius,rail\0"
    "ShowRail(QList<QPointF>,Mesh*,double,RailDesigned*)\0"
    "front,mesh,radius\0"
    "ShowRail(QList<QPointF>,Mesh*,double)\0"
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: on_pushButton_clicked(); break;
        case 1: on_actionCreate_Rail_triggered(); break;
        case 2: on_actionAbout_triggered(); break;
        case 3: on_animatePB_clicked(); break;
        case 4: on_actionShow_Curves_triggered(); break;
        case 5: on_spacingSp_valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 6: on_secNumbSp_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: on_computeButt_clicked(); break;
        case 8: on_actionClose_Rail_triggered(); break;
        case 9: on_zoomOut_clicked(); break;
        case 10: on_zoomIn_clicked(); break;
        case 11: on_actionClose_triggered(); break;
        case 12: on_actionOpen_Rail_triggered(); break;
        case 13: Cancelled((*reinterpret_cast< datawindow*(*)>(_a[1]))); break;
        case 14: ShowDisplayMenu((*reinterpret_cast< bool(*)>(_a[1])),(*reinterpret_cast< datawindow*(*)>(_a[2]))); break;
        case 15: ShowRail((*reinterpret_cast< QList<QPointF>(*)>(_a[1])),(*reinterpret_cast< Mesh*(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])),(*reinterpret_cast< RailDesigned*(*)>(_a[4])),(*reinterpret_cast< bool(*)>(_a[5]))); break;
        case 16: ShowRail((*reinterpret_cast< QList<QPointF>(*)>(_a[1])),(*reinterpret_cast< Mesh*(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])),(*reinterpret_cast< RailDesigned*(*)>(_a[4]))); break;
        case 17: ShowRail((*reinterpret_cast< QList<QPointF>(*)>(_a[1])),(*reinterpret_cast< Mesh*(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3]))); break;
        default: ;
        }
        _id -= 18;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
