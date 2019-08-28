/****************************************************************************
** Meta object code from reading C++ file 'datawindow.h'
**
** Created: Mon 26. Aug 10:42:40 2019
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../datawindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'datawindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_datawindow[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       4,       // signalCount

 // signals: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x05,
      32,   11,   11,   11, 0x05,
      57,   55,   11,   11, 0x05,
      97,   88,   11,   11, 0x05,

 // slots: signature, parameters, type, tag, flags
     128,   11,   11,   11, 0x08,
     137,   11,   11,   11, 0x08,
     160,   11,   11,   11, 0x08,
     189,  185,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_datawindow[] = {
    "datawindow\0\0SetProgressBar(int)\0"
    "Cancelled(datawindow*)\0,\0"
    "CurvesToShow(bool,datawindow*)\0val,max,\0"
    "SetHidden(int,int,datawindow*)\0finish()\0"
    "on_cancelBtn_clicked()\0closeEvent(QCloseEvent*)\0"
    "pos\0setProgressBarPos(int)\0"
};

const QMetaObject datawindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_datawindow,
      qt_meta_data_datawindow, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &datawindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *datawindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *datawindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_datawindow))
        return static_cast<void*>(const_cast< datawindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int datawindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: SetProgressBar((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: Cancelled((*reinterpret_cast< datawindow*(*)>(_a[1]))); break;
        case 2: CurvesToShow((*reinterpret_cast< bool(*)>(_a[1])),(*reinterpret_cast< datawindow*(*)>(_a[2]))); break;
        case 3: SetHidden((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])),(*reinterpret_cast< datawindow*(*)>(_a[3]))); break;
        case 4: finish(); break;
        case 5: on_cancelBtn_clicked(); break;
        case 6: closeEvent((*reinterpret_cast< QCloseEvent*(*)>(_a[1]))); break;
        case 7: setProgressBarPos((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 8;
    }
    return _id;
}

// SIGNAL 0
void datawindow::SetProgressBar(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void datawindow::Cancelled(datawindow * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void datawindow::CurvesToShow(bool _t1, datawindow * _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void datawindow::SetHidden(int _t1, int _t2, datawindow * _t3)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)), const_cast<void*>(reinterpret_cast<const void*>(&_t3)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}
QT_END_MOC_NAMESPACE
