/****************************************************************************
** Meta object code from reading C++ file 'graphics.h'
**
** Created: Mon 26. Aug 10:42:43 2019
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../graphics.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'graphics.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_graphics[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      23,   10,    9,    9, 0x05,

 // slots: signature, parameters, type, tag, flags
      53,    9,    9,    9, 0x08,
      72,    9,    9,    9, 0x08,
     107,    9,    9,    9, 0x08,
     132,    9,    9,    9, 0x08,
     157,    9,    9,    9, 0x08,
     193,    9,    9,    9, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_graphics[] = {
    "graphics\0\0row,col,type\0"
    "SelectionChanged(int,int,int)\0"
    "selectionChanged()\0"
    "NewSelected(QCPAbstractPlottable*)\0"
    "mousePress(QMouseEvent*)\0"
    "mouseWheel(QWheelEvent*)\0"
    "graphClicked(QCPAbstractPlottable*)\0"
    "graphMouseMove(QMouseEvent*)\0"
};

const QMetaObject graphics::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_graphics,
      qt_meta_data_graphics, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &graphics::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *graphics::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *graphics::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_graphics))
        return static_cast<void*>(const_cast< graphics*>(this));
    return QObject::qt_metacast(_clname);
}

int graphics::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: SelectionChanged((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])),(*reinterpret_cast< int(*)>(_a[3]))); break;
        case 1: selectionChanged(); break;
        case 2: NewSelected((*reinterpret_cast< QCPAbstractPlottable*(*)>(_a[1]))); break;
        case 3: mousePress((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        case 4: mouseWheel((*reinterpret_cast< QWheelEvent*(*)>(_a[1]))); break;
        case 5: graphClicked((*reinterpret_cast< QCPAbstractPlottable*(*)>(_a[1]))); break;
        case 6: graphMouseMove((*reinterpret_cast< QMouseEvent*(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 7;
    }
    return _id;
}

// SIGNAL 0
void graphics::SelectionChanged(int _t1, int _t2, int _t3)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)), const_cast<void*>(reinterpret_cast<const void*>(&_t3)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
