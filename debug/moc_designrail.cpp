/****************************************************************************
** Meta object code from reading C++ file 'designrail.h'
**
** Created: Mon 26. Aug 10:42:46 2019
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../designrail.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'designrail.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_DesignRail[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      38,   12,   11,   11, 0x05,

 // slots: signature, parameters, type, tag, flags
      98,   11,   11,   11, 0x08,
     122,   11,   11,   11, 0x08,
     149,   11,   11,   11, 0x08,
     180,   11,   11,   11, 0x08,
     201,   11,   11,   11, 0x08,
     223,   11,   11,   11, 0x08,
     248,  243,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_DesignRail[] = {
    "DesignRail\0\0front,mesh,radius,,refine\0"
    "AddThisRail(QList<QPointF>,Mesh*,double,RailDesigned*,bool)\0"
    "on_pushButton_clicked()\0"
    "on_actionClose_triggered()\0"
    "on_actionSave_Rail_triggered()\0"
    "on_refineB_clicked()\0on_drawmesh_clicked()\0"
    "on_drawPB_clicked()\0text\0"
    "on_HDlineEdit_textChanged(QString)\0"
};

const QMetaObject DesignRail::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_DesignRail,
      qt_meta_data_DesignRail, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &DesignRail::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *DesignRail::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *DesignRail::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_DesignRail))
        return static_cast<void*>(const_cast< DesignRail*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int DesignRail::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: AddThisRail((*reinterpret_cast< QList<QPointF>(*)>(_a[1])),(*reinterpret_cast< Mesh*(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])),(*reinterpret_cast< RailDesigned*(*)>(_a[4])),(*reinterpret_cast< bool(*)>(_a[5]))); break;
        case 1: on_pushButton_clicked(); break;
        case 2: on_actionClose_triggered(); break;
        case 3: on_actionSave_Rail_triggered(); break;
        case 4: on_refineB_clicked(); break;
        case 5: on_drawmesh_clicked(); break;
        case 6: on_drawPB_clicked(); break;
        case 7: on_HDlineEdit_textChanged((*reinterpret_cast< QString(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 8;
    }
    return _id;
}

// SIGNAL 0
void DesignRail::AddThisRail(QList<QPointF> _t1, Mesh * _t2, double _t3, RailDesigned * _t4, bool _t5)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)), const_cast<void*>(reinterpret_cast<const void*>(&_t3)), const_cast<void*>(reinterpret_cast<const void*>(&_t4)), const_cast<void*>(reinterpret_cast<const void*>(&_t5)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
