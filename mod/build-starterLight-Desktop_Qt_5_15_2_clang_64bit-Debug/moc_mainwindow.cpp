/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../starterLight/mainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_MainWindow_t {
    QByteArrayData data[13];
    char stringdata0[387];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MainWindow_t qt_meta_stringdata_MainWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "MainWindow"
QT_MOC_LITERAL(1, 11, 31), // "on_pushButton_suppBruit_clicked"
QT_MOC_LITERAL(2, 43, 0), // ""
QT_MOC_LITERAL(3, 44, 39), // "on_pushButton_suppBruitMailla..."
QT_MOC_LITERAL(4, 84, 32), // "on_pushButton_chargement_clicked"
QT_MOC_LITERAL(5, 117, 36), // "on_pushButton_identify_holes_..."
QT_MOC_LITERAL(6, 154, 31), // "on_pushButton_fix_holes_clicked"
QT_MOC_LITERAL(7, 186, 41), // "on_pushButton_identify_overla..."
QT_MOC_LITERAL(8, 228, 36), // "on_pushButton_fix_overlaping_..."
QT_MOC_LITERAL(9, 265, 31), // "on_pushButton_detectFis_clicked"
QT_MOC_LITERAL(10, 297, 28), // "on_pushButton_repFis_clicked"
QT_MOC_LITERAL(11, 326, 31), // "on_pushButton_voirBruit_clicked"
QT_MOC_LITERAL(12, 358, 28) // "on_pushButton_fixall_clicked"

    },
    "MainWindow\0on_pushButton_suppBruit_clicked\0"
    "\0on_pushButton_suppBruitMaillage_clicked\0"
    "on_pushButton_chargement_clicked\0"
    "on_pushButton_identify_holes_clicked\0"
    "on_pushButton_fix_holes_clicked\0"
    "on_pushButton_identify_overlaping_clicked\0"
    "on_pushButton_fix_overlaping_clicked\0"
    "on_pushButton_detectFis_clicked\0"
    "on_pushButton_repFis_clicked\0"
    "on_pushButton_voirBruit_clicked\0"
    "on_pushButton_fixall_clicked"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      11,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   69,    2, 0x08 /* Private */,
       3,    0,   70,    2, 0x08 /* Private */,
       4,    0,   71,    2, 0x08 /* Private */,
       5,    0,   72,    2, 0x08 /* Private */,
       6,    0,   73,    2, 0x08 /* Private */,
       7,    0,   74,    2, 0x08 /* Private */,
       8,    0,   75,    2, 0x08 /* Private */,
       9,    0,   76,    2, 0x08 /* Private */,
      10,    0,   77,    2, 0x08 /* Private */,
      11,    0,   78,    2, 0x08 /* Private */,
      12,    0,   79,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<MainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->on_pushButton_suppBruit_clicked(); break;
        case 1: _t->on_pushButton_suppBruitMaillage_clicked(); break;
        case 2: _t->on_pushButton_chargement_clicked(); break;
        case 3: _t->on_pushButton_identify_holes_clicked(); break;
        case 4: _t->on_pushButton_fix_holes_clicked(); break;
        case 5: _t->on_pushButton_identify_overlaping_clicked(); break;
        case 6: _t->on_pushButton_fix_overlaping_clicked(); break;
        case 7: _t->on_pushButton_detectFis_clicked(); break;
        case 8: _t->on_pushButton_repFis_clicked(); break;
        case 9: _t->on_pushButton_voirBruit_clicked(); break;
        case 10: _t->on_pushButton_fixall_clicked(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

QT_INIT_METAOBJECT const QMetaObject MainWindow::staticMetaObject = { {
    QMetaObject::SuperData::link<QMainWindow::staticMetaObject>(),
    qt_meta_stringdata_MainWindow.data,
    qt_meta_data_MainWindow,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 11)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 11;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 11)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 11;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
