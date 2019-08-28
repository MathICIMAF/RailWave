#-------------------------------------------------
#
# Project created by QtCreator 2016-12-16T15:13:19
#
#-------------------------------------------------

QT       += core gui

TARGET = WRail
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    gvector2.cpp \
    Railneeds.cpp \
    glsection.cpp \
    Layer.cpp \
    glrail3d.cpp \
    gldispcurves.cpp \
    A48/triquad.cpp \
    A48/stellar.cpp \
    A48/mesh.cpp \
    A48/heap.cpp \
    A48/front.cpp \
    A48/face.cpp \
    A48/edge.cpp \
    A48/adapt.cpp \
    gvector3.cpp \
    Surfaces.cpp \
    Needs.cpp \
    Matrix.cpp \
             eigen/zvout.f \
             eigen/zunm2r.f \
             eigen/ztrsyl.f \
             eigen/ztrsv.f \
             eigen/ztrsm.f \
             eigen/ztrsen.f \
             eigen/ztrmm.f \
             eigen/ztrexc.f \
             eigen/ztrevc.f \
             eigen/zswap.f \
             eigen/zstatn.f \
             eigen/zsortc.f \
             eigen/zscal.f \
             eigen/zrot.f \
             eigen/zPackcg.f \
             eigen/zngets.f \
             eigen/zneupd.f \
             eigen/zneigh.f \
             eigen/znaupd.f \
             eigen/znaup2.f \
             eigen/znapps.f \
             eigen/znaitr.f \
             eigen/zmout.f \
             eigen/zlatrs.f \
             eigen/zlaswp.f \
             eigen/zlassq.f \
             eigen/zlaset.f \
             eigen/zlascl.f \
             eigen/zlartg.f \
             eigen/zlarnv.f \
             eigen/zlarfg.f \
             eigen/zlarf.f \
             eigen/zlanhs.f \
             eigen/zlange.f \
             eigen/zlahqr.f \
             eigen/zladiv.f \
             eigen/zlacpy.f \
             eigen/zlacon.f \
             eigen/zlacn2.f \
             eigen/zgetv0.f \
             eigen/zgetrs.f \
             eigen/zgetrf.f \
             eigen/zgetf2.f \
             eigen/zgeru.f \
             eigen/zgerc.f \
             eigen/zgeqr2.f \
             eigen/zgemv.f \
             eigen/zgemm.f \
             eigen/zgecon.f \
             eigen/zdscal.f \
             eigen/zdrscl.f \
             eigen/zdotu.f \
             eigen/zdotc.f \
             eigen/zcopy.f \
             eigen/zaxpy.f \
             eigen/xerbla.f \
             eigen/svout.f \
             eigen/sscal.f \
             eigen/smout.f \
             eigen/slaruv.f \
             eigen/slapy3.f \
             eigen/slapy2.f \
             eigen/slamch.f \
             eigen/sladiv.f \
             eigen/slabad.f \
             eigen/second.f \
             eigen/scsum1.f \
             eigen/scnrm2.f \
             eigen/scasum.f \
             eigen/lsame.f \
             eigen/izmax1.f \
             eigen/izamax.f \
             eigen/ivout.f \
             eigen/isamax.f \
             eigen/ilazlr.f \
             eigen/ilazlc.f \
             eigen/ilaenv.f \
             eigen/idamax.f \
             eigen/icmax1.f \
             eigen/icamax.f \
             eigen/eigensolverq.f \
             eigen/eigensolvermxmaxb.f \
             eigen/eigensolvermxm.f \
             eigen/dzsum1.f \
             eigen/dznrm2.f \
             eigen/dzasum.f \
             eigen/dvout.f \
             eigen/dswap.f \
             eigen/dstev.f \
             eigen/dsterf.f \
             eigen/dsteqr.f \
             eigen/dscal.f \
             eigen/dlassq.f \
             eigen/dlasrt.f \
             eigen/dlasr.f \
             eigen/dlaset.f \
             eigen/dlascl.f \
             eigen/dlaruv.f \
             eigen/dlartg.f \
             eigen/dlapy3.f \
             eigen/dlapy2.f \
             eigen/dlanst.f \
             eigen/dlamch.f \
             eigen/dlaisnan.f \
             eigen/dlaev2.f \
             eigen/dlae2.f \
             eigen/dladiv.f \
             eigen/dlabad.f \
             eigen/disnan.f \
             eigen/ctrsm.f \
             eigen/cpotrs.f \
             eigen/zpotrs.f \
             eigen/zpotrf.f \
             eigen/zpotf2.f \
             eigen/zlacgv.f \
             eigen/zherk.f \
             eigen/iparmq.f \
             eigen/ieeeck.f \
             eigen/zunmqr.f \
             eigen/zunmlq.f \
             eigen/zunml2.f \
             eigen/zunmbr.f \
             eigen/zungqr.f \
             eigen/zunglq.f \
             eigen/zungl2.f \
             eigen/zungbr.f \
             eigen/zung2r.f \
             eigen/ztrmv.f \
             eigen/zlasr.f \
             eigen/zlarft.f \
             eigen/zlarfb.f \
             eigen/zlabrd.f \
             eigen/zgesvd.f \
             eigen/zgeqrf.f \
             eigen/zgelqf.f \
             eigen/zgelq2.f \
             eigen/zgebrd.f \
             eigen/zgebd2.f \
             eigen/zdrot.f \
             eigen/zbdsqr.f \
             eigen/dlasv2.f \
             eigen/dlasq6.f \
             eigen/dlasq5.f \
             eigen/dlasq4.f \
             eigen/dlasq3.f \
             eigen/dlasq2.f \
             eigen/dlasq1.f \
             eigen/dlas2.f \
             eigen/dcopy.f \
    qcustomplot.cpp \
    Utils.cpp \
    WavesUtils.cpp \
    datawindow.cpp \
    graphics.cpp \
    dialog.cpp \
    designrail.cpp \
    gldesign.cpp

HEADERS  += mainwindow.h \
    RailNeeds.h \
    binheap.h \
    gvector2.h \
    glsection.h \
    Layer.h \
    Constants.h \
    glrail3d.h \
    gldispcurves.h \
    A48/vertex.h \
    A48/Mesh.h \
    A48/heap.h \
    A48/face.h \
    A48/edge.h \
    A48/a48.h \
    gvector3.h \
    Surfaces.h \
    Needs.h \
    Matrix.h \
    eigen/stat.h \
    eigen/linker.h \
    eigen/debug.h \
    qcustomplot.h \
    Utils.h \
    WavesUtils.h \
    datawindow.h \
    graphics.h \
    dialog.h \
    designrail.h \
    gldesign.h


FORMS    += mainwindow.ui \
    datawindow.ui \
    dialog.ui \
    designrail.ui

QT        += opengl

LIBS      += -lgfortran

RESOURCES += \
    rsources.qrc
