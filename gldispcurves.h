#ifndef GLDISPCURVES_H
#define GLDISPCURVES_H

#include <QGLWidget>

class GLDispCurves : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLDispCurves(QWidget *qw);
    GLDispCurves();
protected:
    void paintGL();
    virtual void initializeGL();

signals:

public slots:

};

#endif // GLDISPCURVES_H
