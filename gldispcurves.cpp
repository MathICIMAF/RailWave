#include "gldispcurves.h"

GLDispCurves::GLDispCurves(QWidget *qw) :
    QGLWidget(qw)
{
}
void GLDispCurves::initializeGL(){
    glClearColor (1.0, 1.0, 1.0, 200.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
}

void GLDispCurves::paintGL(){
    glClear( GL_COLOR_BUFFER_BIT);
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
}
