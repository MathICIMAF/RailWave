#include "gldesign.h"
//#include "Plane.h"

GLDesign::GLDesign()
{
    bmesh = false;
}

GLDesign::GLDesign(QWidget *qw) :
    QGLWidget(qw)
{
    bmesh = false;
}

void GLDesign::initializeGL(){
    glClearColor (1.0, 1.0, 1.0, 200.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    bcontrol = false;
    level = 0;
    scale = 1;
}

void GLDesign::paintGL(){
    glClear( GL_COLOR_BUFFER_BIT);
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();    


    glOrtho(-radius, radius, -radius, radius, 10, -10);
    glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();    
    glScalef( scale, scale, scale );
    if(bcontrol){        
        glLineWidth(2.0);
        glColor3f(1,0,0);

        glBegin(GL_LINE_LOOP);//GL_LINE_LOOP
            for(int i = 1; i < polyg.size(); i+=2)
                glVertex2d(polyg[i].x(),polyg[i].y());            

        glEnd();
        glLineWidth(1.0);
        glColor3f(1,1,1);
    }
    if(bmesh){
        for (FaceIter f = mesh->faces_begin(); f != mesh->faces_end(); f++) {
            Pvert *v0 = Pvert::cast((*f)->vertex(0));
            Pvert *v1 = Pvert::cast((*f)->vertex(1));
            Pvert *v2 = Pvert::cast((*f)->vertex(2));
            draw_tri(v0->a.g[0], v0->a.g[1], v0->a.g[2],
                      v1->a.g[0], v1->a.g[1], v1->a.g[2],
                      v2->a.g[0], v2->a.g[1], v2->a.g[2]);
          }
    }
}


void GLDesign::resizeGL(int w, int h){
    glViewport(0, 0, w, h);
}

void GLDesign::drawpolygon(RailPolygon *pol){
    bcontrol = true;
    this->polyg = pol->Polygon();
    this->CInitial = pol->C();
    this->radius = pol->GetRadius();
    ratio = radius/10;
    center[0] = pol->center.x();
    center[1] = pol->center.y();
    Center();
    updateGL();
}

void GLDesign::drawmesh(Mesh*mesh,RailDesigned* rail){
   this->rail = rail;      
   this->mesh = mesh;   
   bmesh = true;
   updateGL();
}

void GLDesign::draw_tri(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2){
    glColor3f (0.0, 0.0, 0.0);
      glBegin(GL_LINE_LOOP);
      glVertex3d(x0, y0, z0);
      glVertex3d(x1, y1, z1);
      glVertex3d(x2, y2, z2);
      glEnd();
}

void GLDesign::refine(){
    rail->target_level = level++;
    mesh->adapt_refine(0);
    updateGL();
}

void GLDesign::ZoomIn(){
    scale *=1.2;
    updateGL();
}

void GLDesign::ZoomOut(){
    scale /=1.2;
    updateGL();
}

void GLDesign::wheelEvent(QWheelEvent *wheelEvent){
    (wheelEvent->delta() > 0)?ZoomIn():ZoomOut();
}

void GLDesign::Center(){
    scale = 1;
    obj_pos[0] = center[0];
    obj_pos[1] = center[1];
    obj_pos[2] = center[2];
    updateGL();
}

void GLDesign::keyPressEvent(QKeyEvent *keyEvent){
    switch(keyEvent->key())
    {
        case Qt::Key_Up:
            MoveUp();
            break;
        case Qt::Key_Down:
            MoveDown();
            break;
        case Qt::Key_Right:
            MoveRight();
            break;
        case Qt::Key_Left:
            MoveLeft();
            break;
        case Qt::Key_Plus:
            ZoomIn();
            break;
        case Qt::Key_Minus:
            ZoomOut();
            break;
    }
}

void GLDesign::MoveDown(){
    obj_pos[1] -= ratio;
    updateGL();
}

void GLDesign::MoveLeft(){
    obj_pos[0] -= ratio;
    updateGL();
}

void GLDesign::MoveRight(){
    obj_pos[0] += ratio;
    updateGL();
}

void GLDesign::MoveUp(){
    obj_pos[1] += ratio;
    updateGL();
}
