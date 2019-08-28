#include "glsection.h"

using namespace std;

GLSection::GLSection(QWidget*qw) :
    QGLWidget(qw)
{
    border = false;
    this->dmesh = false;
    this->bcontrol = false;
}
GLSection::GLSection(){
    border = false;
    this->dmesh = false;
}

void GLSection::initializeGL(){
    glClearColor (1.0, 1.0, 1.0, 200.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
}

void GLSection::paintGL(){
    glClear( GL_COLOR_BUFFER_BIT);
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho(-radius, radius, -radius, radius, 10, -10);
    glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();    
    glScalef( scale, scale, scale );

    if(dmesh){
        glPushMatrix();
            displayMesh();
        glPopMatrix();
    }

    if(border){
        glLineWidth(2.0);
        glColor3f(0,0,1);

        glBegin(GL_LINE_LOOP);
        for(int i = 0; i < frontier.count(); i++)
            glVertex2d(frontier[i].x(),frontier[i].y());
        glEnd();
        glLineWidth(1.0);
        glColor3f(1,1,1);
    }
    if(bcontrol){
        glLineWidth(2.0);
        glColor3f(1,0,0);

        glBegin(GL_LINE_LOOP);
        for(int i = 0; i < control.count(); i++)
            glVertex2d(control[i]._x,control[i]._y);
        glEnd();
        glLineWidth(1.0);
        glColor3f(1,1,1);
    }
}

void GLSection::drawBorder(QList<QPointF> pts){
    this->border = true;
    frontier = pts;
    LocateCenter();
    Center();
    updateGL();
}

void GLSection::drawControl(QList<PointS> pts){
    control = pts;
    this->bcontrol = true;
    updateGL();
}

void GLSection::Center(){
    scale = 1;
    obj_pos[0] = center[0];
    obj_pos[1] = center[1];
    obj_pos[2] = center[2];
}

void GLSection::LocateCenter(){
    float minX = std::numeric_limits<float>::max();
    float minY = minX;
    float maxX = -std::numeric_limits<float>::min();
    float maxY = maxX;

    for (int i = 0; i < frontier.count(); i++)
    {
        if(minX > frontier[i].x())
            minX = frontier[i].x();
        if(minY > frontier[i].y())
            minY = frontier[i].y();
        if(maxX < frontier[i].x())
            maxX = frontier[i].x();
        if(maxY < frontier[i].y())
            maxY = frontier[i].y();
    }
    center[0] = (maxX + minX)/2;
    center[1] = -(maxY + minY)/2;
    center[2] = 0;
    this->radius = max((maxX - minX)/2, (maxY - minY)/2)+0.01;
}

void GLSection::resizeGL(int w, int h){
    glViewport(0, 0, w, h);
}

void GLSection::drawMesh(Mesh *mesh){
    //if(!this->mesh.isNull())
        //delete (this->mesh);
    this->mesh = mesh;
    dmesh = true;
    //GetRadius();
    //Center();
    updateGL();
}

void GLSection::displayMesh(){
    glLineWidth(2.0);    
    glColor3f(0,0,0);
    for (FaceIter f = mesh->faces_begin();f != mesh->faces_end(); f++)
    {
        Pvert *v0 = Pvert::cast((*f)->vertex(0));
        Pvert *v1 = Pvert::cast((*f)->vertex(1));
        Pvert *v2 = Pvert::cast((*f)->vertex(2));        
        drawTriangle(v0->a.g[0], v0->a.g[1],
                     v1->a.g[0], v1->a.g[1],
                     v2->a.g[0], v2->a.g[1]);
    }
    glLineWidth(1.0);
    glColor3f(1,1,1);
}

void GLSection::drawTriangle(float x0, float y0, float x1, float y1, float x2, float y2){
    glBegin(GL_LINE_LOOP);
        glVertex2d(x0, y0);
        glVertex2d(x1, y1);
        glVertex2d(x2, y2);        
    glEnd();
}

void GLSection::GetRadius(){
    float minX = std::numeric_limits<float>::max();
    float minY = minX;
    float maxX = -std::numeric_limits<float>::min();
    float maxY = maxX;

    for (VertexIter v = mesh->verts_begin(); v != mesh->verts_end(); v++)
    {
        Pvert *v0 = Pvert::cast((*v));
        if(minX > v0->a.g.ax)
            minX = v0->a.g.ax;
        if(minY > v0->a.g.ay)
            minY = v0->a.g.ay;
        if(maxX < v0->a.g.ax)
            maxX = v0->a.g.ax;
        if(maxY < v0->a.g.ay)
            maxY = v0->a.g.ay;
    }
    center[0] = (maxX + minX)/2;
    center[1] = -(maxY + minY)/2;
    center[2] = 0;
    this->radius = max((maxX - minX)/2, (maxY - minY)/2);
}

void GLSection::GetTrianglesBoundaryRelation(QList<double> &x, QList<double> &y, QList<QList<int> > &triangles, QList<QList<int> > &boundaryIndices, QList<QList<GVector2> > &boundary){    
    QList<QList<Vertex*> > boundVertex = GetBoundary(mesh);    
    x.begin();y.begin();
    int id = 0;

    QList<QList<int> > index;
    for(int i = 0; i < boundVertex.count();i++)
    {
        QList<int> ith;
        for(int j = 0; j < boundVertex[i].count(); j++)
            ith.append(-1);
        index.append(ith);
    }


    for(VertexIter vIter = mesh->verts_begin(); vIter != mesh->verts_end(); vIter ++)
    {
        Vertex *vert = (*vIter);
        vert->setId(id++);
        Pvert *v = Pvert::cast(vert);
        x.push_back(v->a.g.ax);
        y.push_back(v->a.g.ay);

        QPoint i = Contains(boundVertex, vert);
        if(i.x() != -1)
            index[i.x()][i.y()] = vert->getId();
    }


    for(int i = 0; i<boundVertex.count();i++)
    {
        boundaryIndices.append(QList<int>());
        for(int j = 0; j<boundVertex[i].count();j++)
            boundaryIndices[i].append(index[i][j]);
    }

    for(FaceIter fIter = mesh->faces_begin(); fIter != mesh->faces_end(); fIter ++)
    {
            Face *f = *(fIter);
            f->setSaved(-1);
    }

    for(FaceIter fIter = mesh->faces_begin(); fIter != mesh->faces_end(); fIter ++)
    {
            Face *f = *(fIter);
            if(f->getSaved() == 1)continue;

            QList<int> ti;
            ti.append(Pvert::cast(f->vertex(1))->getId());
            ti.append(Pvert::cast(f->vertex(2))->getId());
            ti.append(Pvert::cast(f->vertex(3))->getId());
            triangles.push_back(ti);
            f->setSaved(1);

            Hedge *hedge = f->hedge(0);
            if(hedge->mate() != NULL)
            {
                Face *fMate = hedge->mate()->face();
                if(fMate!=NULL)
                {
                    if(fMate->getSaved() == 1)
                        continue;
                    QList<int> tiM;
                    tiM.append(Pvert::cast(fMate->vertex(1))->getId());
                    tiM.append(Pvert::cast(fMate->vertex(2))->getId());
                    tiM.append(Pvert::cast(fMate->vertex(3))->getId());
                    triangles.push_back(tiM);
                    fMate->setSaved(1);
                }
            }
        }

    for(int i = 0; i<boundaryIndices.count(); i++)
    {
        boundary.append(QList<GVector2>());
        for(int j = 0; j < boundaryIndices[i].count(); j++)
            boundary[i].push_back(GVector2(x.at(boundaryIndices[i][j]), y.at(boundaryIndices[i][j])));
    }
}

void GLSection::CleanControl(){
    if(!dmesh)
        return;
    radius = 0;
    delete mesh;
    dmesh = border = bcontrol = false;
    frontier.clear();
    control.clear();
    updateGL();
}
