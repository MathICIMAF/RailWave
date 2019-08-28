#include "glrail3d.h"

GLRail3d::GLRail3d(QWidget *qw) :
    QGLWidget(qw)
{
    azimuth_ = 28.8;
    elevation_ = -18.6;
    dmesh = animate = false;
    timer = new QTimer(this);
    row = col = -1;
    connect(timer,SIGNAL(timeout()),this,SLOT(timeElapsed()));
}

void GLRail3d::initializeGL(){
    setFocusPolicy(Qt::StrongFocus);
    glClearColor (1.0, 1.0, 1.0, 1.0);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
    glDepthRange(0.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
}



void GLRail3d::paintGL(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho(-radius, radius, -radius, radius, 10, -10);
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] );
    glRotatef(elevation_,1,0,0);
    glRotatef(azimuth_,0,1,0);

    glScalef( scale, scale, scale );

    if(dmesh)
     {
        glPushMatrix();
            DisplayMesh();
        glPopMatrix();
     }
}

QList<QList<QList<GVector3*> > > GLRail3d::InitDisplacements(int boundaries, int count, int sections){
    QList<QList<QList<GVector3*> > >result;
    //srand ( time(NULL) );

    for(int k = 0; k < boundaries; k++)
    {
        result.append(QList<QList<GVector3*> >());
        for(int i = 0; i < sections; i++)
        {
            QList<GVector3*> iesimo;
            for(int j=0;j<count;j++)
                iesimo.append(new GVector3(0,0,0));
            result[k].append(iesimo);
        }
    }
    return result;
}

void GLRail3d::GetTrianglesBoundaryRelation(){
    triangles.begin();
    boundaryIndices.begin();
    boundary.begin();
    QList<QList<Vertex*> > boundVertex = GetBoundary(initial);
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


    for(VertexIter vIter = initial->verts_begin(); vIter != initial->verts_end(); vIter ++)
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

    for(FaceIter fIter = initial->faces_begin(); fIter != initial->faces_end(); fIter ++)
    {
            Face *f = *(fIter);
            f->setSaved(-1);
    }

    for(FaceIter fIter = initial->faces_begin(); fIter != initial->faces_end(); fIter ++)
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

void GLRail3d::ShowShape(Mesh *sec, int sections, double spacing){    
    initial = sec;
    dmesh = true;    
    this->spacing = spacing;
    this->sections = sections;    
    GetTrianglesBoundaryRelation();
    disp = InitDisplacements(boundary.count(),boundary[0].count(),sections);    
    mesh = new MeshCuad(boundary,boundaryIndices,spacing,sections,disp);
    mesh->NoDisplacements();
    LocateCenter();
    Center();
}

void GLRail3d::Animate(QList<QList<AvalAvect> > disp,QList<QList<AvalAvect> > pfDisplacements, QList<QList<AvalAvect> > plDisplacements, QList<QList<AvalAvect> > ptDisplacements){
    displacements = disp;
    fDisplacements = pfDisplacements;
    lDisplacements = plDisplacements;
    tDisplacements = ptDisplacements;
    epsilon = radius/200;
    if(timer->isActive()){
        timer->stop();
        animate = false;
    }
    else{
        timer->start(250);
        animate = true;
    }

}

void GLRail3d::SetSections(int sections){
    this->sections = sections;
    disp = InitDisplacements(boundary.count(),boundary[0].count(),sections);
    delete mesh;
    mesh = new MeshCuad(boundary,boundaryIndices,spacing,sections,disp);
    updateGL();
}

void GLRail3d::SetSpacing(double value){
    this->spacing = value;
    delete mesh;
    mesh = new MeshCuad(boundary,boundaryIndices,spacing,sections,disp);
    updateGL();
}

void GLRail3d::Center(){
    scale = 1;
    obj_pos[0] = center[0];
    obj_pos[1] = center[1];
    obj_pos[2] = center[2];
    updateGL();
}

void GLRail3d::LocateCenter(){
    float minX = std::numeric_limits<float>::max();
    float minY = minX;
    float maxX = -std::numeric_limits<float>::min();
    float maxY = maxX;

    for(int i = 0; i < mesh->BoundariesCount(); i++)
    {
        Section *section = mesh->BoundaryISectionJ(i,0);
        for(int k = 0; k < section->VerticesCount(); k++)
        {
            if(minX > section->Vertice(k)->X())
                minX = section->Vertice(k)->X();
            if(minY > section->Vertice(k)->Y())
                minY = section->Vertice(k)->Y();
            if(maxX < section->Vertice(k)->X())
                maxX = section->Vertice(k)->X();
            if(maxY < section->Vertice(k)->Y())
                maxY = section->Vertice(k)->Y();
        }
    }
    float maxZ = mesh->DelthaZ()*mesh->SectionsCount();
    center[0] = (maxX + minX) / 2-0.03;
    center[1] = -(maxY + minY) / 2;
    center[2] = maxZ / 2;

    radius = max(max((maxX - minX)/2, (maxY - minY)/2), maxZ);
    ratio = radius/10;
}

void GLRail3d::DisplayMesh(){
    for(int i = 0; i < mesh->BoundariesCount(); i++)
    {
        for(int j = 0; j < mesh->SectionsCount()-1; j++)
        {
            QList<Cuad> cuads = mesh->BoundaryICuadsSectionJ(i,j);
            for(int k = 0; k < cuads.count(); k++)
            {
                double a = cuads[k].Vertice(0).X();
                double b = cuads[k].Vertice(0).Y();
                double c = cuads[k].Vertice(0).Z();
                glColor3f (0.5,0.5,0.5);
                glBegin(GL_QUADS);
                    glVertex3f(cuads[k].Vertice(0).X()+cuads[k].Vertice(0).Displacement()->ax,
                               cuads[k].Vertice(0).Y()+cuads[k].Vertice(0).Displacement()->ay,
                               cuads[k].Vertice(0).Z()+cuads[k].Vertice(0).Displacement()->az);
                    glVertex3f(cuads[k].Vertice(1).X()+cuads[k].Vertice(1).Displacement()->ax,
                               cuads[k].Vertice(1).Y()+cuads[k].Vertice(1).Displacement()->ay,
                               cuads[k].Vertice(1).Z()+cuads[k].Vertice(1).Displacement()->az);
                    glVertex3f(cuads[k].Vertice(2).X()+cuads[k].Vertice(2).Displacement()->ax,
                               cuads[k].Vertice(2).Y()+cuads[k].Vertice(2).Displacement()->ay,
                               cuads[k].Vertice(2).Z()+cuads[k].Vertice(2).Displacement()->az);
                    glVertex3f(cuads[k].Vertice(3).X()+cuads[k].Vertice(3).Displacement()->ax,
                               cuads[k].Vertice(3).Y()+cuads[k].Vertice(3).Displacement()->ay,
                               cuads[k].Vertice(3).Z()+cuads[k].Vertice(3).Displacement()->az);
                 glEnd();
                 glColor3f (0.0,0.0,0.0);
                 glLineWidth(2.0);
                 glBegin(GL_LINE_LOOP);
                    glVertex3f(cuads[k].Vertice(0).X()+cuads[k].Vertice(0).Displacement()->ax,
                               cuads[k].Vertice(0).Y()+cuads[k].Vertice(0).Displacement()->ay,
                               cuads[k].Vertice(0).Z()+cuads[k].Vertice(0).Displacement()->az);
                    glVertex3f(cuads[k].Vertice(1).X()+cuads[k].Vertice(1).Displacement()->ax,
                               cuads[k].Vertice(1).Y()+cuads[k].Vertice(1).Displacement()->ay,
                               cuads[k].Vertice(1).Z()+cuads[k].Vertice(1).Displacement()->az);
                    glVertex3f(cuads[k].Vertice(2).X()+cuads[k].Vertice(2).Displacement()->ax,
                               cuads[k].Vertice(2).Y()+cuads[k].Vertice(2).Displacement()->ay,
                               cuads[k].Vertice(2).Z()+cuads[k].Vertice(2).Displacement()->az);
                    glVertex3f(cuads[k].Vertice(3).X()+cuads[k].Vertice(3).Displacement()->ax,
                               cuads[k].Vertice(3).Y()+cuads[k].Vertice(3).Displacement()->ay,
                               cuads[k].Vertice(3).Z()+cuads[k].Vertice(3).Displacement()->az);
                glEnd();
            }
        }
    }
}

void GLRail3d::ZoomIn(){
    scale *=1.2;
    updateGL();
}

void GLRail3d::ZoomOut(){
    scale /=1.2;
    updateGL();
}

void GLRail3d::resizeGL(int width, int height){
    glViewport(0, 0, width, height);
}

void GLRail3d::CleanControl(){
    if(!dmesh)
        return;
    radius = 0;
    scale = 1;
    dmesh = false;
    delete mesh;
    triangles.clear();
    boundary.clear();
    boundaryIndices.clear();
    x.clear();
    y.clear();
    disp.clear();
    updateGL();
}

void GLRail3d::SelectionChanged(){
    QList<QList<QList<GVector3> > > allDisplacements;
    if(animate){
        for(int z = 0 ; z < mesh->SectionsCount();z++)
        {
            QList<GVector3> newDisplacements;
            if(type == 0)
                newDisplacements = fDisplacements[row][col].GetDisplacements(z,epsilon,timeDiscretization);
            else if(type == 1)
                newDisplacements = lDisplacements[row][col].GetDisplacements(z,epsilon,timeDiscretization);
            else if(type == 2)
                newDisplacements = tDisplacements[row][col].GetDisplacements(z,epsilon,timeDiscretization);
            else
                newDisplacements = displacements[row][col].GetDisplacements(z,epsilon,timeDiscretization);
            QList<QList<GVector3> > result = BoundaryDisplacements(newDisplacements,mesh->indexes);
            for(int i = 0; i < result.count(); i++)
            {
                allDisplacements.append(QList<QList<GVector3> >());
                allDisplacements[i].append(result[i]);
            }
        }
        mesh->ChangeDisplacements(allDisplacements);
    }
}

void GLRail3d::UpdateRowCol(int prow, int pcol, int ptype){
    row = prow;
    col = pcol;
    type = ptype;

    if(row != -1 && col != -1)
        SelectionChanged();
    updateGL();
}

QList<QList<GVector3> > GLRail3d::BoundaryDisplacements(QList<GVector3> allDisplacements,
                                                          QList<QList<int> > boundaryIndexes)
{
    QList<QList<GVector3> > result;
    for(int i = 0 ; i < boundaryIndexes.count();i++)
    {
        result.append(QList<GVector3>());
        for(int j = 0; j < boundaryIndexes[i].count(); j++)
        {
            if(boundaryIndexes[i][j] >= allDisplacements.count())
            {
                result[i].append(GVector3(0,0,0));
            }
            else result[i].append(allDisplacements[boundaryIndexes[i][j]]);
        }
    }
    return result;
}
void GLRail3d::timeElapsed()
{
    if(timeDiscretization > mesh->SectionsCount()-1)
        timeDiscretization = 0;
    if(row != -1 && col != -1)
    {
        SelectionChanged();
    }
    else
        mesh->NoDisplacements();
    timeDiscretization++;
    updateGL();
}

void GLRail3d::wheelEvent(QWheelEvent *wheelEvent){
    (wheelEvent->delta() > 0)?ZoomIn():ZoomOut();
}

void GLRail3d::mouseMoveEvent(QMouseEvent *mouseEvent){
    oldMovementPosition = newMovementPosition;
    newMovementPosition = mouseEvent->pos();

    float realX = ((GLfloat)(newMovementPosition.x() - oldMovementPosition.x())) / (float)width();
    float realY = ((GLfloat)(newMovementPosition.y() - oldMovementPosition.y())) / (float)height();
   // Change distance
    if (leftButton)
    {
        obj_pos[0] += realX*radius;
        obj_pos[1] -= realY*radius;
        azimuth_ += 180*realX;
        elevation_ += 180*realY;
    }
    else
    {



    }
    updateGL();

}

void GLRail3d::mousePressEvent(QMouseEvent *mousePressed){
    newMovementPosition = mousePressed->pos();
    if(mousePressed->button() == Qt::LeftButton)
        leftButton = true;
    else leftButton = false;
}

void GLRail3d::keyPressEvent(QKeyEvent *keyEvent){
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

/*void GLRail3d::keyPressEvent(QKeyEvent *keyEvent)
{
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
*/
void GLRail3d::MoveUp()
{
    obj_pos[1] += ratio;
    updateGL();
}

void GLRail3d::MoveDown()
{
    obj_pos[1] -= ratio;
    updateGL();
}
void GLRail3d::MoveLeft()
{
    obj_pos[0] -= ratio;
    updateGL();
}
void GLRail3d::MoveRight()
{
    obj_pos[0] += ratio;
    updateGL();
}
