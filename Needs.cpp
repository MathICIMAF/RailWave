#include "Needs.h"

bool Equals(Vertex *v0, Vertex *v1)
{
    Pvert *pv0 = Pvert::cast(v0);
    Pvert *pv1 = Pvert::cast(v1);

    double x0 = pv0->a.g[0];
    double y0 = pv0->a.g[1];
    double z0 = pv0->a.g[2];
    double x1 = pv1->a.g[0];
    double y1 = pv1->a.g[1];
    double z1 = pv1->a.g[2];

    return (x0 == x1 && y0 == y1 && z0 == z1);
}
QList<Vertex*> Get1Neighbors(Vertex *v)
{
    QList<Vertex*> result;
    Hedge *e = v->star_first();
    Hedge *el = v->star_first();
    Hedge *er = v->star_first();
    Pvert *eorg = Pvert::cast(e->org());
    Pvert *edtn = Pvert::cast(e->dst());

    if(Equals(v, e->org()))
    {
        e = e->mate();
        el = el->mate();
        er = er->mate();
    }

    //Move Left
    eorg = Pvert::cast(el->org());
    edtn = Pvert::cast(el->dst());
    while(result.count() == 0 || !Equals(result.at(0),el->org()))
    {
        result.append(el->org());
        if(el->mate()->face() != NULL)
        {
            el = el->mate()->prev();
            eorg = Pvert::cast(el->org());
            edtn = Pvert::cast(el->dst());
        }
        else break;
    }
    //Move Right
    er = er->next();
    eorg = Pvert::cast(er->org());
    edtn = Pvert::cast(er->dst());
    while(result.count() == 0 || !Equals(result.at(0),er->dst()))
    {
        result.append(er->dst());
        if(er->mate()->face() != NULL)
        {
            er = er->mate()->next();
            eorg = Pvert::cast(er->org());
            edtn = Pvert::cast(er->dst());
        }
        else break;
    }
    return result;
}
QString GetDirectoryAddress(QString path)
{
    int index = path.lastIndexOf('/');
    QString result = "";
    for(int i = 0; i < index; i++)
        result += path[i];
    return result;
}
void WriteOFFMesh(Mesh* mesh, QString path)
{
    QFile file(path+".off");

    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save .cuad Mesh");
        return;
    }

    QTextStream out(&file);

    //1- write OFF, number of vertices, number of faces
    out << "OFF \n";
    out << mesh->num_verts() << " " << mesh->num_faces() << " " << mesh->num_edges() << "\n";

    int id = 1;
    //2- write vertices coordinates
    for(VertexIter vIter = mesh->verts_begin(); vIter != mesh->verts_end(); vIter ++)
    {
        Vertex* vert = (*vIter);
        vert->setId(id++);
        Pvert* v = Pvert::cast(vert);
        out << v->a.g.ax << " " << v->a.g.ay << "\n";
    }

    for(FaceIter fIter = mesh->faces_begin(); fIter != mesh->faces_end(); fIter ++)
    {
            Face *f = *(fIter);
            f->setSaved(-1);
    }

    //3- write faces
    for(FaceIter fIter = mesh->faces_begin(); fIter != mesh->faces_end(); fIter ++)
    {
            Face *f = *(fIter);
            if(f->getSaved() == 1)continue;

            Pvert *v0 = Pvert::cast(f->vertex(1));
            Pvert *v1 = Pvert::cast(f->vertex(2));
            Pvert *v2 = Pvert::cast(f->vertex(3));

            out << "3 " << v0->getId()<< " ";
            out << v1->getId() << " ";
            out << v2->getId() << "\n";
            f->setSaved(1);

            Hedge *hedge = f->hedge(0);
            if(hedge->mate() != NULL)
            {
                Face *fMate = hedge->mate()->face();
                if(fMate!=NULL)
                {
                    if(fMate->getSaved() == 1)
                        continue;
                    v0 = Pvert::cast(fMate->vertex(1));
                    v1 = Pvert::cast(fMate->vertex(2));
                    v2 = Pvert::cast(fMate->vertex(3));

                    out << "3 " << v0->getId() << " ";
                    out << v1->getId() << " ";
                    out << v2->getId() << "\n";
                    fMate->setSaved(1);
                }
            }
        }

    file.close();

    /*QMessageBox mbx;
    mbx.setText("OFF mesh saved");
    mbx.show();*/
}
QList<QList<Vertex*> > GetBoundary(Mesh *s)
{
    //1- Lista de vertices de la frontera
    QList<QList<Vertex*> > result;
    QList<Vertex*> allBoundVertices;
    //2- Buscar un vertice v0 que sea frontera.
    VertexIter v = s->verts_begin();
    while(v != s->verts_end())
    {
        if((*v)->is_bdry())
            allBoundVertices.append(*v);
        v++;
    }
    //3- Si ningun vertice es frontera devolver la lista vacia.
    if (allBoundVertices.size() == 0) return result;

    //4- Buscar el resto de los vertices
    while(allBoundVertices.count()!=0)
    {
        QList<Vertex*> bound;
        bound.append(allBoundVertices[0]);
        Vertex *v0 = 0;
        do
        {
            v0 = bound[bound.count()-1];
            v0->setId(1);
            QList<Vertex*> neighbors = Get1Neighbors(v0);
            for(int i = 0;i<neighbors.count();i++)
            {
                Vertex *vi = neighbors[i];
                if (vi->is_bdry() && vi->getId()!=1)
                {
                    bound.append(vi);
                    break;
                }
            }
        }
        while (bound[bound.count()-1] != v0);
        result.append(bound);
        allBoundVertices = GetDiff(allBoundVertices,bound);
    }

    if(result.count() > 1)
        result = SortBoundaries(result);
    return result;
}
int Contains(QList<Edge*> vect, Edge *hedge)
{
    for(int i = 0; i < vect.count();i++)
    {
        if(vect[i]->org()->getId() == hedge->org()->getId() && vect[i]->dst()->getId() == hedge->dst()->getId())
            return i;
    }
    return -1;
}
void WriteOterosMesh(Mesh *mesh, QString path)
{
    QString vertexPath = path;
    QString facesPath = path;

    vertexPath += "vertex.txt";
    facesPath += "elems.txt";

    QFile file(vertexPath);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save .cuad Mesh");
        return;
    }

    QTextStream out(&file);
    int id = 1;
    //2- write vertices coordinates
    for(VertexIter vIter = mesh->verts_begin(); vIter != mesh->verts_end(); vIter ++)
    {
        Vertex* vert = (*vIter);
        vert->setId(id++);
        Pvert* v = Pvert::cast(vert);
        out << v->getId() << " " << v->a.g.ax << " " << v->a.g.ay << "\n";
    }
    file.close();

    QFile file2(facesPath);
    if(!file2.open(QFile::WriteOnly))
    {
        printf("Fail To Save .cuad Mesh");
        return;
    }

    QTextStream out2(&file2);
    int count = 1;
    for(FaceIter fIter = mesh->faces_begin(); fIter != mesh->faces_end(); fIter ++)
    {
        Face *f = *(fIter);
        f->setSaved(-1);
    }

    //3- write faces
    for(FaceIter fIter = mesh->faces_begin(); fIter != mesh->faces_end(); fIter ++)
    {
        Face *f = *(fIter);
        if(f->getSaved() == 1)continue;

        Pvert *v0 = Pvert::cast(f->vertex(3));
        Pvert *v1 = Pvert::cast(f->vertex(2));
        Pvert *v2 = Pvert::cast(f->vertex(1));

        out2 << count++ << " " << v2->getId() << " " << v1->getId() << " " << v0->getId() << "\n";
        f->setSaved(1);

        Hedge *hedge = f->hedge(0);
        if(hedge->mate() != NULL)
        {
            Face *fMate = hedge->mate()->face();
            if(fMate!=NULL)
            {
                if(fMate->getSaved() == 1)
                    continue;
                v0 = Pvert::cast(fMate->vertex(3));
                v1 = Pvert::cast(fMate->vertex(2));
                v2 = Pvert::cast(fMate->vertex(1));

                out2 << count++ << " " << v2->getId() << " " << v1->getId() << " " << v0->getId() << "\n";
                fMate->setSaved(1);
            }
        }
    }
    file2.close();

}
int Contains(QList<Pvert*> boundary, Pvert *vert)
{

    for(int i=0;i<boundary.count();i++)
    {
        if(boundary[i]->a.g[0] == vert->a.g[0] && boundary[i]->a.g[1] == vert->a.g[1] && boundary[i]->a.g[2] == vert->a.g[2])
            return i;
    }
    return -1;
}
double GetRadius(QList<double> x, QList<double> y)
{
    double max = x.at(0);
    for(int index = 0; index < x.count(); index++)
    {
        if(max < x.at(index))
            max = x.at(index);
        if(max < y.at(index))
            max = y.at(index);
    }
    return max;
}
int Contains(QList<Vertex*> boundary, Vertex *vert)
{
    Pvert *v1 = Pvert::cast(vert);
    for(int i=0;i<boundary.count();i++)
    {
        Vertex *ver = boundary[i];
        Pvert *v = Pvert::cast(ver);
        if(v->a.g[0] == v1->a.g[0] && v->a.g[1] == v1->a.g[1] && v->a.g[2] == v1->a.g[2])
            return i;
    }
    return -1;
}
QList<Vertex*> Revert(QList<Vertex*> list, int index)
{
    QList<Vertex*> result;
    for(int i = index; i >= 0; i--)
        result.append(list[i]);
    for(int i = list.count()-1; i>index; i--)
        result.append(list[i]);
    return result;
}
QList<Vertex*> Shift(QList<Vertex*> list, int beginIndex)
{
    QList<Vertex*> result;

    for(int i = beginIndex; i<list.count(); i++)
        result.append(list[i]);
    for(int i = 0; i<beginIndex; i++)
        result.append(list[i]);
    return result;
}
QPoint Contains(QList<QList<Vertex*> > boundary, Vertex *vert)
{
    Pvert *v1 = Pvert::cast(vert);
    for(int j = 0; j< boundary.count();j++)
    {
        for(int i=0;i<boundary[j].count();i++)
        {
            Vertex *ver = boundary[j][i];
            Pvert *v = Pvert::cast(ver);
            if(v->a.g[0] == v1->a.g[0] && v->a.g[1] == v1->a.g[1] && v->a.g[2] == v1->a.g[2])
                return QPoint(j,i);
        }
    }
    return QPoint(-1,-1);
}
QList<QList<Vertex*> > SortBoundaries(QList<QList<Vertex*> > list)
{
    //outer-inner boundary Tube case
    if(list.count() == 2)
    {
        Pvert *pv1 = Pvert::cast(list[0][0]);
        Pvert *pv11 = Pvert::cast(list[0][1]);
        Pvert *pv2 = Pvert::cast(list[1][0]);
        Pvert *pv21 = Pvert::cast(list[1][1]);

        if(sqrt(pow(pv1->a.g.ax-pv11->a.g.ax,2) + pow(pv1->a.g.ay - pv11->a.g.ay,2)) <
           sqrt(pow(pv2->a.g.ax-pv21->a.g.ax,2) + pow(pv2->a.g.ay - pv21->a.g.ay,2)))
        {
            QList<Vertex*> temp = list[0];
            list[0] = list[1];
            list[1] = temp;
        }
        //Look for the point that has x = 0 && y >0
        int innerZeroIndex = -1;
        int outerZeroIndex = -1;
        for(int i = 0; i < list[0].count(); i++)
        {
            pv1 = Pvert::cast(list[0][i]);
            pv2 = Pvert::cast(list[1][i]);
            if(pv1->a.g.ax == 0 && pv1->a.g.ay > 0)
                outerZeroIndex = i;
            if(pv2->a.g.ax == 0 && pv2->a.g.ay > 0)
                innerZeroIndex = i;
            if(innerZeroIndex != -1 && outerZeroIndex != -1)
                break;
        }

        //Be sure that boundaries are in the correct order
        Pvert *aux1 = Pvert::cast(list[0][outerZeroIndex]);
        Pvert *aux2;
        if(outerZeroIndex == list[0].count()-1)
            aux2 = Pvert::cast(list[0][0]);
        else aux2 = Pvert::cast(list[0][outerZeroIndex+1]);
        if(aux1->a.g.ax > aux2->a.g.ax)
            list[0] = Revert(list[0], outerZeroIndex);
        else
        {
            if(outerZeroIndex != 0)
                list[0] = Shift(list[0],outerZeroIndex);
        }

        aux1 = Pvert::cast(list[1][innerZeroIndex]);
        if(innerZeroIndex == list[1].count()-1)
            aux2 = Pvert::cast(list[1][0]);
        else aux2 = Pvert::cast(list[1][innerZeroIndex+1]);
        if(aux1->a.g.ax > aux2->a.g.ax)
            list[1] = Revert(list[1],innerZeroIndex);
        else
        {
            if(innerZeroIndex != 0)
                list[1] = Shift(list[1],innerZeroIndex);
        }
    }
    return list;
}
void WriteBoundary(QString path, QList<int>index, QList<Vertex*> boundary)
{
    QFile file(path);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save Boundary");
        return;
    }

    QTextStream out(&file);
    out<<boundary.count()<< "\n";
    for(int i = 0; i < boundary.count();i++)
    {
        Vertex *vert = boundary[i];
        Pvert *v = Pvert::cast(vert);
        out << index[i] << " " << v->a.g.ax << " " << v->a.g.ay << "\n";
    }
    file.close();
}
QList<Vertex*> GetDiff(QList<Vertex*> allBoundVertices, QList<Vertex*> bound)
{
    QList<Vertex*> diff;
    for(int i = 0; i < allBoundVertices.count(); i++)
    {
        if(!bound.contains(allBoundVertices[i]))
            diff.push_back(allBoundVertices[i]);
    }

    return diff;
}
void WriteDATMesh(Mesh *mesh, QString path, ShapeType shape, double ir, double thickness)
{
    QFile file(path+".dat");

    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save .dat Mesh");
        return;
    }
    QList<QList<Vertex*> > boundary = GetBoundary(mesh);
    QList<QList<int> > index;
    for(int i = 0;i<boundary.count();i++)
    {
        index.append(QList<int>());
        for(int j = 0; j<boundary[i].count(); j++)
            index[i].append(0);
    }

    QTextStream out(&file);
    //1- write DAT, shape, ir, thickness, number of vertices, number of faces
    out << "DAT\n";
    out << shape << " " << ir << " " << thickness << "\n";
    out << mesh->num_verts() << " " << mesh->num_faces() << " " << mesh->num_edges() << "\n";

    int id = 1;
    //2- write vertices coordinates
    for(VertexIter vIter = mesh->verts_begin(); vIter != mesh->verts_end(); vIter ++)
    {
        Vertex *vert = (*vIter);
        vert->setId(id++);
        Pvert *v = Pvert::cast(vert);
        out << v->a.g.ax << " " << v->a.g.ay << "\n";

        QPoint aux = Contains(boundary, vert);
        if(aux.x() != -1)
        {
            const int s = vert->getId();
            index[aux.x()][aux.y()] = s;
        }
    }

    for(FaceIter fIter = mesh->faces_begin(); fIter != mesh->faces_end(); fIter ++)
    {
            Face *f = *(fIter);
            f->setSaved(-1);
    }

    //3- write faces
    for(FaceIter fIter = mesh->faces_begin(); fIter != mesh->faces_end(); fIter ++)
    {
            Face *f = *(fIter);
            if(f->getSaved() == 1)continue;

            Pvert *v0 = Pvert::cast(f->vertex(1));
            Pvert *v1 = Pvert::cast(f->vertex(2));
            Pvert *v2 = Pvert::cast(f->vertex(3));
            out << v0->getId()<< " ";
            out << v1->getId() << " ";
            out << v2->getId() << "\n";
            f->setSaved(1);

            Hedge *hedge = f->hedge(0);
            if(hedge->mate() != NULL)
            {
                Face *fMate = hedge->mate()->face();
                if(fMate!=NULL)
                {
                    if(fMate->getSaved() == 1)
                        continue;
                    v0 = Pvert::cast(fMate->vertex(1));
                    v1 = Pvert::cast(fMate->vertex(2));
                    v2 = Pvert::cast(fMate->vertex(3));

                    out << v0->getId() << " ";
                    out << v1->getId() << " ";
                    out << v2->getId() << "\n";
                    fMate->setSaved(1);
                }
            }
        }

    out << "boundary " << index.count() << "\n";
    for(int i = 0; i < index.count();i++)
    {
        out << index.at(i).count() << "\n";
        for(int j = 0; j < index.at(i).count(); j++)
            out << index.at(i).at(j) << "\n";
    }

    file.close();
    QMessageBox mbx;
    mbx.setText("DAT mesh saved");
    mbx.show();
}

