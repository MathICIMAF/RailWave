/*:
**     mesh.cpp - mesh ops
**
**   Copyright (C) 2004 Luiz Velho.
*/

#include "a48.h"
#include <map>
//#include "Mesh.h"

using namespace A48;


/**
* Receives a surface object that determines the mesh.
* @param s Surface pointer.
*/
Mesh::Mesh(Surface* s,QObject* parent) : surf_(s),QObject(parent)
{
  int np, nf, *fcs;
  surf()->base_mesh(&np, &fcs, &nf);
  set_base_mesh(np, fcs, nf);
  if (!is_triquad())
    make_triquad();
  new_front();
}

/**
* Dealocates all data (vertices, faces, edges and half-edges) of the mesh.
*/
Mesh::~Mesh()
{
  for (VertexIter p = vc_.begin(); p != vc_.end(); p++) 
    delete *p;
  vc_.clear();
  for (EdgeIter e = ec_.begin(); e != ec_.end(); e++)
    delete *e;
  ec_.clear();
  for (FaceIter f = fc_.begin(); f != fc_.end(); f++)
    delete *f;
  fc_.clear();
}

/**
* Adds a new face in the mesh
* @param e0 First half-edge.
* @param e1 Second half-edge.
* @param e2 Third half-edge.
* @return pointer to new face.
*/
Face* Mesh::add_face(Hedge *e0, Hedge *e1, Hedge *e2)
{
  Face *f = new Face(e0, e1, e2);
  add_face(f);
  return f;
}

/**
* Adds a new edge in the mesh.
* @param v0 origin vertex pointer.
* @param v1 end vertex pointer.
* @return pointer to half-egde of the origin vertex.
*/
Hedge* Mesh::add_edge(Vertex *v0, Vertex *v1)
{
  Edge *e = new Edge(v0, v1);
  add_edge(e);
  return e->hedge(0);
}

/**
* Adds a new vertex in the mesh.
* @return pointer to new vertex.
*/
Vertex* Mesh::add_vertex(void)
{
  Vertex *v = surf()->new_vertex();
  add_vertex(v);
  return v;
}

/**
* Links stars of mesh
*/
void Mesh::link_mesh()
{
  for (FaceIter f = fc_.begin(); f != fc_.end(); f++)
    (*f)->link_star_verts();
  for (EdgeIter e = ec_.begin(); e != ec_.end(); e++) {
    Hedge* h = (*e)->hedge(0);
    if (h->is_bdry() && h->face() != NULL)
      h->dst()->set_star(h);
    else if (h->face() == NULL)
      h->org()->set_star(h->mate());
  }
}

/**
* Returns a half-edge specified by vertices of index i0 and i1.
* If the vertices don't exist, a new edge is created with two new vertices. 
* @param i0 First vertex index.
* @param i1 Second vertex index.
* @param verts Points to vertices.
* @param hedges Belongs to std::map class. (This class provides efficient search operations.)
* @return half-edge pointer.
*/
Hedge* Mesh::get_hedge(int i0, int i1, Vertex* verts[], HedgeMap *hedges)
{
  bool mate = false;
  if (i0 > i1) {
    std::swap<int>(i0, i1);
    mate = true;
  }
  Hedge* he; HedgeMap::iterator ei  = hedges->find(Ipair(i0,i1));
  if (ei == hedges->end())
    (*hedges)[Ipair(i0,i1)] = he = add_edge(verts[i0], verts[i1]);
  else
    he = (*ei).second;
  return (mate)? he->mate() : he ;
}

/**
* Adds a new face specified by the vertices of index i1, i2, i3. 
* If any vertex doesn't exist, than it is created.
* @param i0 First vertex index.
* @param i1 Second vertex index.
* @param i2 Third vertex index.
* @param verts Points to new vertexes.
* @param hedges Points to std::map class. (This class provide efficient search operations.)
*/
void Mesh::put_face(int i0, int i1, int i2, Vertex* verts[], HedgeMap* hedges)
{
  Hedge* e0 = get_hedge(i1, i2, verts, hedges);
  Hedge* e1 = get_hedge(i2, i0, verts, hedges);
  Hedge* e2 = get_hedge(i0, i1, verts, hedges);
  add_face(e0, e1, e2);
}

/**
* Receives all data information to build the base mesh.
* @param np Number of vertices.
* @param tris Array of indices with mesh topology.
* @param nf Number of faces.
*/
void Mesh::set_base_mesh(int np, int* tris, int nf)
{
  Vertex** verts = new Vertex*[np];
  for (int i=0; i < np; i++) {
    verts[i] = surf()->new_vertex();
    surf()->sample(i, verts[i]);
    add_vertex(verts[i]);
  }
  HedgeMap hedges;
  for (int i=0; i < nf; i++, tris += 3) 
    put_face(tris[2], tris[0],tris[1], verts, &hedges);
  link_mesh();
  delete(verts);
}
