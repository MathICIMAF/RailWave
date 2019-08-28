/*:
**     face.cpp - Face element
**
**  Copyright (C) 2004 Luiz Velho.
*/


#include <algorithm>
#include "a48.h"

/* indices of 3-ring */
#define NEXT3(k) ((k < 2)? (k+1) : 0) ///< Returns next half-edge index.
#define PREV3(k) ((k > 0)? (k-1) : 2) ///< Returns previous half-edge index.

using namespace std;
using namespace A48;

/**
* Updates the half-edges of the face. 
* Make half-edge links to build a face loop.
* Links half-edges to face. 
* @param e0 First half-edge.
* @param e1 Second half-edge.
* @param e2 Third half-edge.
* @return face pointer.
*/
Face* Face::reuse(Hedge *e0, Hedge *e1, Hedge *e2)
{
  e_ = e0;
  e0->set_next(e1); e1->set_next(e2); e2->set_next(e0);
  e0->set_face(this); e1->set_face(this); e2->set_face(this);
  return this;
}
/**
* Returns the half-edge specified by the index k.
* @param k half-edge index. Valid values of k=[0,1,2].
* @return half-edge pointer.
*/
Hedge* Face::hedge(int k)
{
  switch (k) {
  case 0: return e_;
  case 1: return e_->next();
  case 2: return e_->next()->next();
  }
  throw Error("hedge index");
}

/**
* Returns the opposite vertex of the half-edge specified by index k.
* @param k half-edge index. Valid values of k=[0,1,2].
* @return vertex pointer. 
*/
Vertex* Face::vertex(int k)
{
  return hedge(PREV3(k))->org();
}

/**
* Updates the half-edge of index k.
* @param k half-edge index, Valid values of k=[0,1,2].
* @param h Specify a new half-edge.
*/
void  Face::set_hedge(int k, Hedge* h)
{
  Hedge *n = hedge(NEXT3(k));
  Hedge *p = hedge(PREV3(k));
  h->set_next(n);
  p->set_next(h);
  if (k == 0) e_ = h;
}

/**
* Links vertices to half-edges.
* Builds canonical vertex indexing to work with stelar mesh.
*/
void Face::link_star_verts()
{
  for (int k = 0; k < 3; k++) 
    vertex(k)->set_star(hedge(NEXT3(k)));
}

/**
* Updates the opposite vertex of half-edge specified by the index k.
* @param k índex of half-edge. Valid values of k=[0,1,2].
* @param v Specifies the new vertex.
*/
void  Face::set_vertex(int k, Vertex* v)
{
  hedge(PREV3(k))->set_org(v);
}

/**
* Returns the face level. This level is determined by the highest
* resolution level of face vertices.
* @return face's resolution level.
*/
int Face::level()
{
  return max(vertex(0)->level(), max(vertex(1)->level(),vertex(2)->level()));
}

/**
* Returns true if the face belongs to level 0, else return false.
* @return boolean type (true or false).
*/
bool Face::is_inbase()
{
  return (vertex(0)->level() == 0 
	  && vertex(1)->level() == 0
	  && vertex(2)->level() == 0);
}

