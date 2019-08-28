/*:
**     adapt.cpp - Adaptation scheme
**
**   Copyright (C) 2004 Luiz Velho.
*/

#include "a48.h"

using namespace A48;

/**
* Refines a region specified by a half-edge, recursively.
* @param e Specifies half-edge that must be splitted.
* @return vertex of the half-edge splitted.
*/
Vertex* Mesh::refine(Hedge *e)
{
  Face *f; Hedge *r[2]; int n = 0;
  if ((f = e->face()) != NULL && f->subd_edge() != e)
    r[n++] = f->subd_edge();
  if ((f = e->mate()->face()) != NULL && f->subd_edge() != e->mate())
    r[n++] = f->subd_edge();
  for (int k=0; k < n; k++)
    refine(r[k]);
  update_ref_front(e->edge());
  Vertex* w = split(e);
  update_ref_front(w);
  return w;
}

/**
* Simplifies a region specified by a vertex.
* @param w Specifies a vertex that must be removed from the mesh.
* @return half-edge created by welding two half-edges.
*/
Hedge* Mesh::simplify(Vertex *w)
{
  int n = 0, weld_deg = (w->is_bdry()) ? 3 : 4;
  do {
    int lmax = w->level(); Hedge *e; Vertex *u, *v;
    for (e = w->star_first(), n = 0; e != NULL; e = w->star_next(e), n++) {
      u = e->org();
      if (u->level() > lmax) {
	lmax = u->level(); v = u;
      }
    }
    if (lmax > w->level())
      simplify(v);
  } while (n > weld_deg);
  update_simpl_front(w);
  Hedge* e = weld(w);
  update_simpl_front(e->edge());
  return e;
}

/**
* Updates the refinement queue and inserts a half-edge based on rank.
* The adge is inserted in the queue if its rank is less than a threshold.
* @param t Specifies a priority threshold.
*/
void Mesh::adapt_refine(double t)
{
  for (EdgeIter ei = ec_.begin(); ei != ec_.end(); ei++)
    if ((*ei)->is_in_heap())
      rf_.update((MxHeapable*)(*ei), surf()->ref_rank(*ei));
  Edge *e;
  while ((e = (Edge*) rf_.extract()) != NULL) {
   if (e->heap_key() < t) {
     rf_.insert((MxHeapable*)e, surf()->ref_rank(e));
     return;
   }
   refine(e->hedge(0));
  }
} 


/**
* Updates the simplifying queue and inserts a vertex based on rank.
* The edge is inserted in the queue if its rank is less than a threshold.
* @param t Specifies a priority threshold.
*/
void Mesh::adapt_simplify(double t)
{
  for (VertexIter vi = vc_.begin(); vi != vc_.end(); vi++) 
    if ((*vi)->is_in_heap())
      sf_.update((MxHeapable*)(*vi), surf()->simpl_rank(*vi));
  Vertex *v;
  while ((v = (Vertex*) sf_.extract()) != NULL) {
    if (v->heap_key() < t) {
      sf_.insert((MxHeapable*)v, surf()->simpl_rank(v));
      return;
    }
    simplify(v);
  }
}
