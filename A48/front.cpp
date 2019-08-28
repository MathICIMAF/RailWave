/*:
**     front.cpp - Adaptation Front
**
**   Copyright (C) 2004 Luiz Velho.
*/
 
#include "a48.h"
  
using namespace A48;

/**
* Creates simplification and refinement fronts.
* Insert vertices that can be removed from mesh 
* and edges that can be splitted.
*/
void Mesh::new_front()
{
  while (rf_.extract() != NULL)
    ;
  for (EdgeIter ei = edges_begin(); ei != edges_end(); ei++) {
    if ((*ei)->is_subdiv())
      rf_.insert((MxHeapable *) *ei, surf()->ref_rank(*ei));
  }
  while (sf_.extract() != NULL)
    ;
  for (VertexIter w = verts_begin(); w != verts_end(); w++)
    if ((*w)->is_weld())
      sf_.insert((MxHeapable *) *w, surf()->simpl_rank(*w));
}

/**
 * Updates the refinement front (deletion phase)
* @param e Pointer to splitted edge.
*/
void Mesh::update_ref_front(Edge *e)
{
  rf_.remove(e);
  if (e->hedge(0)->next() != NULL) {
    Face *f = e->hedge(0)->face();
    Vertex *v = e->hedge(0)->next()->dst();
    if (!v->is_weld(f)) 
      sf_.remove(v);
  }
  if (e->hedge(1)->next() != NULL) {
    Face *f = e->hedge(1)->face();
    Vertex *v = e->hedge(1)->next()->dst();
    if (!v->is_weld(f))
      sf_.remove(v);
  }
}

/**
* Updates the refinement front (insertion phase)
* @param v Pointer to a vertex added.
*/
void Mesh::update_ref_front(Vertex *v)
{
  sf_.insert(v, surf()->simpl_rank(v));
  for (Hedge* s = v->star_first(); s != NULL; s = v->star_next(s))
    if (s->prev() != NULL) {
      Edge *e = s->prev()->edge();
      if (e->is_subdiv()) 
	rf_.insert(e, surf()->ref_rank(e));
    }
}

/**
* Updates the siplification front (deletion phase)
* @param v Pointer to vertex that must be removed from mesh.
*/
void Mesh::update_simpl_front(Vertex *v)
{
  sf_.remove(v);
  for (Hedge* s = v->star_first(); s != NULL; s = v->star_next(s))
    if (s->prev() != NULL) {
      Hedge *m = s->prev()->mate();
      if (!m->is_subdiv_itself()) {
	Edge *e = s->prev()->edge();
	rf_.remove(e);
      }
    }
}

/**
* Updates the siplification front (insertion phase)
* @param e Weld edge pointer.
*/
void Mesh::update_simpl_front(Edge *e)
{
  rf_.insert(e, surf()->ref_rank(e));
  if (e->hedge(0)->next() != NULL) {
    Vertex* w = e->hedge(0)->next()->dst();
    if (w->is_weld())
      sf_.insert(w, surf()->simpl_rank(w));
  }
  if (e->hedge(1)->next() != NULL) {
    Vertex* w = e->hedge(1)->next()->dst();
    if (w->is_weld())
      sf_.insert(w, surf()->simpl_rank(w));
  }
}
