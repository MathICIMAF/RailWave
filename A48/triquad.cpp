/*
**   triquad.cpp - Make triquad mesh
**
**   Copyright (C) 2004 Luiz Velho.
*/

#include <vector>
#include <stack>
#include <queue>

#include "a48.h"

using namespace A48;

/**
* @struct node
* @brief Associates the edge to a node. 
*
* This association is necessary to make possible evaluate the edge size. 
* The estimation is made by the surface class.
*/
struct node {
    double lenght; ///< Storage the edge length.
    Edge *edge; ///< Edge pointer.

   /** @brief Initializes the variables of the class node.
     * @param e Edge pointer.
     * @param m Mesh pointer.*/
    node( Edge *e, Mesh *m ) {
        edge = e;
        lenght = m->surf()->elenght(e);
    }

   /** @brief Returns true if the edge length is smaller 
    * than the received egde legth.
    * @param a Node structure.*/
    bool operator<( const node &a ) const {
        return lenght < a.lenght;
    }
};

/**
* Updates mark field of the next and previus edges.
* @param e Points to edge that was processed.
*/
void mark_link(Edge *e)
{
  Hedge *h0 = e->hedge(0);
  Hedge *h1 = e->hedge(1);
  if (h0->face() != NULL) {
    h0->next()->edge()->set_mark(true);
    h0->prev()->edge()->set_mark(true);
  }
  if (h1->face() != NULL) {
    h1->next()->edge()->set_mark(true);
    h1->prev()->edge()->set_mark(true);
  }
}

/**
* Transforms a general mesh into a triangulated-quadrangulation.
*/
void Mesh::make_triquad(void)
{
  std::priority_queue< node, std::vector< node >, std::less<node> > q;

  for (VertexIter p = vc_.begin(); p != vc_.end(); p++)
    (*p)->set_level(0);
  for (EdgeIter e = ec_.begin(); e != ec_.end(); e++) {
    (*e)->set_mark(false);
    q.push(node(*e, this));
  }
  while ( !q.empty() ) {
    node n = q.top(); q.pop();
    if (!n.edge->is_marked()) {
      mark_link(n.edge);
      if (n.edge->hedge(0)->face() != NULL)
	split(n.edge->hedge(0));
      else
	split(n.edge->hedge(1));
    }
  }
  for (FaceIter f = fc_.begin(); f != fc_.end(); f++)
    if ((*f)->level() == 0)
      split(*f);
  for (VertexIter p = vc_.begin(); p != vc_.end(); p++)
    (*p)->set_level(0);
}


/**
* Verifies if a mesh is a triangulated-quadrangulation.
*/
bool Mesh::is_triquad(void)
{
  for (EdgeIter ei = ec_.begin(); ei != ec_.end(); ei++) {
    Edge *e = *ei;
    if (!e->is_bdry()) {
      bool f0s = (e->hedge(0)->face()->subd_edge() == e->hedge(0));
      bool f1s = (e->hedge(1)->face()->subd_edge() == e->hedge(1));
      if ((f0s && !f1s) || (!f0s && f1s))
	return false;
    }
  }
  return true;
}

