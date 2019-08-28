/** 
* @file edge.h
* @brief Definition of half-edge and edge classes.
*/

/* Copyright (C) 2004 Luiz Velho. */

#ifndef EDGE_H
#define EDGE_H

#include <QObject>

namespace A48 {

/**
* @class Hedge
* @brief The Hedge class represents Half-edge topological entity of mesh 4-8.
*/
  class Hedge:public QObject {
  friend class Mesh;
  friend class Face;
  friend class Edge;

  Vertex    *o_; ///< Origin vertex pointer.
  Hedge     *n_; ///< Next half-edge pointer.
  Face      *f_; ///< Face pointer.
  Edge      *e_; ///< Edge pointer.

 public:

  Hedge(QObject *parent=0):QObject(parent){}
  /** @brief Returns the half-egde's face.
    * @return face pointer.*/
  Face* face() { return f_; };

  /** @brief Returns the half-edge's origin vertex.
    * @return vertex pointer.*/
  Vertex* org() { return o_; };

  /** @brief Returns the end vertex of the half-edge.
    * @return vertex pointer.*/
  Vertex* dst() { return mate()->org(); };

  /** @brief Returns the half-edges's parent edge.
    * @return edge pointer.*/
  Edge*  edge() { return e_; };

  /** @brief Returns the previous half-edge.
    * @return half-edge pointer.*/
  Hedge* prev() { return (n_)? n_->next() : NULL; };

  /** @brief Returns the next half-edge.
    * @return half-edge pointer.*/
  Hedge* next() { return n_; };

  /** @brief Returns the opposite half-edge.
    * @return half-edge pointer.*/
  Hedge* mate();

  /** @brief Returns half-edge resolution level.*/
  int level();

  /** @brief Returns true if a half-edge belongs to mesh's boundary, else return false.*/
  bool is_bdry();

  /** @brief Returns true if the half-edge is a split half-edge, else return false. 
    * @return true if split half-edge, else return false.*/
  bool is_subdiv_itself() { 
    return (face() != NULL && face()->subd_edge() == this); 
  }

 private:

  /** @brief Links half-edge to face.*/
  void set_face(Face* f) { f_ = f; };

  /** @brief Links half-edge to vertex.*/
  void set_org(Vertex* v) { o_ = v; };

  /** @brief Updates next half-edge.*/
  void set_next(Hedge* h) { n_ = h; };

  /** @brief Links half-edge to edge.*/
  void set_edge(Edge* e) { e_ = e; };

  /** @brief Redefines the vertices of the half-edge.*/
  Hedge* reuse(Vertex *v0, Vertex *v1);
};


/**
 * @class Edge
 * @brief The Edge class represents the edge topological entity of the 4-8 mesh.
 *
 * This class inherits the Markable class properties (used to convert a general mesh into
 * a triangulated quadrangulation). The Edge class also inherits the MxHeapable class properties. 
 * The MxHeapable class is subclass of MxHeap that is a component of the priority_queue class. 
 * The MxHeap class is used to refine and simplify schemes of the mesh 4-8.
*/

class Edge : public Markable, public MxHeapable {
  friend class Mesh;

  Hedge h_[2]; ///< Two half-edges. One for each face use.

 public:

  /* @brief The construtctor initializes the vertices of the edge.*/
  Edge(Vertex *p0, Vertex *p1);

  /* @brief Returns a half-edge specified by the index i.*/
  Hedge* hedge(int i);

  /** @brief Returns the origin vertex.*/
  Vertex* org() { return h_[0].org(); }

  /** @brief Returns the end vertex.*/
  Vertex* dst() { return h_[1].org(); }

  /** @brief Returns the edge resolution level.*/
  int level() { return h_[0].level(); }

  /** @brief Returns true if the edge belongs to mesh border, else return false.
    * @return boolean type (true or false).*/
  bool is_bdry() { return (h_[0].face() == NULL || h_[1].face() == NULL); }

   /** @brief Returns true if the edge is a subdivision edge, else return false.
     *
     * The edge is a subdiv edge if its half-edges are subdiv half-edges.*/
  bool is_subdiv() { return (h_[0].is_subdiv_itself() || h_[1].is_subdiv_itself()); }
};

}

#endif
