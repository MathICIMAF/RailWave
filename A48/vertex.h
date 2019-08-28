/** 
* @file vertex.h
* @brief Definition of Vertex class.
*/

/* Copyright (C) 2004 Luiz Velho. */

#ifndef VERTEX_H
#define VERTEX_H
#include <QObject>

namespace A48 {

/**
* @class Vertex
* @brief The Vertex class represents the vertex entity,
* This class stores geometry and attribute information.
*/

class Vertex : public  Markable, public MxHeapable,public QObject {
  friend class Mesh;
  friend class Face;

  int         l_;    ///< vertex resolution level.
  Hedge*      s_;    ///< vertex incoming half-edge (star handle).
  int		  idProperty_;   ///<id que tiene en la malla

 public:

	 void setId(int value){idProperty_ = value;};
	 int getId(){return idProperty_;};
  /** @brief The constructor initializates the class variables to null values.*/
         Vertex(QObject* parent=0):QObject(parent) {l_ = 0; s_ = NULL;}

  /** @brief Returns vertex resolution level.
    * @return Vertex resolution level.*/
  int  level() { return l_; };

  /** @brief Returns true if the vertex belongs to mesh boundary, else return false
    * @return boolean type (true or false).*/
  bool is_bdry() { return s_->edge()->is_bdry(); };

  /** @brief Returns true if the vertex belongs to base mesh, else return false.
    *
    * The base mesh represents the resolution level 0 of the mesh 4-8.
    * @return boolean type (true or false).*/
  bool is_inbase() { return (level() == 0); };

  /** @brief Returns true if the vertex is a weld vertex, else return false, 
    *
    * A weld vertex is a vertex that can be simplified.
    * @param exclude Face pointer.
    * @return boolean type (true or false).*/
  bool is_weld(Face *exclude = NULL) {
    int n = 0, k = 0;
    for (Hedge *e = star_first(); e != NULL; e = star_next(e), k++) {
      Face *f = e->face();
      if (f != NULL && f->weld_vertex() == this && f != exclude) n++;
    }
    return (n > 0 && k > 2);
  }

  /** @brief Returns the half-edge of the vertex.
    * @return half-edge pointer.*/
  Hedge* star_first() const { return s_; }; 

  /** @brief Returns pointer to next half-edge in face loop.
    * @param h Half-edge pointer.
    * @return Return half-edge pointer.*/
  Hedge*  star_next(Hedge* h) const {
    if (h->face() == NULL ) return NULL; // other side of boundary
    else { Hedge *n = h->next()->mate(); return (n == s_)? NULL : n; }
  }

 private:

  /** @brief Initializes the resolution level value.
    * @param l resolution level*/
  void set_level(int l) { l_ = l; }; 

  /** @brief Inicializes the vertex half-edge.
    * @param h Half-edge pointer.*/
  void set_star(Hedge* h) { s_ = h; }; 
};

}

#endif
