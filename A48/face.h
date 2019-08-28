/** 
* @file face.h
* @brief Definition of the face topological entity of the mesh 4-8.
*/

/* Copyright (C) 2004 Luiz Velho. */

#ifndef FACE_H
#define FACE_H

#include <QObject>

namespace A48 {

/**
* @class Face
* @brief The class Face implements the face topological entity of the mesh 4-8.
*/
  class Face:public QObject {
  friend class Mesh;

  Hedge  *e_; ///< Points to the first half-egde of face's loop.
  int saved_;

 public:

	 int getSaved(){return saved_;}
	 void setSaved(int value){saved_ = value;}
 /** @brief The constructor creates a triangular face.
   * @param e0 First half-edge.
   * @param e1 Second half-edge.
   * @param e2 Third half-edge.*/
  Face(Hedge* e0, Hedge* e1, Hedge* e2,QObject* parent=0):QObject(parent) { reuse(e0, e1, e2); }

  /** @brief Returns the half-edge of index k.
   * @return half-edge pointer.*/
  Hedge* hedge(int k);

 /** @brief Returns the vertex of index k.
   * @return vertex pointer.*/
  Vertex* vertex(int k);

  /** @brief Returns face level.*/
  int level();

  /** @brief Returns true if the face belongs to the level 0, else return false.*/
  bool is_inbase();

  /** @brief Returns the face's subdivision edge.
   * @return half-edge pointer.*/
  Hedge* subd_edge() { return hedge(0); };

  /** @brief Returns the face vertex that can be used for simplification.*/
  Vertex* weld_vertex() { return vertex(0); };

 private:
	 
  /** @brief Updates the half-edge of index k.*/
  void set_hedge(int k, Hedge* h);

  /** @brief Updates the vertex of index k.*/
  void set_vertex(int k, Vertex* v);

  /** @brief Links the vertices to the half-edges.*/
  void link_star_verts();

  /** @brief Redefines a face by specifing replacement half-edges.*/
  Face* reuse(Hedge *e0, Hedge *e1, Hedge *e2);
};

}

#endif
