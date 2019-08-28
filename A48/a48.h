/** 
* @file a48.h
* @brief Definition of types and classes to representing the mesh 48.
*/

/* Copyright (C) 2004 Luiz Velho. */

#ifndef A48_H
#define A48_H

#include <cstdio>
#include <iostream>
#include <set>
#include <map>
#include <iterator>
#include <QObject>
#include "heap.h"

/** 
* \defgroup A48namespace Namespace A48.
* @{
*/

/**
* @namespace A48
* @brief 
*/

namespace A48 {

class Mesh;
class Surface;
class Face;
class Edge;
class Hedge;
class Vertex;

class Markable;
class Error;


/**
* @typedef std::set<Face*>  FaceContainer
* @brief Use the class set to store the faces of mesh 4-8.
*/
typedef std::set<Face*>  FaceContainer; 


/**
* @typedef FaceContainer::iterator FaceIter
* @brief Define a iterator to access the faces of mesh 4-8.
*/
typedef FaceContainer::iterator FaceIter;


/**
* @typedef std::set<Edge*>  EdgeContainer
* @brief Use the class set to store the edge of mesh 4-8.
*/
typedef std::set<Edge*>  EdgeContainer;


/**
* @typedef EdgeContainer::iterator EdgeIter
* @brief Define a iterator to access the edge of mesh 4-8.
*/
typedef EdgeContainer::iterator EdgeIter;


/**
* @typedef std::set<Vertex*>  VertexContainer
* @brief Use the class set to store the vertices of mesh 4-8.
*/
typedef std::set<Vertex*>  VertexContainer;


/**
* @typedef VertexContainer::iterator VertexIter
* @brief Define a iterator to access the vertices of mesh 4-8.
*/
typedef VertexContainer::iterator VertexIter;


/**
* @class Ipair
* @brief This class is used to identify a half-edge of the mesh 4-8.
*/
typedef std::pair<int,int> Ipair;


/**
* @typedef std::map<const Ipair, A48::Hedge*> HedgeMap
* @brief Use the class map to improve fast access under half-edge storage.
*/
typedef std::map<const Ipair, A48::Hedge*> HedgeMap;


/**
* @class Markable
* @brief The Markable class adds a mark to an object
*
* This class is used to store the information about vertex and edge state.
* The state specifies if a vertex or a edge was processed or not processed.
*/

class Markable {
  bool     mark_; ///< Store the edge or vertex state (processed or not processed).

 public:

  /** @brief The constructor initialize the state to false.*/
  Markable(void) : mark_(false) {}

  /** @brief Update the state of variable class.
   * @param m Specify the new state of the variable class (true or false).*/
  void set_mark(bool m) { mark_ = m; };

  /** @brief Return the variable class state.
    * @return Return a boolean type (true or false).*/
  bool is_marked(void) const { return mark_; };
};


/**
* @class Error
* @brief The class Error is used to handle errors of 4-8 mesh operations.
*/
class Error {
 public:
  /** @brief The constructor allows specify a textual error.*/
  Error(char *s=""){ 
    std::cerr<<"A48 error: "<<s<<std::endl; exit(0); 
  }
};

}

#include "face.h"
#include "edge.h"
#include "vertex.h"
#include "mesh.h"

/** @} */ //end of group class.

#endif
