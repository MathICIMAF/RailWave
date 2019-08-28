/** 
* @file mesh.h
* @brief Definition of Surface and mesh classes.
*/

/* Copyright (C) 2004 Luiz Velho. */

#ifndef MESH_H
#define MESH_H

#include <QObject>

namespace A48 {

/** 
* @class Surface
* @brief The the abstract class Surface defines a interface to multi-resolution surface.
*
* The Surface class provides information to the adative mesh.
* The interface describes a set of operations to build the base mesh, 
* control the surface resolution and provide geometric information.
*/
 class Surface {

 public:
     Surface(){}
  /** @brief Specifies base mesh geometry and topology (the coarsest resolution level)
    * @param np Pointer to number of vertices.
    * @param fcs Pointer to array of indices that describes the mesh connectivity.
    * @param nf Pointer to number of faces.*/
  virtual void base_mesh(int *np, int **fcs, int *nf) = 0;

  /** @brief Specifies vertex attributes (geometry, color, normal, prioriry) of index i.
    * @param i Index of the vertex.
    * @param v Vertex pointer.*/
  virtual void sample(int i, Vertex* v) = 0;

  /** @brief Receives edge pointer and calculates vertex coordinate.
    * @param e Edge pointer.
    * @param v Vertex pointer.*/
  virtual void sample(Edge* e, Vertex* v) = 0;

  /** @brief Receives a face pointer and calculates vertex coordinate.
    * @param f Face pointer.
    * @param v Vertex pointer.*/
  virtual void sample(Face* f, Vertex* v) = 0;

  /** @brief Returns the edge length.
    * @param e Edge pointer.
    * @return Return the edge length.*/
  virtual float elenght(Edge* e) = 0;

  /** @brief Returns the rank for split edge operation.
    * @param e Edge pointer.
    * @return Return the split priority.*/
  virtual float ref_rank(Edge* e) = 0;

  /** @brief Return the rank for weld vertex operation.
    * @param v Vertex pointer.
    * @return Return the the weld priority.*/
  virtual float simpl_rank(Vertex* v) = 0;

  /** @brief Creates new vertex.
    * @return Vertex pointer.*/
  virtual Vertex* new_vertex(void) = 0;

  /** @brief Removes vertex.
    * @param v Vertex pointer.*/
  virtual void del_vertex(Vertex *v) = 0;
};


/**
* @class Mesh
* @brief The class Mesh defines data structures and methods that represent a 4-8 mesh.
*/

 class Mesh:public QObject {

     Surface      *surf_; ///< Surface pointer
     FaceContainer   fc_; ///< Set of faces
     EdgeContainer   ec_; ///< Set of edges
     VertexContainer vc_; ///< Set of vertices
     MxHeap          sf_; ///< Priority queue for simplification front
     MxHeap          rf_; ///< Priority queue for refinement front

 public:
  /** @brief The construtor specifies a surface class.*/
     Mesh(Surface *s,QObject* parent=0);

  /** @brief The mesh class destructor. Destroys the mesh data structures.*/
  ~Mesh();

  /** @brief Returns the surface pointer.
    * @return Surface pointer.*/
  Surface* surf() { return surf_; };
 
  /** @brief Returns iterator that points to the first vertex of the mesh.
    * @return Vertex iterator.*/
  VertexIter verts_begin(){ return vc_.begin(); };

  /** @brief Returns  iterador that points to the last vertex of the mesh.
    * @return Vertex iterator.*/
  VertexIter verts_end(){ return vc_.end(); };

  /** @brief Returns the number of vertices.
    * @return Number of faces.*/
  int num_verts() { return vc_.size(); };

  /** @brief Returns iterator that points to the first edge of the mesh.
    * @return Edge iterator.*/
  EdgeIter edges_begin() { return ec_.begin(); };

  /** @brief Returns iterador that points to the last vertex of the mesh.
    * @return Edge iterator.*/
  EdgeIter edges_end() { return ec_.end(); };

  /** @brief Returns the number of egdes.
    * @return Number of edges.*/
  int num_edges() { return ec_.size(); };

  /** @brief Returns iterator that points to the first face of the mesh.
    * @return Face iterator.*/
  FaceIter faces_begin() { return fc_.begin(); };

  /** @brief Returns iterador that points to the last vertex of the mesh.
    * @return Face iterator.*/
  FaceIter faces_end() { return fc_.end(); };

  /** @brief Returns the number of faces.
    * @return Number of faces.*/
  int num_faces() { return fc_.size(); };

  /** @brief Refines the mesh based on threshold t.*/
  void adapt_refine(double t);

  /** @brief Simplifies the mesh based on threshold t.*/
  void adapt_simplify(double t);

  /** @brief Splits an edge.*/
  Vertex* refine(Hedge* e);

  /** @brief Removes a vertex from mesh.*/
  Hedge* simplify(Vertex* w);

  /** @brief Divides a face.*/
  Vertex* split(Face* f);

  /** @brief Splits an edge.*/
  Vertex* split(Hedge* e);

  /** @brief Removes a vertex.*/
  Hedge* weld(Vertex* w);

  /** @brief Flips internal edge.*/
  Hedge* flip(Hedge *h);   

 private:
	
  /** @brief Adds a new vertex.*/
  Vertex *add_vertex(void);

  /** @brief Adds an existing vertex.
    * @param v Vertex pointer.
     * @return Return true if add operation was completed, else return false.*/
  bool add_vertex(Vertex *v) 
    { std::pair<VertexIter, bool> r = vc_.insert(v); return r.second; };

  /** @brief Removes a vertex from mesh.
    * @param v Vertex pointer.*/
  void del_vertex(Vertex* v) { vc_.erase(v); surf()->del_vertex(v); };

  /** @brief Adds a new edge.*/
  Hedge *add_edge(Vertex *v0, Vertex *v1);

  /** @brief Add an existing edge.
    * @param e Edge pointer. 
    * @return true if add operation was completed, else return false.*/
  bool add_edge(Edge *e) 
    { std::pair<EdgeIter, bool> r = ec_.insert(e);  return r.second;};

  /** @brief Removes an edge from mesh.
    * @param e Edge pointer.*/
  void del_edge(Edge* e) { ec_.erase(e); delete(e); }

  /** @brief Adds a new face.*/
  Face* add_face(Hedge* e0, Hedge* e1, Hedge* e2);

  /** @brief Adds an existing face.
    * @param f Face pointer.
    * @return true if add operation was completed, else return false.*/
  bool add_face(Face* f) 
    { std::pair<FaceIter, bool> r = fc_.insert(f); return r.second;}

  /** @brief Removes a face from mesh.*/
  void del_face(Face* f) { fc_.erase(f); delete(f); }

  /** @brief Creates links among faces, egdes and verticess.*/
  void link_mesh();

  /** @brief Splits an edge in two sub-edges.*/
  Vertex* bisect(Hedge* e, Hedge** el, Hedge** er);

  /** @brief Splits a face in two sub-faces.*/
  Hedge* bisect(Face *f, Hedge *e1, Hedge *e2, Hedge *el, Hedge *er);

  /** @brief Creates the simplifying and refining queues.*/
  void new_front();

  /** @brief Inserts an edge in the refining priority queue.*/
  void update_ref_front(Edge* e);

  /** @brief Inserts a vertex in the refining priority queue.*/
  void update_ref_front(Vertex* v);

  /** @brief Inserts a vertex in the simplifying priority queue.*/
  void update_simpl_front(Vertex* v);

  /** @brief Inserts an edge in the simplifying priority queue.*/
  void update_simpl_front(Edge* e); 

  /** @brief Returns true if the mesh is a triangulated-quadrangulation, else return false.*/
  bool is_triquad(void);

  /** @brief Converts the mesh to a triagulated-quadrangulation.*/
  void make_triquad(void);

  /** @brief Specifies a base mesh.*/
  void set_base_mesh(int np, int* tris, int nf);

  /** @brief Returns half-edge pointer identified by the vertices of index i0 and i1.*/
  Hedge* get_hedge(int i0, int i1, Vertex* verts[], HedgeMap* edges);

  /** @brief Creates a face from half-edges. The half-edges are specified by vertices of index i0, i1 and i2.*/
  void put_face(int i0, int i1, int i2, Vertex* verts[], HedgeMap* edges);
};


}

#endif
