#ifndef POLYHEDRON3D_H
#define POLYHEDRON3D_H

#include "GraphBLAS.h"
#include "vector_Points.h"
#include "vector_int.h"
#include "my_real.h"

typedef struct {
    Vector_points3D* vertices;
    GrB_Matrix* edges;
    GrB_Matrix* faces;
    GrB_Matrix* volumes;
    Vector_int* status_face;
} Polyhedron3D;

Polyhedron3D* new_Polyhedron3D();
//Only does shallow copy!
Polyhedron3D* new_Polyhedron3D_vefs(Vector_points3D* vertices, GrB_Matrix* edges, GrB_Matrix* faces, Vector_int* status_face);
//Only does shallow copy!
Polyhedron3D* new_Polyhedron3D_vefvs(Vector_points3D* vertices, GrB_Matrix* edges, GrB_Matrix* faces, GrB_Matrix* volumes, Vector_int* status_face);
Polyhedron3D* Polyhedron3D_from_vertices(const my_real* x_v, unsigned long int n_x, const my_real* y_v, unsigned long int n_y, const my_real* z_v, unsigned long int n_z);
//Deep copy
void copy_Polyhedron3D(const Polyhedron3D *src, Polyhedron3D *dest);
void dealloc_Polyhedron3D(Polyhedron3D* p);

#endif
