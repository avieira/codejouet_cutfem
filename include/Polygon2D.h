#ifndef POLYGON2D_H
#define POLYGON2D_H

#include "GraphBLAS.h"
#include "vector_Points.h"
#include "vector_int.h"

typedef struct {
    Vector_points2D* vertices;
    GrB_Matrix* edges;
    GrB_Matrix* faces;
    Vector_int* status_edge;
} Polygon2D;

Polygon2D* new_Polygon2D();
Polygon2D* new_Polygon2D_ves(Vector_points2D* vertices, GrB_Matrix* edges, Vector_int* status_edge);
Polygon2D* new_Polygon2D_vefs(Vector_points2D* vertices, GrB_Matrix* edges, GrB_Matrix* faces, Vector_int* status_edge);
Polygon2D* polygon2D_from_vertices(const double* x_v, unsigned long int n_x, const double* y_v, unsigned long int n_y);
void copy_Polygon2D(const Polygon2D *src, Polygon2D *dest);
void dealloc_Polygon2D(Polygon2D* p);

#endif
