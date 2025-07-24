#ifndef POLYGON2D_H
#define POLYGON2D_H

#include "GraphBLAS.h"
#include "vector_Points.h"
#include "vector_int.h"
#include "vector_double.h"
#include "my_real.h"

typedef struct {
    Vector_points2D* vertices;
    GrB_Matrix* edges;
    GrB_Matrix* faces;
    Vector_int* status_edge;
    //Vector_double* pressure_edge;
} Polygon2D;

Polygon2D* new_Polygon2D();
Polygon2D* new_Polygon2D_ves(Vector_points2D* vertices, GrB_Matrix* edges, Vector_int* status_edge);
Polygon2D* new_Polygon2D_vefs(Vector_points2D* vertices, GrB_Matrix* edges, GrB_Matrix* faces, Vector_int* status_edge);
Polygon2D* polygon2D_from_vertices(const my_real* x_v, unsigned long int n_x, const my_real* y_v, unsigned long int n_y);
Polygon2D* polygon_from_consecutive_points(const my_real *x_v, const my_real* y_v, unsigned long int nb_pts);
void copy_Polygon2D(const Polygon2D *src, Polygon2D *dest);
void dealloc_Polygon2D(Polygon2D* p);
void retrieve_ith_edge2D(Vector_points2D *vertices, GrB_Matrix* edges, int i, Point2D** e_ext1, Point2D** e_ext2);
void dropzeros(GrB_Matrix* M);
Polygon2D* fuse_polygons(Polygon2D* p1, Polygon2D* p2);

#endif
