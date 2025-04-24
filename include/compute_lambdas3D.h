#ifndef COMPUTE_LAMBDAS3D_H
#define COMPUTE_LAMBDAS3D_H

#include "general_clipping.h"
#include "vector_double.h"
#include "array_double.h"
#include "array_Points.h"

Polyhedron3D* clip3D(const Polyhedron3D *clipper, const Polyhedron3D *clipped);
void compute_lambdas2D_time(const Polyhedron3D* grid, const Polyhedron3D *initial_p, Vector_points3D **lambdas_vec, Point3D *mean_normal, bool *is_narrowband);
void compute_lambdas2D(const Polygon2D* grid, const Polygon2D *clipped, const Vector_points2D *vec_move_solid, const double dt, \
                        Array_double **lambdas_arr, Vector_double** big_lambda_n, Vector_double** big_lambda_np1, Point3D *mean_normal, bool *is_narrowband);
Polyhedron3D* build_space2D_time_cell(const Polygon2D *fn, const Vector_points2D *vs, const double dt, bool split);

#endif