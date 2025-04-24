#ifndef BRIDGE_C_FORTRAN_H
#define BRIDGE_C_FORTRAN_H

#include "Polygon2D.h"


void launch_grb_();
void end_grb_();

/// @brief Generates a polygon with given vertices (x_v, y_v) composed of rectangular faces.
/// @details The polygon is composed of rectangles and has (n_x-1)*(n_y-1) faces, forming a grid of rectangles.
///          If n_x==2 and n_y==2, it produces a single rectangle, with the edges sorted in the (W, S, E, N) order, meaning
///          1. edge((x_v(1), y_v(1))->(x_v(1), y_v(2)))
///          2. edge((x_v(1), y_v(1))->(x_v(2), y_v(1)))
///          3. edge((x_v(2), y_v(1))->(x_v(2), y_v(2)))
///          4. edge((x_v(1), y_v(2))->(x_v(2), y_v(2)))
/// @param x_v 
/// @param n_x 
/// @param y_v 
/// @param n_y 
/// @return 
Polygon2D* polygon2d_from_vertices_fortran_(const double* x_v, const long int *n_x, const double* y_v, const long int *n_y);

/// @brief 
/// @param grid 
/// @param clipped 
/// @param vec_move_solid 
/// @param size_vec_move_solid 
/// @param dt 
/// @param ptr_lambdas_arr 
/// @param height_lambdas_arr 
/// @param width_lambdas_arr 
/// @param ptr_big_lambda_n 
/// @param size_big_lambda_n 
/// @param ptr_big_lambda_np1 
/// @param size_big_lambda_np1 
/// @param mean_normal 
/// @param is_narrowband 
void compute_lambdas2d_fortran_(const Polygon2D** grid, const Polygon2D **clipped, const Point2D *vec_move_solid, const int64_t *size_vec_move_solid, const double *dt, \
                        double *ptr_lambdas_arr,    \
                        double *ptr_big_lambda_n,   \
                        double *ptr_big_lambda_np1, \
                        Point3D *mean_normal, bool *is_narrowband);
GrB_Info grb_matrix_new_fp_fortran_(GrB_Matrix** M, int64_t *signed_nrows, int64_t *signed_ncols);
GrB_Info grb_matrix_setelement_(GrB_Matrix** M, int64_t *signed_irows, int64_t *signed_jcols, my_real* val);
GrB_Info grb_matrix_getelement_(GrB_Matrix** M, int64_t *signed_irows, int64_t *signed_jcols, my_real* val);
void polygon2d_free_fortran_(Polygon2D** p);
long polygon2d_get_nb_edges_(const Polygon2D** p);

#endif