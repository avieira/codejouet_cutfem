#define _USE_MATH_DEFINES

#include "my_real.h"
#include "bridge_c_fortran.h"
#include "array_double.h"
#include "vector_double.h"
#include "compute_lambdas3D.h"
#include "solver_for_graphblas.h"
#include "update_clipped.h"
# define M_PI           3.14159265358979323846

Polygon2D* grid = NULL;
Polygon2D* clipped = NULL;
Polyhedron3D* clipped3D = NULL;

void launch_grb_(){
    //GrB_init(GrB_NONBLOCKING);
    GrB_init(GrB_BLOCKING);
}

void end_grb_(){
    dealloc_Polygon2D(grid);
    dealloc_Polygon2D(clipped);
    dealloc_Polyhedron3D(clipped3D);

    GrB_finalize();
}

///  @brief Builds grid given consecutive points with coordinates (x_v, y_v).
///  @details If grid is already built, it only changes the coordinates of the points (since this is always the same structure of edges and faces).
///  @warning It supposes that all grid cells have the same number of points (usually 3 or 4).
void build_grid_from_points_fortran_(const my_real* x_v, const my_real* y_v, const long long int *signed_nb_pts){
    const unsigned long nb_pts = (unsigned long) *signed_nb_pts;
    Point2D p;
    unsigned long i;
    
    if (grid){
        for(i=0; i<nb_pts; i++){
            p = (Point2D){x_v[i], y_v[i]};
            set_ith_elem_vec_pts2D(grid->vertices, i, &p);
        }
    } else {
        grid = polygon_from_consecutive_points(x_v, y_v, nb_pts);
    }
}

/// @brief 
/// @param x_v 
/// @param y_v 
/// @param signed_nb_pts 
void build_clipped_from_pts_fortran_(const my_real* x_v, const my_real* y_v, const long long int *signed_nb_pts){
    const unsigned long nb_pts = (unsigned long) *signed_nb_pts;
    
    clipped = polygon_from_consecutive_points(x_v, y_v, nb_pts);
}

void compute_lambdas2d_fortran_(const my_real *dt, \
                        my_real *ptr_lambdas_arr,   \
                        my_real *ptr_big_lambda_n,  \
                        my_real *ptr_big_lambda_np1,\
                        Point3D *mean_normal, bool *is_narrowband)
{
    Array_double *lambdas = alloc_empty_arr_double();
    Vector_double *Lambda_n = alloc_empty_vec_double();
    Vector_double *Lambda_np1 = alloc_empty_vec_double();

    compute_lambdas2D(grid, clipped3D, *dt, &lambdas, &Lambda_n, &Lambda_np1, mean_normal, is_narrowband);

    memcpy(ptr_lambdas_arr, lambdas->data, lambdas->ncols*lambdas->nrows*sizeof(my_real));

    memcpy(ptr_big_lambda_n, Lambda_n->data, Lambda_n->size*sizeof(my_real));
    memcpy(ptr_big_lambda_np1, Lambda_np1->data, Lambda_np1->size*sizeof(my_real));

    dealloc_arr_double(lambdas);
    dealloc_vec_double(Lambda_n);
    dealloc_vec_double(Lambda_np1);
}



void nb_pts_clipped_fortran_(long long int* signed_nb_pts_solid){
    *signed_nb_pts_solid = ((long long) clipped->vertices->size);
}

void nb_edge_clipped_fortran_(long long int* signed_nb_edges_solid){
    unsigned long nb;
    GrB_Matrix_ncols(&nb, *(clipped->edges));
    *signed_nb_edges_solid = ((long long) nb);
}

void compute_normals_clipped_fortran_(my_real* normalVecx, my_real* normalVecy, long long* signed_nb_pts, \
                                        my_real* normalVecEdgex, my_real* normalVecEdgey, long long* signed_nb_edges, \
                                        my_real* min_pos_Se){

    unsigned long j;
    GrB_Info infogrb;
    GrB_Vector ej, extr_vals_ej;
    GrB_Vector nz_ej;
    GrB_Index size_nz_ej;
    unsigned long pt_index1, pt_index2;
    Point2D *pt1, *pt2;
    int8_t sign;
    my_real normVec;

    unsigned long nb_pts = ((unsigned long) *signed_nb_pts);
    unsigned long nb_edges = ((unsigned long) *signed_nb_edges);

    infogrb = GrB_Vector_new(&ej, GrB_INT8, nb_pts);
    infogrb = GrB_Vector_new(&nz_ej, GrB_UINT64, nb_pts);
    infogrb = GrB_Vector_new(&extr_vals_ej, GrB_INT8, nb_edges);

    *min_pos_Se = -1.0;

    for (j=0; j<nb_edges; j++){
        infogrb = GrB_extract(ej, GrB_NULL, GrB_NULL, *(clipped->edges), GrB_ALL, 1, j, GrB_NULL); 
        infogrb = GxB_Vector_extractTuples_Vector(nz_ej, extr_vals_ej, ej, GrB_NULL);
        infogrb = GrB_Vector_size(&size_nz_ej, nz_ej);

        if (size_nz_ej > 1){
            infogrb = GrB_Vector_extractElement(&pt_index1, nz_ej, 0); 
            infogrb = GrB_Vector_extractElement(&pt_index2, nz_ej, 1); 
            pt1 = get_ith_elem_vec_pts2D(clipped->vertices, pt_index1);
            pt2 = get_ith_elem_vec_pts2D(clipped->vertices, pt_index2);

            //Choice of normal: take the normal to the "edge" containing the point, sum of the adjacent edges
            GrB_Matrix_extractElement(&sign, *(clipped->edges), pt_index1, j);
            normalVecx[pt_index1] -= sign * pt1->y;
            normalVecy[pt_index1] += sign * pt1->x;
            normalVecx[pt_index2] -= sign * pt1->y;
            normalVecy[pt_index2] += sign * pt1->x;
            normalVecEdgex[j] -= sign * pt1->y;
            normalVecEdgey[j] += sign * pt1->x;

            GrB_Matrix_extractElement(&sign, *(clipped->edges), pt_index2, j);
            normalVecx[pt_index1] -= sign * pt2->y;
            normalVecy[pt_index1] += sign * pt2->x;
            normalVecx[pt_index2] -= sign * pt2->y;
            normalVecy[pt_index2] += sign * pt2->x;
            normalVecEdgex[j] -= sign * pt2->y;
            normalVecEdgey[j] += sign * pt2->x;

            normVec = sqrt(normalVecEdgex[j]*normalVecEdgex[j] + normalVecEdgey[j]*normalVecEdgey[j]);
            if (normVec > 0.){
                normalVecEdgex[j] /= normVec;
                normalVecEdgey[j] /= normVec;

                if (*min_pos_Se > 0.){
                    if (normVec < *min_pos_Se)
                        *min_pos_Se = normVec;
                } else {
                    *min_pos_Se = normVec;
                }
            }
        }
    }
    GrB_Vector_free(&ej);
    GrB_Vector_free(&extr_vals_ej);
    GrB_Vector_free(&nz_ej);
}

void smooth_vel_clipped_fortran_(my_real* vec_move_clippedx, my_real* vec_move_clippedy, my_real* min_pos_Se, my_real *dt){
    my_real eps;
    GrB_Matrix smoothing_op, id;
    unsigned long nb_pts = clipped->vertices->size;
    GrB_Index ind_vector[2] = {1, nb_pts};
    
    eps = 0.001 / ((*min_pos_Se)*(*min_pos_Se));
    GrB_Matrix_new(&smoothing_op, GrB_FP64, nb_pts, nb_pts);
    GrB_Matrix_new(&id, GrB_FP64, nb_pts, nb_pts);

    //smoothing_op = Id + eps*dt*(edges*edges')
    GrB_mxm(smoothing_op, GrB_NULL, GrB_NULL, GrB_NULL, *(clipped->edges), *(clipped->edges), GrB_DESC_T1);
    GrB_assign(id, GrB_NULL, GrB_NULL, 1.0, ind_vector, GxB_RANGE, ind_vector, GxB_RANGE, GrB_NULL);
    GrB_apply(smoothing_op, GrB_NULL, GrB_NULL, GrB_TIMES_FP64, eps*(*dt), smoothing_op, GrB_NULL);
    GrB_eWiseAdd(smoothing_op, GrB_NULL, GrB_NULL, GrB_PLUS_FP64, id, smoothing_op, GrB_NULL);

    solver_for_graphblas_FP_sym(smoothing_op, vec_move_clippedx);
    solver_for_graphblas_FP_sym(smoothing_op, vec_move_clippedy);

    GrB_free(&smoothing_op);
    GrB_free(&id);
}

void get_clipped_ith_vertex_fortran_(long long *k, Point2D *pt, long long *signed_eR, long long *signed_eL){
    uint64_t eR, eL;
    GrB_Matrix e_k;
    GrB_Vector nz_e_k, extr_vals_e_k, I_vec_e_k;
    GrB_Index size_nz_e_k;
    GrB_Index nb_pts;
    GrB_Info infogrb;

    nb_pts = clipped->vertices->size;
    *pt = *get_ith_elem_vec_pts2D(clipped->vertices, *k); 

    infogrb = GrB_Matrix_new(&e_k, GrB_INT8, 1, nb_pts);
    infogrb = GrB_Vector_new(&nz_e_k, GrB_UINT64, 2);
    infogrb = GrB_Vector_new(&extr_vals_e_k, GrB_INT8, 2);
    infogrb = GrB_Vector_new(&I_vec_e_k, GrB_INT8, 2);
    GrB_extract(e_k, GrB_NULL, GrB_NULL, *(clipped->edges), k, 1, GrB_ALL, 1, GrB_NULL); //Get indices of edges connected with point k
    GxB_Matrix_extractTuples_Vector(I_vec_e_k, nz_e_k, extr_vals_e_k, e_k, GrB_NULL);
    GrB_Vector_size(&size_nz_e_k, nz_e_k);
    if (size_nz_e_k == 2){ //Exactly two edges connected to the point k
        GrB_Vector_extractElement(&eR, nz_e_k, 0);
        GrB_Vector_extractElement(&eL, nz_e_k, 1);
        *signed_eR = (long long)eR;
        *signed_eL = (long long)eL;
    } else {                   
        *signed_eR = -1;       
        *signed_eL = -1;       
    }                          
    GrB_free(&e_k);            
    GrB_free(&nz_e_k);         
    GrB_free(&extr_vals_e_k);  
    GrB_free(&I_vec_e_k);      
}                              
                               
void update_clipped_fortran_(const my_real* vec_move_clippedy, const my_real* vec_move_clippedz, const my_real* dt, const my_real *dx){
    const my_real minimal_length = *dx * 0.1;
    const my_real maximal_length = *dx * 2;
    const my_real minimal_angle = 0.1*M_PI;

    update_solid(&clipped, &clipped3D, vec_move_clippedy, vec_move_clippedz, *dt, \
                    minimal_length, maximal_length, minimal_angle);
}                              
                               
                               
                               
                               
                               