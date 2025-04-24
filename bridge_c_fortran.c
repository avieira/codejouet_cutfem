#define my_real double

#include "bridge_c_fortran.h"
#include "array_double.h"
#include "vector_double.h"
#include "compute_lambdas3D.h"

void launch_grb_(){
    GrB_init(GrB_NONBLOCKING);
}

void end_grb_(){
    GrB_finalize();
}

Polygon2D* polygon2d_from_vertices_fortran_(const my_real* x_v, const long int *signed_n_x, const my_real* y_v, const long int *signed_n_y){
    const unsigned long n_x = (unsigned long) *signed_n_x;
    const unsigned long n_y = (unsigned long) *signed_n_y;
    
    Polygon2D *p = polygon2D_from_vertices(x_v, n_x, y_v, n_y);
    return p;
}

void compute_lambdas2d_fortran_(const Polygon2D** grid, const Polygon2D **clipped, const Point2D *ptr_vec_move_solid, const int64_t *signed_size_vec_move_solid, const my_real *dt, \
                        my_real *ptr_lambdas_arr,   \
                        my_real *ptr_big_lambda_n,  \
                        my_real *ptr_big_lambda_np1,\
                        Point3D *mean_normal, bool *is_narrowband)
{
    uint64_t size_vec_move_solid = (uint64_t) *signed_size_vec_move_solid;
    Vector_points2D *vec_move_solid = alloc_with_init_vec_pts2D(ptr_vec_move_solid, size_vec_move_solid);

    Array_double *lambdas = (Array_double*)malloc(sizeof(Array_double));
    Vector_double *Lambda_n = (Vector_double*)malloc(sizeof(Vector_double));
    Vector_double *Lambda_np1 = (Vector_double*)malloc(sizeof(Vector_double));

    compute_lambdas2D(*grid, *clipped, vec_move_solid, *dt, &lambdas, &Lambda_n, &Lambda_np1, mean_normal, is_narrowband);

    memcpy(ptr_lambdas_arr, lambdas->data, lambdas->ncols*lambdas->nrows*sizeof(my_real));

    memcpy(ptr_big_lambda_n, Lambda_n->data, Lambda_n->size*sizeof(my_real));
    memcpy(ptr_big_lambda_np1, Lambda_np1->data, Lambda_np1->size*sizeof(my_real));

    dealloc_vec_pts2D(vec_move_solid);
    dealloc_arr_double(lambdas);
    dealloc_vec_double(Lambda_n);
    dealloc_vec_double(Lambda_np1);
}


GrB_Info grb_matrix_new_fp_fortran_(GrB_Matrix** M, int64_t *signed_nrows, int64_t *signed_ncols){
    GrB_Index nrows = (uint64_t) *signed_nrows;
    GrB_Index ncols = (uint64_t) *signed_ncols;

    *M = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    #if my_real == double
        return GrB_Matrix_new(*M, GrB_FP64, nrows, ncols);
    #else
        return GrB_Matrix_new(*M, GrB_FP32, nrows, ncols);
    #endif
}

GrB_Info grb_matrix_setelement_(GrB_Matrix** M, int64_t *signed_irows, int64_t *signed_jcols, my_real* val){
    GrB_Index irows = (uint64_t) (*signed_irows - 1);
    GrB_Index jcols = (uint64_t) (*signed_jcols - 1);

    return GrB_Matrix_setElement(**M, *val, irows, jcols);
}

GrB_Info grb_matrix_getelement_(GrB_Matrix** M, int64_t *signed_irows, int64_t *signed_jcols, my_real* val){
    GrB_Index irows = (uint64_t) (*signed_irows - 1);
    GrB_Index jcols = (uint64_t) (*signed_jcols - 1);

    *val = 0.;
    return GrB_Matrix_extractElement(val, **M, irows, jcols); //if M(irows, jcols) is not in the matrix, val is not modified ; 
                                                              //thus it will still be 0.
}

void polygon2d_free_fortran_(Polygon2D** p){
    dealloc_Polygon2D(*p);
    *p = NULL;
}

long polygon2d_get_nb_edges_(const Polygon2D** p){
    unsigned long nb;
    GrB_Matrix_ncols(&nb, *((*p)->edges));
    return ((long) nb);
}