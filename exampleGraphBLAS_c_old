#include "GraphBLAS.h"

#include "Point.h"
#include "vector_Points.h"
#include "Polygon2D.h"
#include "Polyhedron3D.h"
#include "compute_lambdas3D.h"
#include <stdio.h>

void test1(){
    int i;
    const int nb_points = 7;
    Point2D *p1, *p2, *p3, *p4;
    Point2D points[] = {{0.0, 0.0},
                        {1.0, 0.0},
                        {1.0, 1.0},
                        {0.0, 1.0},
                        {2.0, 0.0},
                        {2.0, 1.0},
                        {1.0, 2.0}};
    

    Vector_points2D *v1, *v2, *v3, *v4;

    v1 = alloc_empty_vec_pts2D();
    for (i=0; i<nb_points; i++){
        push_back_vec_pts2D(v1, points + i);
        printf("Capacity of v1 = %lu\n", v1->capacity);
    }

    v2 = alloc_with_capacity_vec_pts2D(nb_points);
    for (i=0; i<nb_points; i++){
        push_back_vec_pts2D(v2, points + i);
        printf("Capacity of v2 = %lu\n", v2->capacity);
    }
    
    v3 = alloc_with_init_vec_pts2D(points, nb_points);
    printf("Capacity of v3 = %lu\n", v3->capacity);

    v4 = alloc_empty_vec_pts2D();
    copy_vec_pts2D(v1, v4);

    for (i=0; i<nb_points; i++){
        p1 = get_ith_elem_vec_pts2D(v1, i);
        p2 = get_ith_elem_vec_pts2D(v2, i);
        p3 = get_ith_elem_vec_pts2D(v3, i);
        p4 = get_ith_elem_vec_pts2D(v4, i);
       printf("point #%d : v1 = (%lf, %lf), v2 = (%lf, %lf), v3 = (%lf, %lf), v4=(%lf, %lf)\n", i, p1->x, p1->y, p2->x, p2->y, p3->x, p3->y, p4->x, p4->y); 
    }

    printf("Scalar product between p1(%lf, %lf) and p2(%lf, %lf) = %lf\n", points[1].x, points[1].y, points[2].x, points[2].y, scalProd2D(points[1], points[2]));
    
    *p1 = points[1];
    *p2 = points[2];
    inplace_axpy_points2D(p1, 2.0, p2);
    printf("%f * p1(%lf, %lf) + p2(%lf, %lf) = (%lf, %lf)\n", 2.0, points[1].x, points[1].y, points[2].x, points[2].y, p1->x, p1->y);

    *p1 = points[1];
    *p2 = points[2];
    inplace_xpay_points2D(p1, 2.0, p2);
    printf("p1(%lf, %lf) + %f * p2(%lf, %lf) = (%lf, %lf)\n", points[1].x, points[1].y, 2.0, points[2].x, points[2].y, p1->x, p1->y);


    dealloc_vec_pts2D(v1);
    dealloc_vec_pts2D(v2);
    dealloc_vec_pts2D(v3);
}

void test2(){
    Polygon2D *p1, *p2, *p3, *p4, *p5;
    Point2D *pt1, *pt2, *pt3, *pt4, *pt5;
    long int *se1, *se2, *se3, *se4, *se5;
    unsigned int i;
    const unsigned int n_x = 3;
    const double x_v[] = {0.0,0.5,1.0};
    const unsigned int n_y = 3;
    const double y_v[] = {0.0,0.5,1.0};

    p1 = polygon2D_from_vertices(x_v, n_x, y_v, n_y);
    
    p2 = new_Polygon2D();
    copy_vec_pts2D(p1->vertices, p2->vertices);
    GrB_Matrix_dup(p2->edges, *(p1->edges));
    GrB_Matrix_dup(p2->faces, *(p1->faces));
    copy_vec_int(p1->status_edge, p2->status_edge);

    p3 = new_Polygon2D_ves(p1->vertices, p1->edges, p1->status_edge);

    p4 = new_Polygon2D_vefs(p1->vertices, p1->edges, p1->faces, p1->status_edge);

    p5 = new_Polygon2D();
    copy_Polygon2D(p1, p5);

    for (i=0; i<p1->vertices->size; i++){
        pt1 = get_ith_elem_vec_pts2D(p1->vertices, i);
        pt2 = get_ith_elem_vec_pts2D(p2->vertices, i);
        pt3 = get_ith_elem_vec_pts2D(p3->vertices, i);
        pt4 = get_ith_elem_vec_pts2D(p4->vertices, i);
        pt5 = get_ith_elem_vec_pts2D(p4->vertices, i);
       printf("point #%d : v1 = (%1.1lf, %1.1lf), v2 = (%1.1lf, %1.1lf), v3 = (%1.1lf, %1.1lf), v4=(%1.1lf, %1.1lf), v5=(%1.1lf, %1.1lf)\n", i, \
                pt1->x, pt1->y, pt2->x, pt2->y, pt3->x, pt3->y, pt4->x, pt4->y, pt5->x, pt5->y); 
    }

    GxB_print(*(p1->edges), GxB_COMPLETE);
    GxB_print(*(p2->edges), GxB_COMPLETE);
    GxB_print(*(p3->edges), GxB_COMPLETE);
    GxB_print(*(p4->edges), GxB_COMPLETE);
    GxB_print(*(p5->edges), GxB_COMPLETE);

    GxB_print(*(p1->faces), GxB_COMPLETE);
    GxB_print(*(p2->faces), GxB_COMPLETE);
    GxB_print(*(p3->faces), GxB_COMPLETE);
    GxB_print(*(p4->faces), GxB_COMPLETE);
    GxB_print(*(p5->faces), GxB_COMPLETE);

    for (i=0; i<p1->status_edge->size; i++){
        se1 = get_ith_elem_vec_int(p1->status_edge, i);
        se2 = get_ith_elem_vec_int(p2->status_edge, i);
        se3 = get_ith_elem_vec_int(p3->status_edge, i);
        se4 = get_ith_elem_vec_int(p4->status_edge, i);
        se5 = get_ith_elem_vec_int(p4->status_edge, i);
       printf("status_edge #%d : se1 = %ld, se2 = %ld, v3 = %ld, v4= %ld, v5= %ld\n", i, \
                *se1, *se2, *se3, *se4, *se5); 
    }
}

void test3(){
    Polyhedron3D *p1, *p2, *p3, *p4, *p5;
    Point3D *pt1, *pt2, *pt3, *pt4, *pt5;
    long int *se1, *se2, *se3, *se4, *se5;
    unsigned int i;
    const unsigned int n_x = 3;
    const double x_v[] = {0.0,0.5,1.0};
    const unsigned int n_y = 3;
    const double y_v[] = {0.0,0.5,1.0};
    const unsigned int n_z = 3;
    const double z_v[] = {0.0,0.5,1.0};

    p1 = Polyhedron3D_from_vertices(x_v, n_x, y_v, n_y, z_v, n_z);
    
    p2 = new_Polyhedron3D();
    copy_vec_pts3D(p1->vertices, p2->vertices);
    GrB_Matrix_dup(p2->edges, *(p1->edges));
    GrB_Matrix_dup(p2->faces, *(p1->faces));
    GrB_Matrix_dup(p2->volumes, *(p1->volumes));
    copy_vec_int(p1->status_face, p2->status_face);

    p3 = new_Polyhedron3D_vefs(p1->vertices, p1->edges, p1->faces, p1->status_face);

    p4 = new_Polyhedron3D_vefvs(p1->vertices, p1->edges, p1->faces, p1->volumes, p1->status_face);

    p5 = new_Polyhedron3D();
    copy_Polyhedron3D(p1, p5);

    for (i=0; i<p1->vertices->size; i++){
        pt1 = get_ith_elem_vec_pts3D(p1->vertices, i);
        pt2 = get_ith_elem_vec_pts3D(p2->vertices, i);
        pt3 = get_ith_elem_vec_pts3D(p3->vertices, i);
        pt4 = get_ith_elem_vec_pts3D(p4->vertices, i);
        pt5 = get_ith_elem_vec_pts3D(p4->vertices, i);
       printf("point #%d : v1 = (%1.1lf, %1.1lf, %1.1lf), v2 = (%1.1lf, %1.1lf, %1.1lf), \
        v3 = (%1.1lf, %1.1lf, %1.1lf), v4=(%1.1lf, %1.1lf, %1.1lf), v5=(%1.1lf, %1.1lf, %1.1lf)\n", i, \
                pt1->x, pt1->y, pt1->t, pt2->x, pt2->y, pt2->t, pt3->x, pt3->y, pt3->t, pt4->x, pt4->y, pt4->t, pt5->x, pt5->y, pt5->t); 
    }

    GxB_print(*(p1->edges), GxB_COMPLETE);
    GxB_print(*(p2->edges), GxB_COMPLETE);
    GxB_print(*(p3->edges), GxB_COMPLETE);
    GxB_print(*(p4->edges), GxB_COMPLETE);
    GxB_print(*(p5->edges), GxB_COMPLETE);

    GxB_print(*(p1->faces), GxB_COMPLETE);
    GxB_print(*(p2->faces), GxB_COMPLETE);
    GxB_print(*(p3->faces), GxB_COMPLETE);
    GxB_print(*(p4->faces), GxB_COMPLETE);
    GxB_print(*(p5->faces), GxB_COMPLETE);

    GxB_print(*(p1->volumes), GxB_COMPLETE);
    GxB_print(*(p2->volumes), GxB_COMPLETE);
    GxB_print(*(p3->volumes), GxB_COMPLETE);
    GxB_print(*(p4->volumes), GxB_COMPLETE);
    GxB_print(*(p5->volumes), GxB_COMPLETE);

    for (i=0; i<p1->status_face->size; i++){
        se1 = get_ith_elem_vec_int(p1->status_face, i);
        se2 = get_ith_elem_vec_int(p2->status_face, i);
        se3 = get_ith_elem_vec_int(p3->status_face, i);
        se4 = get_ith_elem_vec_int(p4->status_face, i);
        se5 = get_ith_elem_vec_int(p4->status_face, i);
       printf("status_face #%d : se1 = %ld, se2 = %ld, v3 = %ld, v4= %ld, v5= %ld\n", i, \
                *se1, *se2, *se3, *se4, *se5); 
    }
}

void test4(){
    Polyhedron3D *clipper, *clipped;
    Polyhedron3D *result;
    Point3D *pt;
    long int *se = (long int*)malloc(sizeof(long int));
    unsigned int i;
    const unsigned int n_x = 2;
    const double x_v[] = {0.0,1.0};
    const double x_v_clipped[] = {0.5,1.5};
    const unsigned int n_y = 2;
    const double y_v[] = {0.0,1.0};
    const double y_v_clipped[] = {0.5,1.5};
    const unsigned int n_z = 2;
    const double z_v[] = {0.0,1.0};
    const double z_v_clipped[] = {0.5,1.5};

    clipper = Polyhedron3D_from_vertices(x_v, n_x, y_v, n_y, z_v, n_z);
    for (i=0; i<clipper->status_face->size; i++){
        *se = 3;
        set_ith_elem_vec_int(clipper->status_face, i, se);
    }

    clipped = Polyhedron3D_from_vertices(x_v_clipped, n_x, y_v_clipped, n_y, z_v_clipped, n_z);
    printf("Before clipping\n");
    for (i=0; i<clipped->vertices->size; i++){
        pt = get_ith_elem_vec_pts3D(clipped->vertices, i);
       printf("point #%2d : v = (%1.3lf, %1.3lf, %1.3lf)\n", i, pt->x, pt->y, pt->t); 
    }

    result = clip3D(clipper, clipped);

    printf("After clipping\n");
    for (i=0; i<result->vertices->size; i++){
        pt = get_ith_elem_vec_pts3D(result->vertices, i);
       printf("point #%2d : v = (%1.3lf, %1.3lf, %1.3lf)\n", i, pt->x, pt->y, pt->t); 
    }
    GxB_print(*(result->edges), GxB_COMPLETE);
    GxB_print(*(result->faces), GxB_COMPLETE);
    GxB_print(*(result->volumes), GxB_COMPLETE);
    for (i=0; i<result->status_face->size; i++){
        se = get_ith_elem_vec_int(result->status_face, i);
        printf("status_face #%2d : se = %ld\n", i, *se); 
    }
}

void test4bis(){
    Polyhedron3D *clipper, *clipped;
    Polyhedron3D *result;
    Point3D *pt;
    long int *se = (long int*)malloc(sizeof(long int));
    unsigned int i;
    const unsigned int n_x = 2;
    const double x_v[] = {0.0,1.0};
    const double x_v_clipped[] = {0.5,1.5};
    const unsigned int n_y = 2;
    const double y_v[] = {0.0,1.0};
    const double y_v_clipped[] = {0.5,1.5};
    const unsigned int n_z = 2;
    const double z_v[] = {0.0,1.0};
    const double z_v_clipped[] = {0.0,1.0};

    clipper = Polyhedron3D_from_vertices(x_v, n_x, y_v, n_y, z_v, n_z);
    for (i=0; i<clipper->status_face->size; i++){
        *se = 3;
        set_ith_elem_vec_int(clipper->status_face, i, se);
    }

    clipped = Polyhedron3D_from_vertices(x_v_clipped, n_x, y_v_clipped, n_y, z_v_clipped, n_z);
    *se = 1;
    set_ith_elem_vec_int(clipped->status_face, 0, se);
    *se = 2;
    set_ith_elem_vec_int(clipped->status_face, 1, se);
    for (i=2; i<clipped->status_face->size; i++){
        *se = -1;
        set_ith_elem_vec_int(clipped->status_face, i, se);
    }


    printf("Before clipping\n");
    for (i=0; i<clipped->vertices->size; i++){
        pt = get_ith_elem_vec_pts3D(clipped->vertices, i);
       printf("point #%2d : v = (%1.3lf, %1.3lf, %1.3lf)\n", i, pt->x, pt->y, pt->t); 
    }

    result = clip3D(clipper, clipped);

    printf("After clipping\n");
    for (i=0; i<result->vertices->size; i++){
        pt = get_ith_elem_vec_pts3D(result->vertices, i);
       printf("point #%2d : v = (%1.3lf, %1.3lf, %1.3lf)\n", i, pt->x, pt->y, pt->t); 
    }
    GxB_print(*(result->edges), GxB_COMPLETE);
    GxB_print(*(result->faces), GxB_COMPLETE);
    GxB_print(*(result->volumes), GxB_COMPLETE);
    for (i=0; i<result->status_face->size; i++){
        se = get_ith_elem_vec_int(result->status_face, i);
        printf("status_face #%2d : se = %ld\n", i, *se); 
    }
}

void test5(){
    Polygon2D *clipper2D, *clipped2D;
    Vector_points2D *vs;
    Polyhedron3D *clipper, *clipped;
    Point2D pt = (Point2D){0.,0.};
    long unsigned i;
    const unsigned int n_x = 2;
    const double x_v[] = {0.0,1.0};
    const double x_v_clipped[] = {0.5,1.5};
    const unsigned int n_y = 2;
    const double y_v[] = {0.0,1.0};
    const double y_v_clipped[] = {0.5,1.5};

    clipper2D = polygon2D_from_vertices(x_v, n_x, y_v, n_y);
    clipped2D = polygon2D_from_vertices(x_v_clipped, n_x, y_v_clipped, n_y);
    vs = alloc_with_capacity_vec_pts2D(clipper2D->vertices->size);
    for(i=0; i<clipper2D->vertices->size; i++){
        set_ith_elem_vec_pts2D(vs, i, &pt);
    }

    clipper = build_space2D_time_cell(clipper2D, vs, 1.0, false);
    clipped = build_space2D_time_cell(clipped2D, vs, 1.0, true);

    print_vec_pt3D(*clipper->vertices);
    GrB_wait (*clipper->edges, GrB_MATERIALIZE);
    GxB_print(*clipper->edges, GxB_COMPLETE);
    GrB_wait (*clipper->faces, GrB_MATERIALIZE);
    GxB_print(*clipper->faces, GxB_COMPLETE);
    GrB_wait (*clipper->volumes, GrB_MATERIALIZE);
    GxB_print(*clipper->volumes, GxB_COMPLETE);
    print_vec_int(clipper->status_face);

    print_vec_pt3D(*clipped->vertices);
    GrB_wait (*clipped->edges, GrB_MATERIALIZE);
    GxB_print(*clipped->edges, GxB_COMPLETE);
    GrB_wait (*clipped->faces, GrB_MATERIALIZE);
    GxB_print(*clipped->faces, GxB_COMPLETE);
    GrB_wait (*clipped->volumes, GrB_MATERIALIZE);
    GxB_print(*clipped->volumes, GxB_COMPLETE);
    print_vec_int(clipped->status_face);
}

void test6(){
    Polygon2D *clipper2D, *clipped2D;
    Vector_points2D *vs;
    Polyhedron3D *clipper, *clipped, *clipped_in;
    Point2D pt = (Point2D){0.,0.};
    Vector_points3D *lambdas;
    Point3D normal;
    bool is_narrowband;
    long unsigned i;
    long *sfi;
    const unsigned int n_x = 2;
    const double x_v[] = {0.0,1.0};
    const double x_v_clipped[] = {0.5,1.5};
    const unsigned int n_y = 2;
    const double y_v[] = {0.0,1.0};
    const double y_v_clipped[] = {0.5,1.5};

    clipper2D = polygon2D_from_vertices(x_v, n_x, y_v, n_y);
    clipped2D = polygon2D_from_vertices(x_v_clipped, n_x, y_v_clipped, n_y);
    vs = alloc_with_capacity_vec_pts2D(clipper2D->vertices->size);
    for(i=0; i<clipper2D->vertices->size; i++){
        set_ith_elem_vec_pts2D(vs, i, &pt);
    }

    clipper = build_space2D_time_cell(clipper2D, vs, 1.0, false);
    clipped = build_space2D_time_cell(clipped2D, vs, 1.0, true);

    for (i=0; i<clipped->status_face->size; i++){
        sfi = get_ith_elem_vec_int(clipped->status_face, i);
        if (*sfi > 2){
            *sfi = -1;
        }
    }

    clipped_in = clip3D(clipper, clipped);

    //printf("clipped_in->vertices = [\n");
    //for(i=0; i<clipped_in->vertices->size; i++){
    //    pt3D = get_ith_elem_vec_pts3D(clipped_in->vertices, i);
    //    printf("pt #%ld ", i);
    //    print_pt3D(pt3D);
    //    printf("\n");
    //}
    //printf("]\n");
    //GrB_wait (*clipped_in->edges, GrB_MATERIALIZE);
    //GxB_print(*clipped_in->edges, GxB_COMPLETE);
    //GrB_wait (*clipped_in->faces, GrB_MATERIALIZE);
    //GxB_print(*clipped_in->faces, GxB_COMPLETE);
    //GrB_wait (*clipped_in->volumes, GrB_MATERIALIZE);
    //GxB_print(*clipped_in->volumes, GxB_COMPLETE);
    //printf("clipped_in->status_face = ");
    //print_vec_int(clipped_in->status_face);

    compute_lambdas2D_time(clipper, clipped, &lambdas, &normal, &is_narrowband);
    printf("lambdas = ");
    print_vec_pt3D(*lambdas);
    printf("\n");
    printf("normal = ");
    print_pt3D(&normal);
    printf("\n");
    printf("is_narrowband = %d\n", is_narrowband);
}

/*
void test7_(){
    Polygon2D *clipper2D, *clipped2D;
    Vector_points2D *vs;
    Point2D pt = (Point2D){0.,0.};
    Array_double *lambdas;
    Vector_double *Lambda_n, *Lambda_np1;
    Point3D normal;
    bool is_narrowband;
    long unsigned i;
    const unsigned int n_x = 2;
    const double x_v[] = {0.0,1.0};
    const double x_v_clipped[] = {0.5,1.5};
    const unsigned int n_y = 2;
    const double y_v[] = {0.0,1.0};
    const double y_v_clipped[] = {0.5,1.5};

    clipper2D = polygon2D_from_vertices(x_v, n_x, y_v, n_y);
    clipped2D = polygon2D_from_vertices(x_v_clipped, n_x, y_v_clipped, n_y);
    vs = alloc_with_capacity_vec_pts2D(clipper2D->vertices->size);
    for(i=0; i<clipper2D->vertices->size; i++){
        set_ith_elem_vec_pts2D(vs, i, &pt);
    }

    lambdas = (Array_double*)malloc(sizeof(Array_double));
    Lambda_n = (Vector_double*)malloc(sizeof(Vector_double));
    Lambda_np1 = (Vector_double*)malloc(sizeof(Vector_double));
    compute_lambdas2D(clipper2D, clipped2D, vs, 1.0, \
                        &lambdas, &Lambda_n, &Lambda_np1, &normal, &is_narrowband);

    printf("lambdas = ");
    print_arr_double(lambdas);
    printf("\n");
    printf("Lambda_n = ");
    print_vec_double(Lambda_n);
    printf("\n");
    printf("Lambda_np1 = ");
    print_vec_double(Lambda_np1);
    printf("\n");
    printf("normal = ");
    print_pt3D(&normal);
    printf("\n");
    printf("is_narrowband = %d\n", is_narrowband);
}
*/

void test7(
    const int *signed_n_x,// = 2;
    const double x_v[],// = {0.0,1.0};
    const double x_v_clipped[],// = {0.5,1.5};
    const int *signed_n_y,// = 2;
    const double y_v[],// = {0.0,1.0};
    const double y_v_clipped[]// = {0.5,1.5};
){
    Polygon2D *clipper2D, *clipped2D;
    Vector_points2D *vs;
    Point2D pt = (Point2D){0.,0.};
    Array_double *lambdas;
    Vector_double *Lambda_n, *Lambda_np1;
    Point3D normal;
    bool is_narrowband;
    long unsigned i;
    const unsigned int n_x = (unsigned int) *signed_n_x;
    const unsigned int n_y = (unsigned int) *signed_n_y;

    clipper2D = polygon2D_from_vertices(x_v, n_x, y_v, n_y);
    clipped2D = polygon2D_from_vertices(x_v_clipped, n_x, y_v_clipped, n_y);
    vs = alloc_with_capacity_vec_pts2D(clipper2D->vertices->size);
    for(i=0; i<clipper2D->vertices->size; i++){
        set_ith_elem_vec_pts2D(vs, i, &pt);
    }

    lambdas = (Array_double*)malloc(sizeof(Array_double));
    Lambda_n = (Vector_double*)malloc(sizeof(Vector_double));
    Lambda_np1 = (Vector_double*)malloc(sizeof(Vector_double));
    compute_lambdas2D(clipper2D, clipped2D, vs, 1.0, \
                        &lambdas, &Lambda_n, &Lambda_np1, &normal, &is_narrowband);

    printf("lambdas = ");
    print_arr_double(lambdas);
    printf("\n");
    printf("Lambda_n = ");
    print_vec_double(Lambda_n);
    printf("\n");
    printf("Lambda_np1 = ");
    print_vec_double(Lambda_np1);
    printf("\n");
    printf("normal = ");
    print_pt3D(&normal);
    printf("\n");
    printf("is_narrowband = %d\n", is_narrowband);
}


/*
int main(){
    GrB_init(GrB_NONBLOCKING) ;
    //test1();
    //test2();
    //test3();
    //test4();
    //test4bis();
    //test5();
    //test6();
    test7();
    GrB_finalize();
}*/