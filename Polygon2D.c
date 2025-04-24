#include "Polygon2D.h"

Polygon2D* new_Polygon2D(){
    GrB_Info infogrb;
    Polygon2D* p = (Polygon2D*) malloc(sizeof(Polygon2D));

    p->vertices = alloc_empty_vec_pts2D();
    p->edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    p->faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    infogrb = GrB_Matrix_new(p->edges, GrB_INT8, 0,0);
    infogrb = GrB_Matrix_new(p->faces, GrB_INT8, 0,0);
    p->status_edge = alloc_empty_vec_int();
    
    return p;
}

Polygon2D* new_Polygon2D_ves(Vector_points2D* vertices, GrB_Matrix* edges, Vector_int* status_edge){
    GrB_Info infogrb;
    GrB_Index nb_rows;
    Polygon2D* p = (Polygon2D*) malloc(sizeof(Polygon2D));

    p->vertices = vertices;
    p->edges = edges;
    p->faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));

    infogrb = GrB_Matrix_ncols(&nb_rows, *(p->edges));
    infogrb = GrB_Matrix_new(p->faces, GrB_INT8, nb_rows, 1);
    p->status_edge = status_edge;
    
    return p;
}

Polygon2D* new_Polygon2D_vefs(Vector_points2D* vertices, GrB_Matrix* edges, GrB_Matrix* faces, Vector_int* status_edge){
    Polygon2D* p = (Polygon2D*) malloc(sizeof(Polygon2D));

    p->vertices = vertices;
    p->edges = edges;
    p->faces = faces;
    p->status_edge = status_edge;
    
    return p;
}

Polygon2D* polygon2D_from_vertices(const double* x_v, unsigned long int n_x, const double* y_v, unsigned long int n_y){
    GrB_Info infogrb;
    GrB_Index nb_pts, nb_edges, nb_faces;
    uint64_t curr_face, curr_pt, curr_edge, iy, ix;
    uint64_t ind_pt_SW, ind_pt_SE, ind_pt_NW;
    uint64_t ind_e_W, ind_e_S, ind_e_E, ind_e_N;
    double yS, xW;
    Point2D ptSW;
    Vector_points2D* vertices;
    Vector_int* status_edge;
    GrB_Matrix* edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix* faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));

    if (n_x<2 || n_y<2)
        return new_Polygon2D();

    nb_pts = n_x*n_y;
    nb_edges = (2*(n_x-1) + 1)*(n_y-1) + (n_y-1);
    nb_faces = (n_x-1)*(n_y-1);

    vertices = alloc_with_capacity_vec_pts2D(nb_pts);
    infogrb = GrB_Matrix_new(edges, GrB_INT8, nb_pts, nb_edges);
    infogrb = GrB_Matrix_new(faces, GrB_INT8, nb_edges, nb_faces);

    ptSW = (Point2D){x_v[0], y_v[0]};
    push_back_vec_pts2D(vertices, &ptSW);
    curr_face = 0;
    curr_pt = 0;
    curr_edge = 0;
    for (iy = 0; iy < n_y-2; iy++){
        yS = y_v[iy];
        for (ix = 0; ix < n_x-1; ix++){
            ind_pt_SW = curr_pt;
            ind_pt_SE = curr_pt + 1;
            ind_pt_NW = curr_pt + n_x;

            ind_e_W = curr_edge;
            ind_e_S = curr_edge + 1;
            ind_e_E = curr_edge + 2;
            ind_e_N = curr_edge + (2*(n_x-1) + 2);

            xW = x_v[ix];
            ptSW = (Point2D){xW,yS};

            set_ith_elem_vec_pts2D(vertices, ind_pt_SW, &ptSW);

            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_W);
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW, ind_e_W);
            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_S);
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE, ind_e_S);

            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S, curr_face);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_E, curr_face);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_N, curr_face);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W, curr_face);
        
            curr_pt += 1;
            curr_edge += 2;
            curr_face += 1;
        }
        //Treat east-most point and edge
        xW = x_v[n_x-1];
        ptSW = (Point2D){xW,yS};
        ind_pt_SW = curr_pt;
        ind_pt_NW = curr_pt + n_x;
        ind_e_W = curr_edge;

        set_ith_elem_vec_pts2D(vertices, ind_pt_SW, &ptSW);

        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_W);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW, ind_e_W);

        curr_pt += 1;
        curr_edge += 1;
    }

    //Treat north-most cells
    iy = n_y-2;
    yS = y_v[iy];
    ind_e_N = curr_edge + (2*(n_x-1));
    for (ix = 0; ix < n_x-1; ix++){
        ind_pt_SW = curr_pt;
        ind_pt_SE = curr_pt + 1;
        ind_pt_NW = curr_pt + n_x;

        ind_e_W = curr_edge;
        ind_e_S = curr_edge + 1 ;
        ind_e_E = curr_edge + 2;
        ind_e_N += 1;

        xW = x_v[ix];
        ptSW = (Point2D){xW,yS};
        set_ith_elem_vec_pts2D(vertices, ind_pt_SW, &ptSW);

        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_W);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW, ind_e_W);
        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_S);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE, ind_e_S);

        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S, curr_face);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_E, curr_face);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_N, curr_face);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W, curr_face);

        curr_pt += 1;
        curr_edge += 2;
        curr_face += 1;
    }
    //Treat north-east-most point and edge
    xW = x_v[n_x - 1];
    ptSW = (Point2D){xW,yS};
    ind_pt_SW = curr_pt;
    ind_pt_NW = curr_pt + n_x;
    ind_e_W = curr_edge;

    set_ith_elem_vec_pts2D(vertices, ind_pt_SW, &ptSW);

    infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_W);
    infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW, ind_e_W);

    curr_pt += 1;
    curr_edge += 1;
    

    //Treat north-most points and edges
    yS = y_v[n_y - 1];
    for (ix = 0; ix < n_x-1; ix++){
        ind_pt_SW = curr_pt;
        ind_pt_SE = curr_pt + 1;

        ind_e_S = curr_edge;

        xW = x_v[ix];
        ptSW = (Point2D){xW,yS};

        set_ith_elem_vec_pts2D(vertices, ind_pt_SW, &ptSW);

        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_S);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE, ind_e_S);

        curr_pt += 1;
        curr_edge += 1;
    }
    xW = x_v[n_x-1];
    ptSW = (Point2D){xW,yS};
    ind_pt_SW = curr_pt;
    set_ith_elem_vec_pts2D(vertices, ind_pt_SW, &ptSW);

    status_edge = alloc_empty_vec_int();
    status_edge->data = (long int*) calloc(nb_edges, sizeof(long int)); //all zeros
    status_edge->capacity = nb_edges;
    status_edge->size = nb_edges;
    return new_Polygon2D_vefs(vertices, edges, faces, status_edge);
}

void copy_Polygon2D(const Polygon2D *src, Polygon2D *dest){
    if (src == NULL){
        if (dest != NULL) {
            dealloc_vec_pts2D(dest->vertices);
            GrB_free(dest->edges);
            GrB_free(dest->faces);
            dealloc_vec_int(dest->status_edge);
            dest->vertices = NULL;
            dest->edges = NULL;
            dest->faces = NULL;
            dest->status_edge = NULL;
        }
    } else {
        if (dest != NULL){
            GrB_free(dest->edges);
            GrB_free(dest->faces);
            copy_vec_pts2D(src->vertices, dest->vertices);
            GrB_Matrix_dup(dest->edges, *(src->edges));
            GrB_Matrix_dup(dest->faces, *(src->faces));
            copy_vec_int(src->status_edge, dest->status_edge);
        }
    }
}

void dealloc_Polygon2D(Polygon2D* p){
    dealloc_vec_pts2D(p->vertices);
    GrB_free(p->edges);
    GrB_free(p->faces);
    dealloc_vec_int(p->status_edge);
}