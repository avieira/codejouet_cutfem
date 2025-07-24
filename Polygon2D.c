#include "Polygon2D.h"

Polygon2D* new_Polygon2D(){
    GrB_Info infogrb;
    Polygon2D* p = (Polygon2D*) malloc(sizeof(Polygon2D));

    p->vertices = alloc_empty_vec_pts2D();
    p->edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    p->faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    infogrb = GrB_Matrix_new(p->edges, GrB_INT8, 1,1);
    infogrb = GrB_Matrix_new(p->faces, GrB_INT8, 1,1);
    p->status_edge = alloc_empty_vec_int();
    //p->pressure_edge = alloc_empty_vec_double();
    
    return p;
}

Polygon2D* new_Polygon2D_ves(Vector_points2D* vertices, GrB_Matrix* edges, Vector_int* status_edge){
    GrB_Info infogrb;
    GrB_Index nb_rows;
    Polygon2D* p = (Polygon2D*) malloc(sizeof(Polygon2D));

    p->vertices = alloc_empty_vec_pts2D();
    copy_vec_pts2D(vertices, p->vertices);
    p->edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix_dup(p->edges, *edges);
    p->faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));

    infogrb = GrB_Matrix_ncols(&nb_rows, *(p->edges));
    infogrb = GrB_Matrix_new(p->faces, GrB_INT8, nb_rows, 1);

    p->status_edge = alloc_empty_vec_int();
    copy_vec_int(status_edge, p->status_edge);

    //p->pressure_edge = alloc_empty_vec_double();
    //p->pressure_edge->data = (my_real*) calloc(nb_rows, sizeof(my_real)); //all zeros
    //p->pressure_edge->capacity = nb_rows;
    //p->pressure_edge->size = nb_rows;

    return p;
}

Polygon2D* new_Polygon2D_vefs(Vector_points2D* vertices, GrB_Matrix* edges, GrB_Matrix* faces, Vector_int* status_edge){
    //GrB_Index nb_rows;
    Polygon2D* p = (Polygon2D*) malloc(sizeof(Polygon2D));

    p->vertices = alloc_empty_vec_pts2D();
    copy_vec_pts2D(vertices, p->vertices);
    p->edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix_dup(p->edges, *edges);
    p->faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix_dup(p->faces, *faces);
    p->status_edge = alloc_empty_vec_int();
    copy_vec_int(status_edge, p->status_edge);

    //GrB_Matrix_ncols(&nb_rows, *(p->edges));
    //p->pressure_edge = alloc_empty_vec_double();
    //p->pressure_edge->data = (my_real*) calloc(nb_rows, sizeof(my_real)); //all zeros
    //p->pressure_edge->capacity = nb_rows;
    //p->pressure_edge->size = nb_rows;
    
    return p;
}

Polygon2D* polygon2D_from_vertices(const my_real* x_v, unsigned long int n_x, const my_real* y_v, unsigned long int n_y){
    GrB_Info infogrb;
    GrB_Index nb_pts, nb_edges, nb_faces;
    uint64_t curr_face, curr_pt, curr_edge, iy, ix;
    uint64_t ind_pt_SW, ind_pt_SE, ind_pt_NW;
    uint64_t ind_e_W, ind_e_S, ind_e_E, ind_e_N;
    my_real yS, xW;
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

Polygon2D* polygon_from_consecutive_points(const my_real *x_v, const my_real* y_v, unsigned long int nb_pts){
    GrB_Info infogrb;
    unsigned long int nb_edges, nb_faces;
    Vector_points2D* vertices;
    Vector_int* status_edge;
    GrB_Matrix* edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix* faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    Point2D pt;
    unsigned long int i;
    long int zero = 0;

    if (nb_pts < 3)
        return new_Polygon2D();

    nb_edges = nb_pts;
    nb_faces = 1;

    vertices = alloc_with_capacity_vec_pts2D(nb_pts);
    status_edge = alloc_with_capacity_vec_int(nb_pts);
    infogrb = GrB_Matrix_new(edges, GrB_INT8, nb_pts, nb_edges);
    infogrb = GrB_Matrix_new(faces, GrB_INT8, nb_edges, nb_faces);

    for(i=0; i<nb_pts; i++){
        pt = (Point2D){x_v[i], y_v[i]};
        push_back_vec_pts2D(vertices, &pt);
        push_back_vec_int(status_edge, &zero);
    }

    for(i=0; i<(nb_edges-1); i++){
        infogrb = GrB_Matrix_setElement(*edges, -1, i, i);
        infogrb = GrB_Matrix_setElement(*edges, 1, i+1, i);
    }
    infogrb = GrB_Matrix_setElement(*edges, -1, nb_pts-1, nb_pts-1);
    infogrb = GrB_Matrix_setElement(*edges, 1, 0, nb_pts-1);

    for(i=0; i<nb_edges; i++)
        infogrb = GrB_Matrix_setElement(*faces, 1, i, 0);
    
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
            dest->edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
            dest->faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
            copy_vec_pts2D(src->vertices, dest->vertices);
            GrB_Matrix_dup(dest->edges, *(src->edges));
            GrB_Matrix_dup(dest->faces, *(src->faces));
            copy_vec_int(src->status_edge, dest->status_edge);
        }
    }
}

void dealloc_Polygon2D(Polygon2D* p){
    if (p){
        dealloc_vec_pts2D(p->vertices);
        GrB_free(p->edges);
        GrB_free(p->faces);
        dealloc_vec_int(p->status_edge);
    }
}

void retrieve_ith_edge2D(Vector_points2D *vertices, GrB_Matrix* edges, int i, Point2D** e_ext0, Point2D** e_ext1){
    GrB_Info infogrb;
    GrB_Index size_nz_ei;
    GrB_Vector ei;
    GrB_Vector nz_ei;
    GrB_Vector extr_vals_ei;
    GrB_Index i_pt0, i_pt1;
    GrB_Index nb_edges;

    GrB_Matrix_ncols(&nb_edges, *edges);
    infogrb = GrB_Vector_new(&ei, GrB_INT8, nb_edges);
    infogrb = GrB_Vector_new(&nz_ei, GrB_UINT64, 2);
    infogrb = GrB_Vector_new(&extr_vals_ei, GrB_INT8, 2);
    
    infogrb = GrB_extract(ei, GrB_NULL, GrB_NULL, *edges, GrB_ALL, 1, i, GrB_NULL); //Get indices of points composing edge i
    infogrb = GxB_Vector_extractTuples_Vector(nz_ei, extr_vals_ei, ei, GrB_NULL);
    infogrb = GrB_Vector_size(&size_nz_ei, nz_ei);
    if (size_nz_ei>0){
        infogrb = GrB_Vector_extractElement(&i_pt0, nz_ei, 0);
        infogrb = GrB_Vector_extractElement(&i_pt1, nz_ei, 1);
        *e_ext0 = get_ith_elem_vec_pts2D(vertices, i_pt0);
        *e_ext1 = get_ith_elem_vec_pts2D(vertices, i_pt1);
    }

    GrB_free(&ei);
    GrB_free(&nz_ei);
    GrB_free(&extr_vals_ei);
}

void dropzeros(GrB_Matrix* M){
    //GrB_assign(*M, *M, GrB_NULL, *M, GrB_ALL, 0, GrB_ALL, 0, GrB_DESC_R); //In case GxB_select gets changed in the future, this line will always work
    GxB_select(*M, NULL, NULL, GxB_NONZERO, *M, NULL, NULL) ;
}

Polygon2D* fuse_polygons(Polygon2D* p1, Polygon2D* p2){
    Vector_points2D* fused_vertices;
    GrB_Matrix *fused_edges = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    GrB_Matrix *fused_faces = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    Vector_int *fused_status;
    GrB_Index nrows1, ncols1, nrows2, ncols2;
    GrB_Matrix zeros1, zeros2;

    fused_vertices = cat_vec_pts2D(p1->vertices, p2->vertices);
    fused_status = cat_vec_int(p1->status_edge, p2->status_edge);

    GrB_Matrix_nrows(&nrows1, *(p1->edges));
    GrB_Matrix_ncols(&ncols1, *(p1->edges));
    GrB_Matrix_nrows(&nrows2, *(p2->edges));
    GrB_Matrix_ncols(&ncols2, *(p2->edges));
    GrB_Matrix_new(&zeros1, GrB_INT8, nrows1, ncols2);
    GrB_Matrix_new(&zeros2, GrB_INT8, nrows2, ncols1);
    GrB_Matrix_new(fused_edges, GrB_INT8, nrows1+nrows2, ncols1+ncols2);
    GxB_Matrix_concat(*fused_edges, (GrB_Matrix[]){*(p1->edges), zeros1, zeros2, *(p2->edges)}, 2, 2, GrB_NULL);
    //fused_edges = [[p1.edges spzeros(Int8, size(p1.edges, 1), size(p2.edges, 2))];
    //                [spzeros(Int8, size(p2.edges, 1), size(p1.edges, 2)) p2.edges]];

    GrB_Matrix_nrows(&nrows1, *(p1->faces));
    GrB_Matrix_ncols(&ncols1, *(p1->faces));
    GrB_Matrix_nrows(&nrows2, *(p2->faces));
    GrB_Matrix_ncols(&ncols2, *(p2->faces));
    GrB_Matrix_resize(zeros1, nrows1, ncols2);
    GrB_Matrix_resize(zeros2, nrows2, ncols1);
    GrB_Matrix_new(fused_faces, GrB_INT8, nrows1+nrows2, ncols1+ncols2);
    GxB_Matrix_concat(*fused_faces, (GrB_Matrix[]){*(p1->faces), zeros1, zeros2, *(p2->faces)}, 2, 2, GrB_NULL);
    //fused_faces = [[p1.faces  spzeros(Int8, size(p1.faces, 1), size(p2.faces, 2))];
    //                [spzeros(Int8, size(p2.faces, 1), size(p1.faces, 2)) p2.faces]];
    //fused_pressure = vcat(p1.pressure_edge, p2.pressure_edge) //TODO Report this
    
    GrB_free(&zeros1);
    GrB_free(&zeros2);
    
    return new_Polygon2D_vefs(fused_vertices, fused_edges, fused_faces, fused_status);
}