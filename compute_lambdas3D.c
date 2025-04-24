#include "compute_lambdas3D.h"
#include <math.h>

static double norm(const Point3D* p){
    return sqrt(p->x*p->x + p->y*p->y + p->t*p->t);
}

static Polygon2D* extract_ith_face(const Polygon2D *p, uint64_t extr_index){
    GrB_Matrix *new_faces, *intermediate_faces, *new_edges;
    Vector_points2D *new_vertices;
    Vector_int *new_status_edge;
    GrB_Vector val_extr, J_vect, pt_indices;
    GrB_Index ncols_faces, ncols_edges;
    GrB_Index nrows_faces, nrows_edges;
    GrB_Index length_edge_ind;
    GrB_Index i, j, e_i;
    unsigned long int indpt1, indpt2;
    GrB_Index *I_;
    GrB_Index *edge_indices;
    uint8_t *nm;
    Vector_uint *pts_kept_inds;
    bool found;

    GrB_Matrix_ncols(&ncols_faces, *(p->faces));
    if (extr_index>=ncols_faces){
        printf("Can't extract face %ld of a polygon with only %ld face(s).", extr_index, ncols_faces);
        return NULL;
    }

    GrB_Matrix_nrows(&nrows_faces, *(p->faces));
    GrB_Matrix_ncols(&ncols_edges, *(p->edges));
    GrB_Matrix_nrows(&nrows_edges, *(p->edges));

    new_faces = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    intermediate_faces = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    new_edges = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    
    //First, extract new face.
    I_ = (GrB_Index*)malloc((nrows_faces>nrows_edges ? nrows_faces : nrows_edges)*sizeof(GrB_Index));
    for (i=0; i<nrows_faces; i++){
        I_[i] = i;
    }
    GrB_Matrix_new(intermediate_faces, GrB_INT8, nrows_faces, 1);
    GrB_Matrix_extract(*intermediate_faces, GrB_NULL, GrB_NULL, *(p->faces), I_, nrows_faces, (GrB_Index[]){extr_index}, 1, GrB_NULL);

    //Second, take all edge indices in the extracted face.
    GrB_Matrix_nvals(&length_edge_ind, *intermediate_faces);
    edge_indices = (GrB_Index*) calloc(length_edge_ind, sizeof(GrB_Index));
    nm = (uint8_t*) malloc(length_edge_ind*sizeof(uint8_t));
    GrB_Matrix_extractTuples(edge_indices, I_, nm, &length_edge_ind, *intermediate_faces);

    //Third, build a unique list of indices of points involved in kept edges
    GrB_Vector_new(&pt_indices, GrB_UINT64, 2);
    GrB_Vector_new(&J_vect, GrB_UINT64, nrows_edges);
    GrB_Vector_new(&val_extr, GrB_INT8, 2);
    pts_kept_inds = alloc_empty_vec_uint();
    for (i=0; i<length_edge_ind; i++){
        e_i = edge_indices[i];
        GrB_extract(J_vect, GrB_NULL, GrB_NULL, *(p->edges), GrB_ALL, 1, e_i, GrB_NULL);
        GrB_Vector_extractTuples(pt_indices, val_extr, J_vect, GrB_NULL);
        GrB_Vector_extractElement(&indpt1, pt_indices, 0);
        GrB_Vector_extractElement(&indpt2, pt_indices, 1);

        found = false;
        for(j=0;j<pts_kept_inds->size; j++){
            if (indpt1 == *get_ith_elem_vec_uint(pts_kept_inds, j)){
                found = true;
                break;
            }
        }
        if (!found){
            push_back_vec_uint(pts_kept_inds, &indpt1);
        }

        found = false;
        for(j=0;j<pts_kept_inds->size; j++){
            if (indpt2 == *get_ith_elem_vec_uint(pts_kept_inds, j)){
                found = true;
                break;
            }
        }
        if (!found){
            push_back_vec_uint(pts_kept_inds, &indpt2);
        }
    }

    //Build the new list of vertices
    new_vertices = alloc_with_capacity_vec_pts2D(pts_kept_inds->size);
    for (i=0; i<pts_kept_inds->size; i++){
        indpt1 = *get_ith_elem_vec_uint(pts_kept_inds, i);
        set_ith_elem_vec_pts2D(new_vertices, i, get_ith_elem_vec_pts2D(p->vertices, indpt1));
    }
    //Build the new list of status of each edge
    new_status_edge = alloc_with_capacity_vec_int(length_edge_ind);
    for (i=0; i<pts_kept_inds->size; i++){
        indpt1 = edge_indices[i];
        set_ith_elem_vec_int(new_status_edge, i, get_ith_elem_vec_int(p->status_edge, indpt1));
    }

    //Build the new array describing the edges
    GrB_Matrix_new(new_edges, GrB_INT8, pts_kept_inds->size, length_edge_ind);
    GrB_Matrix_extract(*new_edges, GrB_NULL, GrB_NULL, *(p->edges), pts_kept_inds->data, pts_kept_inds->size, edge_indices, length_edge_ind, GrB_NULL);

    //Re-extract new_faces in order to get rid of un-necessary zeros
    GrB_Matrix_new(new_faces, GrB_INT8, length_edge_ind, 1);
    GrB_Matrix_extract(*new_faces, GrB_NULL, GrB_NULL, *intermediate_faces, edge_indices, length_edge_ind, (GrB_Index[]){extr_index}, 1, GrB_NULL);

    free(I_);
    free(edge_indices);
    free(nm);
    dealloc_vec_uint(pts_kept_inds);

    return new_Polygon2D_vefs(new_vertices, new_edges, new_faces, new_status_edge);
}

static void build_clipped_in_partial(Polyhedron3D* p, Point3D *n_face, Point3D *pt_face, long int mark_edge, int8_t sign_taken){
    Vector_points3D *pts_copy = alloc_empty_vec_pts3D();
    Vector_int *status_face = alloc_empty_vec_int();
    GrB_Matrix *edges_in, *faces_in, *volumes_in;
    
    edges_in = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    faces_in = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    volumes_in = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));

    cut_edges3D(p, n_face, pt_face, sign_taken, pts_copy, edges_in);

    close_cells(edges_in, p->faces, NULL, -1, faces_in);

    copy_vec_int(p->status_face, status_face);
    close_cells(faces_in, p->volumes, status_face, mark_edge, volumes_in);

    //Copy new polyhedron in p.
    copy_vec_pts3D(pts_copy, p->vertices);
    dealloc_vec_pts3D(pts_copy);
    free(pts_copy);
    copy_vec_int(status_face, p->status_face);
    dealloc_vec_int(status_face);
    free(status_face);
    //free(p->edges);
    //p->edges = edges_in;
    //free(p->faces);
    //p->faces = faces_in;
    //free(p->volumes);
    //p->volumes = volumes_in;
    GrB_Matrix_dup(p->edges, *edges_in);
    GrB_Matrix_dup(p->faces, *faces_in);
    GrB_Matrix_dup(p->volumes, *volumes_in);
    GrB_free(edges_in);
    GrB_free(faces_in);
    GrB_free(volumes_in);
}

static Point3D* first_point_of_ith_face(const Polyhedron3D *p, GrB_Index i){
    Point3D* pt;
    GrB_Index nb_faces, e0, pt_index;
    GrB_Info infogrb;
    GrB_Vector whole_vec, nz, extr_vals;
    GrB_Index nrows;

    GrB_Matrix_ncols(&nb_faces, *(p->faces));
    if (i>=nb_faces){
        pt = (Point3D*)malloc(sizeof(Point3D));
        pt->x = 0./0.;
        pt->y = 0./0.;
        pt->t = 0./0.;
    } else {
        infogrb = GrB_Matrix_nrows(&nrows, *(p->faces));
        infogrb = GrB_Vector_new(&whole_vec, GrB_INT8, nrows);
        infogrb = GrB_Vector_new(&nz, GrB_UINT64, nrows);
        infogrb = GrB_Vector_new(&extr_vals, GrB_INT8, nrows);
        infogrb = GrB_extract(whole_vec, GrB_NULL, GrB_NULL, *(p->faces), GrB_ALL, 1, i, GrB_NULL); //Get indices of edges composing face i
        infogrb = GxB_Vector_extractTuples_Vector(nz, extr_vals, whole_vec, GrB_NULL);
        infogrb = GrB_Vector_extractElement(&e0, nz, 0);//Get index of first edge composing face i

        infogrb = GrB_Matrix_nrows(&nrows, *(p->edges));
        infogrb = GrB_Vector_resize(whole_vec, nrows);
        infogrb = GrB_Vector_resize(nz, nrows);
        infogrb = GrB_Vector_resize(extr_vals, nrows);
        infogrb = GrB_extract(whole_vec, GrB_NULL, GrB_NULL, *(p->edges), GrB_ALL, 1, e0, GrB_NULL); //Get indices of points of first edge composing face i
        infogrb = GxB_Vector_extractTuples_Vector(nz, extr_vals, whole_vec, GrB_NULL);
        infogrb = GrB_Vector_extractElement(&pt_index, nz, 0);//Get index of first point of first edge composing face i

        pt = get_ith_elem_vec_pts3D(p->vertices, pt_index);
    }

    return pt;
}

//    We suppose that clipper.status_face < 0 when the corresponding face is at t==t^n or t==t^{n+1}.
//    Return: Polyhedron3D
Polyhedron3D* clip3D(const Polyhedron3D *clipper, const Polyhedron3D *clipped){
    Polyhedron3D* clipped_in = new_Polyhedron3D();
    Vector_points3D* normal_vectors = points3D_from_matrix(surfaces_poly3D(clipper));
    GrB_Index i;//, size_nv;
    Point3D *n_face, *pt_face;
    long int *stf_i;
    int8_t vol_i;

    copy_Polyhedron3D(clipped, clipped_in);
    //for i in eachindex(normal_vectors)
    for (i = 0; i<normal_vectors->size; i++){
        stf_i = get_ith_elem_vec_int(clipper->status_face, i);
        if (*stf_i > 2){
            n_face = get_ith_elem_vec_pts3D(normal_vectors, i);
            pt_face = first_point_of_ith_face(clipper, i);
            GrB_Matrix_extractElement(&vol_i, *(clipper->volumes), i, 0);
            build_clipped_in_partial(clipped_in, n_face, pt_face, *stf_i, -vol_i);
        }
    }

    return clipped_in;
}

//    Compute the effective area on each face of grid.
//    The effective area is the area of each face minus the area intersected of initial_p.
//    The result consists of an ordered list of the faces of `grid` and 
//    λs, the effective area in a list ordered in the same order.
//    In λs, we store the effective area and the "complementarity", which is the area of the face minus the effective area.
//    Return : Vector{Vector{Number}}.
// *lambdas_vec is allocated with correct size inside the function.
void compute_lambdas2D_time(const Polyhedron3D* grid, const Polyhedron3D *initial_p, Vector_points3D **lambdas, Point3D *mean_normal, bool *is_narrowband){
    Polyhedron3D *p = clip3D(grid, initial_p);
    GrB_Index nb_edge, nb_cols_vol, nb_cols_fac, i, j;//, e;
    Point3D pt, *nvpi, *lam;
    int8_t pvij_int;
    long int *psfi;
    double pvij;
    Vector_points3D *norm_vec_poly;

    GrB_Matrix_ncols(&nb_edge, *(grid->faces)); //actually number of edges + 2 faces at times tn and tn+dt
    GrB_Matrix_ncols(&nb_cols_vol, *(p->volumes)); 
    GrB_Matrix_ncols(&nb_cols_fac, *(p->faces)); 

    *lambdas = alloc_with_capacity_vec_pts3D(nb_edge);
    for (j=0; j<nb_edge; j++){
        pt = (Point3D){0.,0.,0.};
        set_ith_elem_vec_pts3D(*lambdas, j, &pt);
    }

    norm_vec_poly = points3D_from_matrix(surfaces_poly3D(p));
    
    *mean_normal = (Point3D){0.0, 0.0, 0.0};
    *is_narrowband = false;
    for(i=0; i<nb_cols_fac; i++){
        psfi = get_ith_elem_vec_int(p->status_face, i);
        nvpi = get_ith_elem_vec_pts3D(norm_vec_poly, i);
        for (j=0; j<nb_cols_vol; j++){
            pvij_int = 0; //if volumes[i,j] does not exist, it won't change the value of pvij_int.
            GrB_Matrix_extractElement(&pvij_int, *(p->volumes), i, j);
            if (pvij_int != 0){
                if (*psfi<1){
                    pvij = (double) pvij_int;
                    inplace_xpay_points3D(mean_normal, pvij, nvpi);
                    *is_narrowband = true;
                }
                else {
                    lam = get_ith_elem_vec_pts3D(*lambdas, *psfi-1);
                    inplace_axpy_points3D(lam, 1.0, nvpi);
                }
                //for (e = 1; e<=nb_edge; e++){
                //    if (*psfi == ((long int) e)){
                //        lam = get_ith_elem_vec_pts3D(lambdas, e-1);
                //        inplace_axpy_points3D(lam, 1.0, nvpi);
                //    }
                //}
            }
        }
    }

    dealloc_vec_pts3D(norm_vec_poly);
    dealloc_Polyhedron3D(p);
}

//In status_faces, 1 denotes the face created at t^n, 2 the face created at t^{n+1}, and i+2 the face created in time-space from edge i.
Polyhedron3D* build_space2D_time_cell(const Polygon2D *fn, const Vector_points2D *vs, const double dt, bool split){
    uint64_t i;
    int8_t fne_ptind;
    long int val;
    Point2D *pt2D, *v;
    Point3D *pt3D;
    GrB_Index nb_edges, nb_rows_faces, nb_cols_faces, nb_cols_newedges, nb_cols_newfaces;
    GrB_Index curr_edge, v_edges_index, curr_face;
    GrB_Index pt_index1, pt_index2, temp;
    GrB_Index *I_index, *J_index;
    GrB_Info infogrb;
    GrB_Vector ed_i, nz_ei, val_nz_ei;
    GrB_Matrix emptyNE, emptySW, emptyE_N, emptyE_S;
    GrB_Matrix *new_edges = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    GrB_Matrix *new_faces = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    GrB_Matrix *new_volume = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    Vector_int *status_faces = (Vector_int*)malloc(sizeof(Vector_int));
    const uint64_t nb_pts = fn->vertices->size;
    Vector_points3D *new_vertices = alloc_with_capacity_vec_pts3D(2*nb_pts);
    
    pt3D = (Point3D*)malloc(sizeof(Point3D));
    for (i = 0; i<nb_pts; i++){
        pt2D = get_ith_elem_vec_pts2D(fn->vertices, i);
        *pt3D = (Point3D){pt2D->x, pt2D->y, 0.0};
        set_ith_elem_vec_pts3D(new_vertices, i, pt3D);
    }
    for (i = nb_pts; i<2*nb_pts; i++){
        pt2D = get_ith_elem_vec_pts2D(fn->vertices, i-nb_pts);
        v = get_ith_elem_vec_pts2D(vs, i-nb_pts);
        *pt3D = (Point3D){pt2D->x + dt*v->x, pt2D->y + dt*v->y, dt};
        set_ith_elem_vec_pts3D(new_vertices, i, pt3D);
    }

    GrB_Matrix_ncols(&nb_edges, *(fn->edges));

    GrB_Matrix_new(&emptyNE, GrB_INT8, nb_pts, nb_edges);
    GrB_Matrix_new(&emptySW, GrB_INT8, nb_pts, nb_edges);
    if (split){
        GrB_Matrix_new(&emptyE_N, GrB_INT8, nb_pts, 2*nb_pts);
        GrB_Matrix_new(&emptyE_S, GrB_INT8, nb_pts, 2*nb_pts);
        GrB_Matrix_new(new_edges, GrB_INT8, 2*nb_pts, 2*(nb_edges + nb_pts));
    } else {
        GrB_Matrix_new(&emptyE_N, GrB_INT8, nb_pts, nb_pts);
        GrB_Matrix_new(&emptyE_S, GrB_INT8, nb_pts, nb_pts);
        GrB_Matrix_new(new_edges, GrB_INT8, 2*nb_pts, 2*nb_edges + nb_pts);
    }
    GxB_Matrix_concat(*new_edges, (GrB_Matrix[]){*(fn->edges), emptyNE, emptyE_N, emptySW, *(fn->edges), emptyE_S}, 2, 3, GrB_NULL); //new_edges = [fn->edges, 0, 0; 0, fn->edges, 0]

    GrB_Matrix_ncols(&nb_cols_newedges, *new_edges);
    GrB_Matrix_ncols(&nb_cols_faces, *(fn->faces));
    

    if (split){
        GrB_Matrix_new(new_faces, GrB_INT8, nb_cols_newedges, 2*(nb_cols_faces + nb_edges));
    } else {
        GrB_Matrix_new(new_faces, GrB_INT8, nb_cols_newedges, 2*nb_cols_faces + nb_edges);
    }
    I_index = (GrB_Index*)malloc(nb_edges*sizeof(GrB_Index));
    J_index = (GrB_Index*)malloc(nb_cols_faces*sizeof(GrB_Index));
    for (i = 0; i<nb_edges; i++){
        I_index[i] = nb_edges + i;
    }
    for (i = 0; i<nb_cols_faces; i++){
        J_index[i] = nb_cols_faces + i;
    }
    GxB_Matrix_subassign(*new_faces, GrB_NULL, GrB_NULL, *(fn->faces), I_index, nb_edges, J_index, nb_cols_faces, GrB_NULL);
    GrB_apply(*new_faces, GrB_NULL, GrB_NULL, GrB_AINV_INT8, *new_faces, GrB_NULL);
    for (i = 0; i<nb_edges; i++){
        I_index[i] = i;
    }
    for (i = 0; i<nb_cols_faces; i++){
        J_index[i] = i;
    }
    GxB_Matrix_subassign(*new_faces, GrB_NULL, GrB_NULL, *(fn->faces), I_index, nb_edges, J_index, nb_cols_faces, GrB_NULL);

    GrB_Matrix_ncols(&nb_cols_newfaces, *new_faces);
    status_faces = alloc_with_capacity_vec_int(nb_cols_newfaces);
    for(i = 0; i < nb_cols_faces; i++){
        val = 1;
        set_ith_elem_vec_int(status_faces, i, &val);
    }
    for(i = nb_cols_faces; i < 2*nb_cols_faces; i++){
        val = 2;
        set_ith_elem_vec_int(status_faces, i, &val);
    }
    for(i = 2*nb_cols_faces; i < nb_cols_newfaces; i++){
        val = 0;
        set_ith_elem_vec_int(status_faces, i, &val);
    }

    GrB_Matrix_new(new_volume, GrB_INT8, nb_cols_newfaces, 1);
    for(i=0; i<nb_edges; i++){
        infogrb = GrB_Matrix_setElement(*new_volume, -1, i, 0);
    }
    for(i=nb_edges; i<nb_cols_newfaces; i++){
        infogrb = GrB_Matrix_setElement(*new_volume, 1, i, 0);
    }

    curr_edge = 2*nb_edges;
    v_edges_index = curr_edge;
    curr_face = 2*nb_cols_faces;

    for(i=0; i<nb_pts; i++){
        GrB_Matrix_setElement(*new_edges, -1, i, curr_edge);
        GrB_Matrix_setElement(*new_edges, 1, i + nb_pts, curr_edge);
        curr_edge++;
    }

    //At this point, new_edges has the following organisation :
    //First nb_edges: original edges of f at t=tn
    //Following nb_edges: edges of f at t=tn+dt
    //Following nb_pts: edges linking pt at t=tn and pt at t=tn+dt
    GrB_Matrix_nrows(&nb_rows_faces, *(fn->faces));
    infogrb = GrB_Vector_new(&ed_i, GrB_INT8, nb_rows_faces);
    infogrb = GrB_Vector_new(&nz_ei, GrB_UINT64, nb_rows_faces);
    infogrb = GrB_Vector_new(&val_nz_ei, GrB_INT8, nb_rows_faces);
    if (split){ //Each edge in 2D is turned into 2 faces (triangulation)
        //Now, we create the diagonal edges, and create all faces between t=tn and t=tn+dt
        for (i=0; i<nb_edges; i++){
            infogrb = GrB_extract(ed_i, GrB_NULL, GrB_NULL, *(fn->edges), GrB_ALL, 1, i, GrB_NULL); //Get point indices of edge i
            infogrb = GxB_Vector_extractTuples_Vector(nz_ei, val_nz_ei, ed_i, GrB_NULL);

            infogrb = GrB_Vector_extractElement(&pt_index1, nz_ei, 0);
            infogrb = GrB_Vector_extractElement(&pt_index2, nz_ei, 1);
            
            //GrB_Matrix_extractElement(&fne_ptind, *(fn->edges), pt_index1, i);
            GrB_Vector_extractElement(&fne_ptind, val_nz_ei, 0);
            if (fne_ptind > 0){ //pt_index1 should always be the index where the edges starts
                temp = pt_index1;
                pt_index1 = pt_index2;
                pt_index2 = temp;
            }

            GrB_Matrix_setElement(*new_edges,  1, pt_index2 + nb_pts, curr_edge); //diagonal edge
            GrB_Matrix_setElement(*new_edges, -1, pt_index1         , curr_edge); //diagonal edge

            GrB_Matrix_setElement(*new_faces,  1, i                        , curr_face); //edge at time tn
            GrB_Matrix_setElement(*new_faces,  1, v_edges_index + pt_index2, curr_face); //vertical edge
            GrB_Matrix_setElement(*new_faces, -1, curr_edge                , curr_face); //diagonal edge

            GrB_Matrix_setElement(*new_faces, -1, i + nb_edges             , curr_face + 1); //edge at time tn+dt
            GrB_Matrix_setElement(*new_faces, -1, v_edges_index + pt_index1, curr_face + 1); //vertical edge of pt_index1 between t=tn and t=tn+dt
            GrB_Matrix_setElement(*new_faces,  1, curr_edge                , curr_face + 1); //diagonal edge

            val = 3+i;
            set_ith_elem_vec_int(status_faces, curr_face    , &val);
            set_ith_elem_vec_int(status_faces, curr_face + 1, &val);

            GrB_Matrix_extractElement(&fne_ptind, *(fn->faces), i, 0);
            GrB_Matrix_setElement(*new_volume, fne_ptind, curr_face    , 0);
            GrB_Matrix_setElement(*new_volume, fne_ptind, curr_face + 1, 0);

            curr_edge += 1;
            curr_face += 2;
        }
    } else {//Each edge in 2D is turned into 1 single face (no triangulation)
        //Now, we create all faces between t=tn and t=tn+dt
        for (i=0; i<nb_edges; i++){
            infogrb = GrB_extract(ed_i, GrB_NULL, GrB_NULL, *(fn->edges), GrB_ALL, 1, i, GrB_NULL); //Get point indices of edge i
            infogrb = GxB_Vector_extractTuples_Vector(nz_ei, val_nz_ei, ed_i, GrB_NULL);

            infogrb = GrB_Vector_extractElement(&pt_index1, nz_ei, 0);
            infogrb = GrB_Vector_extractElement(&pt_index2, nz_ei, 1);
            
            GrB_Matrix_extractElement(&fne_ptind, *(fn->edges), pt_index1, i);
            if (fne_ptind > 0){ //pt_index1 should always be the index where the edges starts
                temp = pt_index1;
                pt_index1 = pt_index2;
                pt_index2 = temp;
            }

            GrB_Matrix_setElement(*new_faces,  1, i                        , curr_face); //edge at time tn
            GrB_Matrix_setElement(*new_faces, -1, i + nb_edges             , curr_face); //edge at time tn+dt
            GrB_Matrix_setElement(*new_faces,  1, v_edges_index + pt_index2, curr_face); //vertical edge
            GrB_Matrix_setElement(*new_faces, -1, v_edges_index + pt_index1, curr_face); //vertical edge of pt_index1 between t=tn and t=tn+dt

            val = 3+i;
            set_ith_elem_vec_int(status_faces, curr_face, &val);

            GrB_Matrix_extractElement(&fne_ptind, *(fn->faces), i, 0);
            GrB_Matrix_setElement(*new_volume, fne_ptind, curr_face, 0);

            curr_edge += 1;
            curr_face += 1;
        }
    }

    return new_Polyhedron3D_vefvs(new_vertices, new_edges, new_faces, new_volume, status_faces);
}

//    vec_move_solid should be an array of array of size #columns(clipped.faces), 
//    and each member i of this vector should be of size mini_clipped->vertices->size, 
//    where mini_clipped = extract_ith_face(clipped, i)
void compute_lambdas2D(const Polygon2D* grid, const Polygon2D *clipped, const Vector_points2D *vec_move_solid, const double dt, \
                        Array_double **lambdas_arr, Vector_double** big_lambda_n, Vector_double** big_lambda_np1, Point3D *mean_normal, bool *is_narrowband){
    const unsigned int nb_regions = 2;
    double *val = (double*)malloc(sizeof(double));
    double nm;
    Point3D *area, *occupied;
    Vector_points3D *occupied_area;
    Polygon2D *mini_clipped;
    Polyhedron3D *clipped3D, *cell3D;
    long *sfj;
    const Point2D pt2D = (Point2D){0.,0.};
    Point3D *pt3D, local_mean_normal, *local_l;
    bool local_narrowband;
    Vector_points2D *vec_move_grid;
    Vector_points3D *surfaces;
    GrB_Index nb_edge, i, j, k, nb_clipped_faces;
    Vector_points3D *lambdas3D; //TODO : Change this when nb_regions>2
    Array_points3D *local_lambdas;
    
    GrB_Matrix_ncols(&nb_edge, *(grid->edges));
    *lambdas_arr = alloc_with_capacity_arr_double(nb_edge, nb_regions); //All set to 0
    local_lambdas = alloc_with_capacity_arr_pts3D(nb_edge + 2, nb_regions-1); //All set to 0
    occupied_area = alloc_with_capacity_vec_pts3D(nb_edge + 2);

    *big_lambda_n = alloc_with_capacity_vec_double(nb_regions);
    *big_lambda_np1 = alloc_with_capacity_vec_double(nb_regions);
    *val = 0.;
    for (i=0; i<nb_regions; i++){
        set_ith_elem_vec_double(*big_lambda_n, i, val);
        set_ith_elem_vec_double(*big_lambda_np1, i, val);
    }
    
    *mean_normal = (Point3D){0.,0.,0.};
    for (i=0; i<nb_edge + 2; i++){
        set_ith_elem_vec_pts3D(occupied_area, i, mean_normal);
    }

    *is_narrowband = false;
    vec_move_grid = alloc_with_capacity_vec_pts2D(grid->vertices->size);
    for (i=0; i<grid->vertices->size; i++){
        set_ith_elem_vec_pts2D(vec_move_grid, i, &pt2D);
    }
    cell3D = build_space2D_time_cell(grid, vec_move_grid, dt, false);
    surfaces = points3D_from_matrix(surfaces_poly3D(cell3D));

    if (clipped->vertices->size>2){
        GrB_Matrix_ncols(&nb_clipped_faces, *(clipped->faces));
        for (i=0; i<nb_clipped_faces; i++){
            mini_clipped = extract_ith_face(clipped, i);
            k = 0; //should be k = clipped->status_edge[i], or another variable to indicate what region covers face nb i.
            clipped3D = build_space2D_time_cell(mini_clipped, vec_move_solid + i, dt, true);
            for (j=0; j<clipped3D->status_face->size; j++){
                sfj = get_ith_elem_vec_int(clipped3D->status_face, j);
                if (*sfj > 2)   *sfj = -1;
            }
            
            compute_lambdas2D_time(cell3D, clipped3D, &lambdas3D, &local_mean_normal, &local_narrowband);
            dealloc_Polyhedron3D(clipped3D);
            //λ, ni, is_na = compute_lambdas(cell3D, clipped3D)

            mean_normal->x += local_mean_normal.x;
            mean_normal->y += local_mean_normal.y;
            mean_normal->t += local_mean_normal.t;
            *is_narrowband |= local_narrowband;
            
            //λe += λ[2:end]
            for (j=0; j<lambdas3D->size; j++){
                pt3D = get_ith_elem_vec_pts3D(lambdas3D, j); 
                local_l = get_ijth_elem_arr_pts3D(local_lambdas, j, k);
                local_l->x += pt3D->x;
                local_l->y += pt3D->y;
                local_l->t += pt3D->t;

                local_l = get_ith_elem_vec_pts3D(occupied_area, j);
                local_l->x += pt3D->x;
                local_l->y += pt3D->y;
                local_l->t += pt3D->t;
            }
            dealloc_Polygon2D(mini_clipped);
        } 

        i = 0;
        area = get_ith_elem_vec_pts3D(surfaces, i);
        occupied = get_ith_elem_vec_pts3D(occupied_area, i);
        *val = fmax(0., norm(area) - norm(occupied));

        set_ith_elem_vec_double(*big_lambda_n, 0, val);
        for(k = 1; k<nb_regions; k++){
            nm = norm(get_ijth_elem_arr_pts3D(local_lambdas, i, k-1));
            set_ith_elem_vec_double(*big_lambda_n, k, &nm);
        }
        

        i = 1;
        area = get_ith_elem_vec_pts3D(surfaces, i);
        occupied = get_ith_elem_vec_pts3D(occupied_area, i);
        *val = fmax(0., norm(area) - norm(occupied));
        set_ith_elem_vec_double(*big_lambda_np1, 0, val);
        for(k = 1; k<nb_regions; k++){
            nm = norm(get_ijth_elem_arr_pts3D(local_lambdas, i, k-1));
            set_ith_elem_vec_double(*big_lambda_np1, k, &nm);
        }

        for (i=2; i<nb_edge + 2; i++){
            area = get_ith_elem_vec_pts3D(surfaces, i);
            occupied = get_ith_elem_vec_pts3D(occupied_area, i);
            *val = fmax(0., norm(area) - norm(occupied));
            set_ijth_elem_arr_double(*lambdas_arr, i-2, 0, val);
            for(k = 1; k<nb_regions; k++){
                nm = norm(get_ijth_elem_arr_pts3D(local_lambdas, i, k-1));
                set_ijth_elem_arr_double(*lambdas_arr, i-2, k, &nm); 
            }
        }
    } else {
        i = 0;
        *val = norm(get_ith_elem_vec_pts3D(surfaces, i));
        set_ith_elem_vec_double(*big_lambda_n, 0, val);
        nm = 0.;
        for(k = 1; k<nb_regions; k++){
            set_ith_elem_vec_double(*big_lambda_n, k, &nm);
        }

        i = 1;
        *val = norm(get_ith_elem_vec_pts3D(surfaces, i));
        set_ith_elem_vec_double(*big_lambda_np1, 0, val);
        for(k = 1; k<nb_regions; k++){
            set_ith_elem_vec_double(*big_lambda_np1, k, &nm);
        }

        for (i=2; i<nb_edge + 2; i++){
            *val = norm(get_ith_elem_vec_pts3D(surfaces, i));
            set_ijth_elem_arr_double(*lambdas_arr, i-2, 0, val);

            for(k = 1; k<nb_regions; k++){
                set_ijth_elem_arr_double(*lambdas_arr, i-2, k, &nm); 
            }
        }
    }

    dealloc_Polyhedron3D(cell3D);
    dealloc_vec_pts2D(vec_move_grid);
    dealloc_vec_pts3D(surfaces);
}

/*
//Fusion of cells with "too small" fluid region. Searches for neighbour cells with a big enough ratio of fluid and fuses them.\n
//Returns two arrays: \n
//    - one that gives, for each cell, the cell with which it will fuse \n
//    - one that tells the status of the cell with one of the flag in FusingCellType \n
//Returns: Vector{Int}, Vector{FusingCellType}
function fuse_cells(grid::Grid2D, threshold::Real, is_narrowband::Vector{Bool}, λ_per_edge::Matrix{<:Number}, 
                    Λn::Array{Number}, Λnp1::Array{Number})
    target_cells = Vector{Vector{Int}}(undef, size(grid.p.faces, 2))
    cell_type = Vector{FusingCellType}(undef, size(grid.p.faces, 2))

    narrowBand1 = [((0<Λn[i, 1]/grid.areas[i]<threshold)||(0<Λnp1[i,1]/grid.areas[i]<threshold))&&is_narrowband[i] for i in axes(grid.p.faces, 2)]
    narrowBand2 = [((0<Λn[i, 2]/grid.areas[i]<threshold)||(0<Λnp1[i,2]/grid.areas[i]<threshold))&&is_narrowband[i] for i in axes(grid.p.faces, 2)]

    edge_indices = rowvals(grid.p.faces)
    cells_to_be_merged1 = Set{Int}()
    cells_to_be_merged2 = Set{Int}()
    for f in axes(grid.p.faces, 2)
        if narrowBand1[f]
            push!(cells_to_be_merged1, f)
            cell_type[f] = PROBLEMATIC_UNSOLVED
            target_cells[f] = Vector{Int}()
            krange = nzrange(grid.p.faces, f)
            for e in edge_indices[krange]
                if (λ_per_edge[e, 1]>0.0)
                    faces = adjacency_edge(grid, e)
                    if !isnothing(faces) && faces[2]>0
                        adjf = faces[1]
                        if adjf == f
                            adjf = faces[2]
                        end
                        if (Λn[adjf, 1]>0.0) || (Λnp1[adjf, 1]>0.0)
                            push!(target_cells[f], adjf)
                        end
                    end
                end 
            end
        elseif narrowBand2[f]
            push!(cells_to_be_merged2, f)
            cell_type[f] = PROBLEMATIC_UNSOLVED
            target_cells[f] = Vector{Int}()
            krange = nzrange(grid.p.faces, f)
            for e in edge_indices[krange]
                if (λ_per_edge[e,2]>0.0)
                    faces = adjacency_edge(grid, e)
                    if !isnothing(faces) && faces[2]>0
                        adjf = faces[1]
                        if adjf == f
                            adjf = faces[2]
                        end
                        if (Λn[adjf,2]>0.0) || (Λnp1[adjf,2]>0.0)
                            push!(target_cells[f], adjf)
                        end
                    end
                end 
            end
        else
            cell_type[f] = NOT_FUSED
            target_cells[f] = [f]
        end
    end
    
    global_fusion_happened = true 
    while global_fusion_happened
        global_fusion_happened = false

        fusion_happened = true
        while fusion_happened
            fusion_happened = false
            copy_to_be_merged = copy(cells_to_be_merged1)
            for c in copy_to_be_merged
                targ_c = [tc for tc in target_cells[c] if (cell_type[tc]==NOT_FUSED || cell_type[tc]==FUSED)]
                if !isempty(targ_c)
                    Λs = [Λn[f,1]/grid.areas[f] for f in targ_c]
                    i = argmax(Λs)
                    target_cells[c] = [targ_c[i]]
                    delete!(cells_to_be_merged1, c)
                    cell_type[c] = FUSED
                    cell_type[targ_c[i]] = FUSED
                    fusion_happened = true
                    global_fusion_happened = true
                end
            end
        end

        fusion_happened = true
        while fusion_happened
            fusion_happened = false
            copy_to_be_merged = copy(cells_to_be_merged2)
            for c in copy_to_be_merged
                targ_c = [tc for tc in target_cells[c] if (cell_type[tc]==NOT_FUSED || cell_type[tc]==FUSED)]
                if !isempty(targ_c)
                    Λs = [Λn[f,2]/grid.areas[f] for f in targ_c]
                    i = argmax(Λs)
                    target_cells[c] = [targ_c[i]]
                    delete!(cells_to_be_merged2, c)
                    cell_type[c] = FUSED
                    cell_type[targ_c[i]] = FUSED
                    fusion_happened = true
                    global_fusion_happened = true
                end
            end
        end
    end

    while !isempty(cells_to_be_merged1) //Problem : some cells can't be merged with a cell big enough.
        rand_cell = rand(cells_to_be_merged1)
        target_cells[rand_cell] = [rand_cell]
        cell_type[rand_cell] = PROBLEMATIC
        delete!(cells_to_be_merged1, rand_cell)
        fusion_happened = true
        while fusion_happened
            fusion_happened = false
            copy_to_be_merged = copy(cells_to_be_merged1)
            for c in copy_to_be_merged
                targ_c = [tc for tc in target_cells[c] if (cell_type[tc]==PROBLEMATIC)]
                if !isempty(targ_c)
                    Λs = [Λn[f,1]/grid.areas[f] for f in targ_c]
                    i = argmax(Λs)
                    target_cells[c] = [targ_c[i]]
                    delete!(cells_to_be_merged1, c)
                    cell_type[c] = PROBLEMATIC
                    fusion_happened = true
                end
            end
        end
    end
    while !isempty(cells_to_be_merged2) //Problem : some cells can't be merged with a cell big enough.
        rand_cell = rand(cells_to_be_merged2)
        target_cells[rand_cell] = [rand_cell]
        cell_type[rand_cell] = PROBLEMATIC
        delete!(cells_to_be_merged2, rand_cell)
        fusion_happened = true
        while fusion_happened
            fusion_happened = false
            copy_to_be_merged = copy(cells_to_be_merged2)
            for c in copy_to_be_merged
                targ_c = [tc for tc in target_cells[c] if (cell_type[tc]==PROBLEMATIC)]
                if !isempty(targ_c)
                    Λs = [Λn[f,2]/grid.areas[f] for f in targ_c]
                    i = argmax(Λs)
                    target_cells[c] = [targ_c[i]]
                    delete!(cells_to_be_merged2, c)
                    cell_type[c] = PROBLEMATIC
                    fusion_happened = true
                end
            end
        end
    end

    final_target_cells = Vector{Int}(undef, size(grid.p.faces, 2))
    for f in axes(grid.p.faces, 2)
        final_target_cells[f] = target_cells[f][1]
    end

    for f in axes(grid.p.faces, 2)
        ftc = final_target_cells[f]
        while ftc != final_target_cells[ftc]
            ftc = final_target_cells[ftc]
        end
        final_target_cells[f] = ftc
        if f == ftc && cell_type[f] == FUSED
            cell_type[f] = TARGET_FUSED
        end
    end

    return final_target_cells, cell_type
end
*/
