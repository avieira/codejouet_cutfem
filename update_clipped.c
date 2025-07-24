#include "update_clipped.h"
#include <string.h>

static my_real compute_angle_between_edges2D(Point2D *e1_ext1, Point2D *e1_ext2, Point2D *e2_ext1, Point2D *e2_ext2){
        my_real n1x, n1y, n2x, n2y;
        my_real dotp, norm_n1, norm_n2;
        if ((e1_ext1->x == e2_ext1->x) && (e1_ext1->y == e2_ext1->y)){
            n1x = e1_ext2->x - e1_ext1->x;
            n1y = e1_ext2->y - e1_ext1->y;
            n2x = e2_ext2->x - e2_ext1->x;
            n2y = e2_ext2->y - e2_ext1->y;
        } else if ((e1_ext1->x == e2_ext2->x) && (e1_ext1->y == e2_ext2->y)) {
            n1x = e1_ext2->x - e1_ext1->x;
            n1y = e1_ext2->y - e1_ext1->y;
            n2x = e2_ext1->x - e2_ext2->x;
            n2y = e2_ext1->y - e2_ext2->y;
        } else if (e1_ext2 == e2_ext1){
            n1x = e1_ext1->x - e1_ext2->x;
            n1y = e1_ext1->y - e1_ext2->y;
            n2x = e2_ext2->x - e2_ext1->x;
            n2y = e2_ext2->y - e2_ext1->y;
        } else { //e1_ext2 == e2_ext2
            n1x = e1_ext1->x - e1_ext2->x;
            n1y = e1_ext1->y - e1_ext2->y;
            n2x = e2_ext1->x - e2_ext2->x;
            n2y = e2_ext1->y - e2_ext2->y;
        }

        dotp = n1x*n2x + n1y*n2y;
        norm_n1 = sqrt(n1x*n1x + n1y*n1y);
        norm_n2 = sqrt(n2x*n2x + n2y*n2y);
    #if (my_real == double)
        return acos(dotp/(norm_n1*norm_n2));
    #else
       return acosf(dotp/(norm_n1*norm_n2));
    #endif
}

/*
list_del_pts : array of size (nb_pts x 2). For each line, indicates if a point must be deleted or not.
    - if -1 : should not be deleted.
    - otherwise : should be deleted, and the cells of this line register the two other points linked to the one deleted.
list_changed_edges : array of size (nb_pts x 3). Same thing as list_del_pts, but each line registers the two edges connected to this point,
and the face concerned.
*/
static void detect_pts_to_delete(Polygon2D* p, my_real dt, const my_real *vsx, const my_real *vsy, my_real minimal_length, my_real minimal_angle,\
                                 Array_int **list_del_pts, Array_int **list_changed_edges){
    uint64_t nb_pts = p->vertices->size;
    GrB_Index nb_faces, nb_edges;
    GrB_Index size_nz_fj, size_nz_ej, size_nz_e_ipt0;
    uint64_t i_f, i_e, i_e_fj, i_pt0;
    uint64_t i_e1, i_e0;
    uint64_t i_e0_pt0, i_e0_pt1;
    uint64_t i_e1_pt0, i_e1_pt1;
    int64_t val;
    GrB_Info infogrb;
    GrB_Vector fj, ej;
    GrB_Matrix e_ipt0;
    GxB_Iterator iterator;
    GrB_Vector nz_fj, nz_ej;
    GrB_Vector extr_vals_fj, extr_vals_ej;
    GrB_Vector nz_e_ipt0, extr_vals_e_ipt0, I_vec;
    Point2D *e1_ext1, *e1_ext2, *e2_ext1, *e2_ext2;
    my_real alpha, norm_e1, norm_e2;
    GrB_Vector inds_pts_edge0, inds_pts_edge1;
    int64_t ind_pt_del, ind_pt1, ind_pt2, curr_i;
    int64_t signed_ie0, signed_ie1, signed_i_f;
    Vector_int64 *list_inds_del1, *list_inds_del2, *list_inds_del;
    int8_t sign_f, sign_e;
    Point2D barycenter;

    GxB_Iterator_new(&iterator) ;

    *list_del_pts = alloc_with_capacity_arr_int(nb_pts, 2);
    *list_changed_edges = alloc_with_capacity_arr_int(nb_pts, 3);

    val = -1;
    for (i_f=0; i_f<nb_pts; i_f++){
        set_ijth_elem_arr_int(*list_del_pts, i_f, 0, &val);
        set_ijth_elem_arr_int(*list_del_pts, i_f, 1, &val);
        set_ijth_elem_arr_int(*list_changed_edges, i_f, 0, &val);
        set_ijth_elem_arr_int(*list_changed_edges, i_f, 1, &val);
        set_ijth_elem_arr_int(*list_changed_edges, i_f, 2, &val);
    }

    GrB_Matrix_ncols(&nb_faces, *(p->faces));
    GrB_Matrix_ncols(&nb_edges, *(p->edges));

    //rows_f = rowvals(p.faces)
    //rows_e = rowvals(p.edges)
    infogrb = GrB_Vector_new(&fj, GrB_INT8, nb_edges);
    infogrb = GrB_Vector_new(&nz_fj, GrB_UINT64, nb_edges);
    infogrb = GrB_Vector_new(&extr_vals_fj, GrB_INT8, nb_edges);
    infogrb = GrB_Vector_new(&ej, GrB_INT8, nb_pts);
    infogrb = GrB_Vector_new(&nz_ej, GrB_UINT64, 2);
    infogrb = GrB_Vector_new(&extr_vals_ej, GrB_INT8, 2);
    infogrb = GrB_Matrix_new(&e_ipt0, GrB_INT8, 1, nb_edges);
    infogrb = GrB_Vector_new(&nz_e_ipt0, GrB_UINT64, nb_edges);
    infogrb = GrB_Vector_new(&I_vec, GrB_UINT64, nb_edges);
    infogrb = GrB_Vector_new(&extr_vals_e_ipt0, GrB_INT8, nb_edges);
    infogrb = GrB_Vector_new(&inds_pts_edge0, GrB_UINT64, 2);
    infogrb = GrB_Vector_new(&inds_pts_edge1, GrB_UINT64, 2);
    list_inds_del1 = alloc_with_capacity_vec_int64(1);
    list_inds_del2 = alloc_with_capacity_vec_int64(1);
    list_inds_del = alloc_with_capacity_vec_int64(1);

    for (i_f=0 ; i_f<nb_faces ; i_f++){ //for each face
        infogrb = GrB_extract(fj, GrB_NULL, GrB_NULL, *(p->faces), GrB_ALL, 1, i_f, GrB_NULL); //Get indices of edges composing face i_f
        infogrb = GxB_Vector_extractTuples_Vector(nz_fj, extr_vals_fj, fj, GrB_NULL);
        infogrb = GrB_Vector_size(&size_nz_fj, nz_fj);

        for(i_e_fj = 0; i_e_fj < size_nz_fj; i_e_fj++){
            infogrb = GrB_Vector_extractElement(&i_e, nz_fj, i_e_fj);
            infogrb = GrB_extract(ej, GrB_NULL, GrB_NULL, *(p->edges), GrB_ALL, 1, i_e, GrB_NULL); //Get indices of points composing edge i_e
            infogrb = GxB_Vector_extractTuples_Vector(nz_ej, extr_vals_ej, ej, GrB_NULL);
            infogrb = GrB_Vector_size(&size_nz_ej, nz_ej);
            if (size_nz_ej>0){
                infogrb = GrB_Vector_extractElement(&sign_f, extr_vals_fj, i_e_fj);
                infogrb = GrB_Vector_extractElement(&sign_e, extr_vals_ej, 0);
                if (sign_f*sign_e == -1)
                    infogrb = GrB_Vector_extractElement(&i_pt0, nz_ej, 0);
                else
                    infogrb = GrB_Vector_extractElement(&i_pt0, nz_ej, 1);
                infogrb = GrB_extract(e_ipt0, GrB_NULL, GrB_NULL, *(p->edges), (GrB_Index[]){i_pt0}, 1, GrB_ALL, 1, GrB_DESC_T0); //Get indices of edges connected with point i_pt0
                infogrb = GxB_Matrix_extractTuples_Vector(I_vec, nz_e_ipt0, extr_vals_e_ipt0, e_ipt0, GrB_NULL);
                infogrb = GrB_Vector_size(&size_nz_e_ipt0, nz_e_ipt0);
                if (size_nz_e_ipt0 == 2){ //Exactly two edges connected to the point i_pt0
                    infogrb = GrB_Vector_extractElement(&i_e0, nz_e_ipt0, 0);
                    infogrb = GrB_Vector_extractElement(&i_e1, nz_e_ipt0, 1);
                    retrieve_ith_edge2D(p->vertices, p->edges, i_e0, &e1_ext1, &e1_ext2); 
                    retrieve_ith_edge2D(p->vertices, p->edges, i_e1, &e2_ext1, &e2_ext2);
                    if (((e1_ext1->x != e1_ext2->x) || ((e1_ext1->y != e1_ext2->y))) && \
                        ((e2_ext1->x != e2_ext2->x) || ((e2_ext1->y != e2_ext2->y)))){ //The edges do not degenerate to only one point
                        alpha = compute_angle_between_edges2D(e1_ext1, e1_ext2, e2_ext1, e2_ext2);
                        norm_e1 = sqrt((e1_ext1->x - e1_ext2->x)*(e1_ext1->x - e1_ext2->x) + (e1_ext1->y - e1_ext2->y)*(e1_ext1->y - e1_ext2->y));
                        norm_e2 = sqrt((e2_ext1->x - e2_ext2->x)*(e2_ext1->x - e2_ext2->x) + (e2_ext1->y - e2_ext2->y)*(e2_ext1->y - e2_ext2->y));
                        if ((norm_e1 + norm_e2 < minimal_length) || (alpha <= minimal_angle)){ //This point must be deleted!
                            //Find the two edges connected to this point
                            if (i_e == i_e0){
                                inds_pts_edge0 = nz_ej;
                                infogrb = GrB_extract(ej, GrB_NULL, GrB_NULL, *(p->edges), GrB_ALL, 1, i_e1, GrB_NULL);
                                infogrb = GxB_Vector_extractTuples_Vector(inds_pts_edge1, extr_vals_ej, ej, GrB_NULL);
                            } else { //i_e == i_e1
                                inds_pts_edge1 = nz_ej;
                                infogrb = GrB_extract(ej, GrB_NULL, GrB_NULL, *(p->edges), GrB_ALL, 1, i_e0, GrB_NULL);
                                infogrb = GxB_Vector_extractTuples_Vector(inds_pts_edge0, extr_vals_ej, ej, GrB_NULL);
                            }
                            //Find the index of the points composing these two edges
                            infogrb = GrB_Vector_extractElement(&i_e0_pt0, inds_pts_edge0, 0);
                            infogrb = GrB_Vector_extractElement(&i_e0_pt1, inds_pts_edge0, 1);
                            infogrb = GrB_Vector_extractElement(&i_e1_pt0, inds_pts_edge1, 0);
                            infogrb = GrB_Vector_extractElement(&i_e1_pt1, inds_pts_edge1, 1);

                            //The point in common in these two edges is the point to delete.
                            if ((i_e0_pt0 == i_e1_pt0) || (i_e0_pt0 == i_e1_pt1)){
                                ind_pt_del = (int64_t)i_e0_pt0;
                                ind_pt1 = (int64_t)i_e0_pt1;
                            } else {
                                ind_pt_del = (int64_t)i_e0_pt1;
                                ind_pt1 = (int64_t)i_e0_pt0;
                            }
                            if ((int64_t)i_e1_pt0 == ind_pt_del)
                                ind_pt2 = (int64_t)i_e1_pt1;
                            else
                                ind_pt2 = (int64_t)i_e1_pt0;

                            //Store the indices of the points connected to the one to delete
                            set_ijth_elem_arr_int(*list_del_pts, ind_pt_del, 0, &ind_pt1);
                            set_ijth_elem_arr_int(*list_del_pts, ind_pt_del, 1, &ind_pt2);

                            //Store the indices of the edges and face connected to the point to delete
                            signed_ie0 = (int64_t)i_e0;
                            signed_ie1 = (int64_t)i_e1;
                            signed_i_f = (int64_t)i_f;
                            set_ijth_elem_arr_int(*list_changed_edges, ind_pt_del, 0, &signed_ie0);
                            set_ijth_elem_arr_int(*list_changed_edges, ind_pt_del, 1, &signed_ie1);
                            set_ijth_elem_arr_int(*list_changed_edges, ind_pt_del, 2, &signed_i_f);
                        }
                    }
                }
            }
        }
    }

    //Correction on the lists: when several points, forming a path connected by edges, are supposed to be suppressed,
    //then only register the end points of this path. Also, register the edges in the middle to be entirely suppressed.
    i_pt0 = p->vertices->size;
    for (ind_pt_del=0; ind_pt_del<(int64_t)i_pt0 ; ind_pt_del++){ //for each point
        ind_pt1 = *get_ijth_elem_arr_int(*list_del_pts, ind_pt_del, 0);
        if ((ind_pt1 > -1) && (ind_pt1<(int64_t)i_pt0)) { //point should be deleted
            ind_pt2 = *get_ijth_elem_arr_int(*list_del_pts, ind_pt_del, 1);
            set_ith_elem_vec_int64(list_inds_del1, 0, &ind_pt_del);
            list_inds_del1->size = 1;
            set_ith_elem_vec_int64(list_inds_del2, 0, &ind_pt_del);
            list_inds_del2->size = 1;

            curr_i = ind_pt_del;
            
            signed_ie0 = *get_ijth_elem_arr_int(*list_del_pts, ind_pt1, 0);
            while (signed_ie0 > -1) {
                if (signed_ie0 == curr_i){
                    curr_i = ind_pt1;
                    ind_pt1 = *get_ijth_elem_arr_int(*list_del_pts, ind_pt1, 1);
                } else {
                    curr_i = ind_pt1;
                    ind_pt1 = signed_ie0;
                }
                if (curr_i == ind_pt_del){
                    ind_pt1 = curr_i;
                    break;
                }
                push_back_vec_int64(list_inds_del1, &curr_i);
                signed_ie0 = *get_ijth_elem_arr_int(*list_del_pts, ind_pt1, 0);
            }

            if (curr_i != ind_pt_del){ //no loop detected: we can explore the other end.
                curr_i = ind_pt_del;
                signed_ie0 = *get_ijth_elem_arr_int(*list_del_pts, ind_pt2, 0);
                while (signed_ie0 > -1){
                    if (signed_ie0 == curr_i){
                        curr_i = ind_pt2;
                        ind_pt2 = *get_ijth_elem_arr_int(*list_del_pts, ind_pt2, 1);
                    } else {
                        curr_i = ind_pt2;
                        ind_pt2 = signed_ie0;
                    }
                    if (curr_i == ind_pt_del){
                        ind_pt2 = curr_i;
                        break;
                    }
                    push_back_vec_int64(list_inds_del2, &curr_i);
                    signed_ie0 = *get_ijth_elem_arr_int(*list_del_pts, ind_pt2, 0);
                }
            }

            if ((list_inds_del1->size > 1) || (list_inds_del2->size > 1)){ //There are several points linked together to be deleted!
                if ((ind_pt1 == ind_pt_del) && (list_inds_del1->size < 5)){ //we have looped, and it is a triangle or a square!
                    //list_inds_del = unique([list_inds_del1; list_inds_del2]);
                    set_ith_elem_vec_int64(list_inds_del, 0, get_ith_elem_vec_int64(list_inds_del1, 0));
                    list_inds_del->size = 1;
                    for(i_f = 1; i_f<list_inds_del1->size; i_f++){
                        push_back_unique_vec_int64(list_inds_del, get_ith_elem_vec_int64(list_inds_del1, i_f));
                    }
                    for(i_f = 0; i_f<list_inds_del2->size; i_f++){
                        push_back_unique_vec_int64(list_inds_del, get_ith_elem_vec_int64(list_inds_del2, i_f));
                    }

                    barycenter = (Point2D){0., 0.};
                    signed_i_f = (int64_t)(list_inds_del->size + 1);
                    for (i_f = 0; i_f < list_inds_del->size; i_f++){
                        i_e = *get_ith_elem_vec_int64(list_inds_del, i_f);
                        e1_ext1 = get_ith_elem_vec_pts2D(p->vertices, i_e); //p.vertices[list_inds_del1[i_f]]
                        barycenter.x += e1_ext1->x + dt*vsx[i_e];
                        barycenter.y += e1_ext1->y + dt*vsy[i_e];
                        set_ijth_elem_arr_int(*list_del_pts, i_e, 0, &signed_i_f);
                        set_ijth_elem_arr_int(*list_del_pts, i_e, 1, &signed_i_f);
                    }
                    barycenter.x /= list_inds_del->size;
                    barycenter.y /= list_inds_del->size;
                    push_back_vec_pts2D(p->vertices, &barycenter);

                    signed_i_f = *get_ijth_elem_arr_int(*list_changed_edges, *get_ith_elem_vec_int64(list_inds_del, 0), 0);
                    for(i_f = 0; i_f < list_inds_del->size; i_f++){
                        i_e = *get_ith_elem_vec_int64(list_inds_del, i_f);
                        set_ijth_elem_arr_int(*list_changed_edges, i_e, 0, &signed_i_f);
                        set_ijth_elem_arr_int(*list_changed_edges, i_e, 1, &signed_i_f);
                    }
                } else {
                    if (list_inds_del1->size > 1){
                        i_f = 1;
                        while (i_f<list_inds_del1->size){
                            signed_i_f = *get_ith_elem_vec_int64(list_inds_del1, i_f);
                            i_e = (uint64_t)signed_i_f;
                            signed_i_f = -1;
                            set_ijth_elem_arr_int(*list_del_pts, i_e, 0, &signed_i_f);//Actually don't suppress that point...
                            set_ijth_elem_arr_int(*list_del_pts, i_e, 1, &signed_i_f);
                            set_ijth_elem_arr_int(*list_changed_edges, i_e, 0, &signed_i_f);//Actually don't suppress that point...
                            set_ijth_elem_arr_int(*list_changed_edges, i_e, 1, &signed_i_f);
                            set_ijth_elem_arr_int(*list_changed_edges, i_e, 2, &signed_i_f);
                            i_f += 2;
                        }
                    }
                    if (list_inds_del2->size > 1){
                        i_f = 1;
                        while (i_f<list_inds_del2->size){
                            signed_i_f = *get_ith_elem_vec_int64(list_inds_del2, i_f);
                            i_e = (uint64_t)signed_i_f;
                            signed_i_f = -1;
                            set_ijth_elem_arr_int(*list_del_pts, i_e, 0, &signed_i_f);//Actually don't suppress that point...
                            set_ijth_elem_arr_int(*list_del_pts, i_e, 1, &signed_i_f);
                            set_ijth_elem_arr_int(*list_changed_edges, i_e, 0, &signed_i_f);//Actually don't suppress that point...
                            set_ijth_elem_arr_int(*list_changed_edges, i_e, 1, &signed_i_f);
                            set_ijth_elem_arr_int(*list_changed_edges, i_e, 2, &signed_i_f);
                            i_f += 2;
                        }
                    }
                }
            }
        }
    }
    infogrb = GrB_free(&fj);
    infogrb = GrB_free(&nz_fj);
    infogrb = GrB_free(&extr_vals_fj);
    infogrb = GrB_free(&ej);
    infogrb = GrB_free(&nz_ej);
    infogrb = GrB_free(&extr_vals_ej);
    infogrb = GrB_free(&e_ipt0);
    infogrb = GrB_free(&nz_e_ipt0);
    infogrb = GrB_free(&extr_vals_e_ipt0);
    infogrb = GrB_free(&inds_pts_edge0);
    infogrb = GrB_free(&inds_pts_edge1);
    infogrb = GrB_free(&I_vec);
    dealloc_vec_int64(list_inds_del1);
    dealloc_vec_int64(list_inds_del2);
    dealloc_vec_int64(list_inds_del);
    GxB_Iterator_free(&iterator);
}


static Vector_points2D* build_vertices_tnp1(const Polygon2D* p, const my_real* vsx, const my_real* vsy, uint64_t size_vs, my_real dt, Array_int* list_del_pts){
    Vector_points2D* vertices_tnp1;
    uint64_t i;
    int64_t ind_pt1, ind_pt2;
    Point2D pt;
    Point2D pt3D, pt3D_2;

    vertices_tnp1 = alloc_with_capacity_vec_pts2D(p->vertices->size);
    if (!list_del_pts){
        for (i = 0; i<p->vertices->size; i++){
            pt = *get_ith_elem_vec_pts2D(p->vertices, i);
            pt3D = (Point2D){pt.x + dt*vsx[i], pt.y + dt*vsy[i]};
            set_ith_elem_vec_pts2D(vertices_tnp1, i, &pt3D);
        }
    } else {
        for (i = 0; i<size_vs; i++){
            ind_pt1 = *get_ijth_elem_arr_int(list_del_pts, i, 0);
            if (ind_pt1 == -1){
                pt = *get_ith_elem_vec_pts2D(p->vertices, i);
                pt3D = (Point2D){pt.x + dt*vsx[i], pt.y + dt*vsy[i]};
                set_ith_elem_vec_pts2D(vertices_tnp1, i, &pt3D);
            }
        }

        for (i = 0; i<size_vs; i++){ //correct when point must be suppressed.
            ind_pt1 = *get_ijth_elem_arr_int(list_del_pts, i, 0);
            if ((ind_pt1 > -1) && (ind_pt1<(int64_t)size_vs)){
                ind_pt2 = *get_ijth_elem_arr_int(list_del_pts, i, 1);

                pt3D = *get_ith_elem_vec_pts2D(vertices_tnp1, ind_pt1);
                pt3D_2 = *get_ith_elem_vec_pts2D(vertices_tnp1, ind_pt2);

                pt3D.x = 0.5*(pt3D.x + pt3D_2.x);
                pt3D.y = 0.5*(pt3D.y + pt3D_2.y);

                set_ith_elem_vec_pts2D(vertices_tnp1, i, &pt3D);
            } else if (ind_pt1>=(int64_t)size_vs){
                pt = *get_ith_elem_vec_pts2D(p->vertices, ind_pt1);
                set_ith_elem_vec_pts2D(vertices_tnp1, i, &pt);
            }
        }
    }
    
    return vertices_tnp1;
}

static Polyhedron3D* build_space_time_cell_given_tnp1_vertices(const Polygon2D* fn, const Vector_points2D* vertices_tnp1, my_real dt, bool split){
    uint64_t i;
    int8_t fne_ptind;
    long int val;
    Point2D *pt2D;
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
    //for (i = nb_pts; i<2*nb_pts; i++){
    //    pt2D = get_ith_elem_vec_pts2D(fn->vertices, i-nb_pts);
    //    v = get_ith_elem_vec_pts2D(vs, i-nb_pts);
    //    *pt3D = (Point3D){pt2D->x + dt*v->x, pt2D->y + dt*v->y, dt};
    //    set_ith_elem_vec_pts3D(new_vertices, i, pt3D);
    //}
    for(i = 0; i<nb_pts; i++){
        pt2D = get_ith_elem_vec_pts2D(vertices_tnp1, i);
        *pt3D = (Point3D){pt2D->x, pt2D->y, dt};
        set_ith_elem_vec_pts3D(new_vertices, nb_pts + i, pt3D);
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

    //pressure_face = repeat([0.0], size(new_faces, 2)) //TODO : report that
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

            //pressure_face[curr_face] = fn.pressure_edge[i] //TODO Report this
            //pressure_face[curr_face+1] = fn.pressure_edge[i] //TODO Report this


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

            //pressure_face[curr_face] = fn.pressure_edge[i] //TODO Report this

            GrB_Matrix_extractElement(&fne_ptind, *(fn->faces), i, 0);
            GrB_Matrix_setElement(*new_volume, fne_ptind, curr_face, 0);

            curr_edge += 1;
            curr_face += 1;
        }
    }


    free(pt3D);
    GrB_Matrix_free(&emptyNE);
    GrB_Matrix_free(&emptySW);
    GrB_Matrix_free(&emptyE_N);
    GrB_Matrix_free(&emptyE_S);
    GrB_Vector_free(&ed_i);
    GrB_Vector_free(&nz_ei);
    GrB_Vector_free(&val_nz_ei);

    return new_Polyhedron3D_vefvs(new_vertices, new_edges, new_faces, new_volume, status_faces);
}

/// @brief Using a polygon `fn` defined at time t^n and a set of vectors `vs` defining the translation of each point of `fn` at time t^{n+1} = t^n+`dt`,
///        this function builds a 3D polyhedron (2D + time) connecting `fn` and the polygon at time t^{n+1}. 
/// @details In the resulting polyhedron `status_faces`, 1 denotes the face created at t^n, 
///          2 the face created at t^{n+1}, and i+2 the face created in time-space from edge i.
/// @param fn [IN] polygon at time t^n
/// @param vs* [IN] Vector translating points of `fn` over time `dt` to form the polygon at time t^{n+1}
/// @param dt [IN] time-step
/// @param split [IN] flag to triangulate (or not) faces in time linking the edges at time t^n and t^{n+1}.
Polyhedron3D* build_space2D_time_cell(const Polygon2D *fn, const my_real *vsx, const my_real* vsy, uint64_t size_vs, const my_real dt, bool split, Array_int *list_del_pts){
    Vector_points2D* vertices_tnp1 = build_vertices_tnp1(fn, vsx, vsy, size_vs, dt, list_del_pts);
    return build_space_time_cell_given_tnp1_vertices(fn, vertices_tnp1, dt, split);
}

static void split_edge(Polygon2D* p, uint64_t i_e, Point2D *new_pt){
    GrB_Index nb_edge;
    GrB_Info infogrb;
    GrB_Vector ej, extr_vals_ej, inds_pts_edge;
    GrB_Index i_e0_pt0, i_e0_pt1;
    uint64_t ind_new_pt;
    int8_t val;
    GrB_Matrix new_line, new_edge;
    GrB_Index nb_faces;

    GrB_Matrix_ncols(&nb_edge, *(p->edges));
    GrB_Matrix_ncols(&nb_faces, *(p->faces));
    infogrb = GrB_Vector_new(&ej, GrB_INT8, 2);
    infogrb = GrB_Vector_new(&extr_vals_ej, GrB_INT8, 2);
    infogrb = GrB_Vector_new(&inds_pts_edge, GrB_UINT64, 2);

    //Extract the column i_e of p.edges, and take the two first non-zero values.
    infogrb = GrB_extract(ej, GrB_NULL, GrB_NULL, *(p->edges), GrB_ALL, 1, i_e, GrB_NULL);
    infogrb = GxB_Vector_extractTuples_Vector(inds_pts_edge, extr_vals_ej, ej, GrB_NULL);
    infogrb = GrB_Vector_extractElement(&i_e0_pt0, inds_pts_edge, 0);
    infogrb = GrB_Vector_extractElement(&i_e0_pt1, inds_pts_edge, 1);
    infogrb = GrB_Vector_extractElement(&val, extr_vals_ej, 0);

    //pts_indices = rowvals(p.edges)
    //inds_pts_edge = pts_indices[nzrange(p.edges, i_e)]

    push_back_vec_pts2D(p->vertices, new_pt);
    ind_new_pt = p->vertices->size;

    infogrb = GrB_Matrix_setElement(*(p->edges), 0, i_e0_pt0, i_e);

    infogrb = GrB_Matrix_new(&new_line, GrB_INT8, 1, nb_edge);
    infogrb = GrB_Matrix_new(&new_edge, GrB_INT8, ind_new_pt, 1);
    infogrb = GrB_Matrix_setElement(new_line, val, 0, i_e);
    infogrb = GrB_Matrix_setElement(new_edge, val, i_e0_pt0, 0);
    infogrb = GrB_Matrix_setElement(new_edge, -val, ind_new_pt, 0);
    
    infogrb = GxB_Matrix_concat(*(p->edges), (GrB_Matrix[]){*(p->edges), new_line}, 2, 1, GrB_NULL); //edges = [edges ; new_line]
    infogrb = GxB_Matrix_concat(*(p->edges), (GrB_Matrix[]){*(p->edges), new_edge}, 1, 2, GrB_NULL); //edges = [edges   new_edge]

    infogrb = GrB_Matrix_resize(new_line, 1, nb_faces);
    infogrb = GrB_extract(new_line, GrB_NULL, GrB_NULL, *(p->faces), (GrB_Index[]){i_e}, 1, GrB_ALL, 1,GrB_NULL); 
    infogrb = GxB_Matrix_concat(*(p->faces), (GrB_Matrix[]){*(p->faces), new_line}, 2, 1, GrB_NULL); //faces = [faces ; faces[i_e,:]]

    //push!(p.pressure_edge, p.pressure_edge[i_e]) //TODO Report this
    //push!(p.status_edge, p.status_edge[i_e]) //TODO Report this

    dropzeros(p->edges);
    dropzeros(p->faces);

    GrB_Vector_free(&ej);
    GrB_Vector_free(&extr_vals_ej);
    GrB_Vector_free(&inds_pts_edge);
    GrB_Matrix_free(&new_line);
    GrB_Matrix_free(&new_edge);
}

static void fuse_points(Polygon2D *p, int64_t ind_e0_fused, int64_t ind_e1_fused, int64_t ind_face, bool do_dropzeros){
    GrB_Info infogrb;
    GrB_Vector ej, extr_vals_ej, inds_pts_edge0, inds_pts_edge1;
    GrB_Index i_e0_pt0, i_e0_pt1, i_e1_pt0, i_e1_pt1;
    uint64_t ind_pt_del, ind_pt1, ind_pt2;
    int8_t val;

    infogrb = GrB_Vector_new(&ej, GrB_INT8, 2);
    infogrb = GrB_Vector_new(&extr_vals_ej, GrB_INT8, 2);
    infogrb = GrB_Vector_new(&inds_pts_edge0, GrB_UINT64, 2);
    infogrb = GrB_Vector_new(&inds_pts_edge1, GrB_UINT64, 2);

    if (ind_e1_fused != ind_e0_fused) { //one edge to suppress, the other one is modified.
        infogrb = GrB_extract(ej, GrB_NULL, GrB_NULL, *(p->edges), GrB_ALL, 1, ind_e0_fused, GrB_NULL);
        infogrb = GxB_Vector_extractTuples_Vector(inds_pts_edge0, extr_vals_ej, ej, GrB_NULL);
        infogrb = GrB_extract(ej, GrB_NULL, GrB_NULL, *(p->edges), GrB_ALL, 1, ind_e1_fused, GrB_NULL);
        infogrb = GxB_Vector_extractTuples_Vector(inds_pts_edge1, extr_vals_ej, ej, GrB_NULL);
        //Find the index of the points composing these two edges
        infogrb = GrB_Vector_extractElement(&i_e0_pt0, inds_pts_edge0, 0);
        infogrb = GrB_Vector_extractElement(&i_e0_pt1, inds_pts_edge0, 1);
        infogrb = GrB_Vector_extractElement(&i_e1_pt0, inds_pts_edge1, 0);
        infogrb = GrB_Vector_extractElement(&i_e1_pt1, inds_pts_edge1, 1);

        //The point in common in these two edges is the point to delete.
        if ((i_e0_pt0 == i_e1_pt0) || (i_e0_pt0 == i_e1_pt1)){
            ind_pt_del = i_e0_pt0;
            ind_pt1 = i_e0_pt1;
        } else {
            ind_pt_del = i_e0_pt1;
            ind_pt1 = i_e0_pt0;
        }
        if (i_e1_pt0 == ind_pt_del)
            ind_pt2 = i_e1_pt1;
        else
            ind_pt2 = i_e1_pt0;
        
        infogrb = GrB_Matrix_extractElement(&val, *(p->edges), ind_pt_del, ind_e0_fused);
        infogrb = GrB_Matrix_setElement(*(p->edges), val, ind_pt2, ind_e0_fused);
        infogrb = GrB_Matrix_setElement(*(p->edges), 0, ind_pt_del, ind_e0_fused);
        infogrb = GrB_Matrix_setElement(*(p->edges), 0, ind_pt2, ind_e1_fused);
        infogrb = GrB_Matrix_setElement(*(p->edges), 0, ind_pt_del, ind_e1_fused);

        infogrb = GrB_Matrix_setElement(*(p->faces), 0, ind_e1_fused, ind_face);
    } else { //only one edge: means it must be suppressed.
        infogrb = GrB_extract(ej, GrB_NULL, GrB_NULL, *(p->edges), GrB_ALL, 1, ind_e0_fused, GrB_NULL);
        infogrb = GxB_Vector_extractTuples_Vector(inds_pts_edge0, extr_vals_ej, ej, GrB_NULL);
        infogrb = GrB_Vector_extractElement(&i_e0_pt0, inds_pts_edge0, 0);
        infogrb = GrB_Vector_extractElement(&i_e0_pt1, inds_pts_edge0, 1);

        //inds_pts_edge = pts_indices[nzrange(p.edges, ind_e1_fused)]
        infogrb = GrB_Matrix_setElement(*(p->edges), 0, i_e0_pt0, ind_e1_fused);
        infogrb = GrB_Matrix_setElement(*(p->edges), 0, i_e0_pt1, ind_e1_fused);

        infogrb = GrB_Matrix_setElement(*(p->faces), 0, ind_e1_fused, ind_face);
    }

    if (do_dropzeros){
        dropzeros(p->edges);
        dropzeros(p->faces);
    }

    GrB_Vector_free(&ej);
    GrB_Vector_free(&extr_vals_ej);
    GrB_Vector_free(&inds_pts_edge0);
    GrB_Vector_free(&inds_pts_edge1);
}

static void refine_interface(Polygon2D *p, my_real maximal_length){
    uint64_t i, nb_edges;
    Point2D *ext1, *ext2;
    my_real norm_e;

    GrB_Matrix_ncols(&nb_edges, *(p->edges));
    for(i=0; i<nb_edges; i++){
        retrieve_ith_edge2D(p->vertices, p->edges, i, &ext1, &ext2);
        if ((ext1->x != ext2->x) || (ext1->y != ext2->y)){
            norm_e = sqrt((ext1->x-ext2->x)*(ext1->x-ext2->x) + (ext1->y-ext2->y)*(ext1->y-ext2->y));
            if (norm_e > maximal_length){
                ext1->x = 0.5*(ext1->x + ext2->x);
                ext1->y = 0.5*(ext1->y + ext2->y);
                split_edge(p, i, ext1);
            }
        }
    }
}

static void coarsen_interface(Polygon2D *p, Array_int *list_changed_edges){
    uint64_t i;
    int64_t ind_e1_fused, ind_e2_fused, ind_f;

    for (i=0; i<p->vertices->size; i++){
        ind_e1_fused = *get_ijth_elem_arr_int(list_changed_edges, i, 0);
        if (ind_e1_fused > -1){
            ind_e2_fused = *get_ijth_elem_arr_int(list_changed_edges, i, 1);
            ind_f = *get_ijth_elem_arr_int(list_changed_edges, i, 2);
            fuse_points(p, ind_e1_fused, ind_e2_fused, ind_f, false);
        }
    }
    dropzeros(p->edges);
    dropzeros(p->faces);
}

void update_solid(Polygon2D **solid, Polyhedron3D** solid3D, const my_real* vec_move_solidx, const my_real* vec_move_solidy, my_real dt, \
                    my_real minimal_length, my_real maximal_length, my_real minimal_angle){
    Array_int *list_del_pts,  *list_changed_edges;
    Vector_points2D* vertices_tnp1;
    Polygon2D* solid_new;
    uint64_t i, nb_faces, ind_pt, j_f, j;
    Vector_int64* ind_kept_pts;
    Vector_points2D* new_vertices;
    GrB_Matrix new_edges, new_faces;
    GrB_Vector grb_ind_kept_pts, justone;
    Vector_int *new_status_edge;
    GrB_Info infogrb;
    GrB_Vector fj, extr_vals_fj, edge_indices;
    GrB_Vector ej, extr_vals_ej, pt_indices;
    GrB_Index nb_edges, size_edge_indices, size_pt_indices, val;
    GrB_Index ncols_new_edges, ell;
    int64_t nrows_new_edges;
    const long int zero = 0;
    uint64_t size_vs = (*solid)->vertices->size;
    
    detect_pts_to_delete(*solid, dt, vec_move_solidx, vec_move_solidy, minimal_length, minimal_angle, &list_del_pts, &list_changed_edges);

    //Create 3D space-time solid interface: change this when self-intersection occurs!
    vertices_tnp1 = build_vertices_tnp1(*solid, vec_move_solidx, vec_move_solidy, size_vs, dt, list_del_pts);
    *solid3D = build_space_time_cell_given_tnp1_vertices(*solid, vertices_tnp1, dt, true);
    for(i=0; i<(*solid3D)->status_face->size; i++){
        if ((*solid3D)->status_face->data[i] > 2)
            (*solid3D)->status_face->data[i] = -(*solid3D)->status_face->data[i]; //Face status numbered negatively for solid space-time reconstructions
    }

    solid_new = new_Polygon2D();
    copy_Polygon2D(*solid, solid_new);
    copy_vec_pts2D(vertices_tnp1, solid_new->vertices);
    if (minimal_length > 0)
        coarsen_interface(solid_new, list_changed_edges);
    if (maximal_length<INFINITY)
        refine_interface(solid_new, maximal_length);
    
    //Clean the result if some points were suppressed
    GrB_Matrix_ncols(&nb_faces, *(solid_new->faces));
    GrB_Matrix_ncols(&nb_edges, *(solid_new->edges));
    ind_kept_pts = alloc_with_capacity_vec_int64(1);
    new_vertices = alloc_with_capacity_vec_pts2D(1);
    GrB_Matrix_new(&new_edges, GrB_INT8, 1, 1);
    GrB_Vector_new(&grb_ind_kept_pts, GrB_UINT64, 1);
    GrB_Matrix_new(&new_faces, GrB_INT8, 1, 1);
    GrB_Vector_new(&justone, GrB_UINT64, 1);
    new_status_edge = alloc_with_capacity_vec_int(1);
    GrB_Vector_setElement(justone, 0, 0);
    GrB_Vector_new(&fj, GrB_INT8, nb_edges);
    GrB_Vector_new(&edge_indices, GrB_UINT64, nb_edges);
    GrB_Vector_new(&extr_vals_fj, GrB_INT8, nb_edges);
    GrB_Vector_new(&ej, GrB_INT8, nb_edges);
    GrB_Vector_new(&pt_indices, GrB_UINT64, nb_edges);
    GrB_Vector_new(&extr_vals_ej, GrB_INT8, nb_edges);

    for (i=0; i<nb_faces; i++){
        //krange = nzrange(solid_new.faces, i)
        //edge_indices = rowvals(solid_new.faces)[krange]
        infogrb = GrB_extract(fj, GrB_NULL, GrB_NULL, *(solid_new->faces), GrB_ALL, 1, i, GrB_NULL); //Get indices of edges composing face i
        infogrb = GxB_Vector_extractTuples_Vector(edge_indices, extr_vals_fj, fj, GrB_NULL);
        infogrb = GrB_Vector_size(&size_edge_indices, edge_indices);

        ind_kept_pts->size = 0;
        for(j_f=0; j_f<size_edge_indices; j_f++){
            infogrb = GrB_Vector_extractElement(&j, edge_indices, j_f);
            infogrb = GrB_extract(ej, GrB_NULL, GrB_NULL, *(solid_new->edges), GrB_ALL, 1, j, GrB_NULL); //Get indices of points composing edge j
            infogrb = GxB_Vector_extractTuples_Vector(pt_indices, extr_vals_ej, ej, GrB_NULL);
            infogrb = GrB_Vector_size(&size_pt_indices, pt_indices);
            if(size_pt_indices>0){
                infogrb = GrB_Vector_extractElement(&val, pt_indices, 0);
                if(ind_kept_pts->size == 0)
                    ind_kept_pts->data[0] = val;
                else
                    push_back_unique_vec_int64(ind_kept_pts, &val);
                
                infogrb = GrB_Vector_extractElement(&val, pt_indices, 1);
                push_back_unique_vec_int64(ind_kept_pts, &val);
            }
        }
        sort_vec_int64(ind_kept_pts);
        if (ind_kept_pts->size > 0){
            ind_pt = *get_ith_elem_vec_int64(ind_kept_pts, 0);
            set_ith_elem_vec_pts2D(new_vertices, 0, get_ith_elem_vec_pts2D(solid_new->vertices, ind_pt));
            for (j_f=1; j_f<ind_kept_pts->size; j_f++){
                ind_pt = *get_ith_elem_vec_int64(ind_kept_pts, j_f);
                set_ith_elem_vec_pts2D(new_vertices, j_f, get_ith_elem_vec_pts2D(solid_new->vertices, ind_pt));
            }
        }
        
        //new_edges = solid_new.edges[ind_kept_pts, edge_indices]
        GrB_reduce(&ncols_new_edges, GrB_NULL, GrB_MAX_MONOID_UINT64, edge_indices, GrB_NULL);
        if (ind_kept_pts->size>1){
            nrows_new_edges = *get_ith_elem_vec_int64(ind_kept_pts, ind_kept_pts->size-1);
            GrB_Vector_resize(grb_ind_kept_pts, ind_kept_pts->size);
            for(ell=0; ell<ind_kept_pts->size; ell++){
                GrB_Vector_setElement(grb_ind_kept_pts, *get_ith_elem_vec_int64(ind_kept_pts, ell), ell);
            }
            GrB_Matrix_resize(new_edges, nrows_new_edges, ncols_new_edges);
            GrB_extract(new_edges, GrB_NULL, GrB_NULL, *(solid_new->edges), grb_ind_kept_pts, edge_indices, GrB_NULL);

            //new_faces = solid_new.faces[edge_indices, i]
            GrB_Matrix_resize(new_faces, nrows_new_edges, 1);
            GrB_extract(new_faces, GrB_NULL, GrB_NULL, *(solid_new->faces), edge_indices, justone, GrB_NULL);
        

            //new_status_edge = zeros(nrows_new_edges)
            if(new_status_edge->size <= (uint64_t)nrows_new_edges){
                for (j_f=new_status_edge->size; j_f<(uint64_t)nrows_new_edges; j_f++){
                    set_ith_elem_vec_int(new_status_edge, j_f, &zero);
                }
            }
            new_status_edge->size = nrows_new_edges;
            //new_pressure_edge = zeros(Int, size(new_edges, 2))
            if (i == 0){ //First face created
                *solid = new_Polygon2D_vefs(new_vertices, &new_edges, &new_faces, new_status_edge);
            } else { //Faces already exist
                *solid = fuse_polygons(*solid, new_Polygon2D_vefs(new_vertices, &new_edges, &new_faces, new_status_edge));
            }
        } else {
            if (i == 0){ //First face created
                GrB_Matrix_resize(new_edges, 0, 0);
                GrB_Matrix_resize(new_faces, 0, 0);
                new_status_edge->size = 0;
                *solid = new_Polygon2D_vefs(new_vertices, &new_edges, &new_faces, new_status_edge);
                //*solid = new_Polygon2D();
            }
        }
    }

    dealloc_vec_pts2D(vertices_tnp1);
    dealloc_Polygon2D(solid_new);
    dealloc_vec_int64(ind_kept_pts);
    dealloc_vec_int(new_status_edge);
    dealloc_vec_pts2D(new_vertices);
    GrB_free(&new_edges);
    GrB_free(&new_faces);
    GrB_free(&grb_ind_kept_pts);
    GrB_free(&justone);
    GrB_free(&fj);
    GrB_free(&extr_vals_fj);
    GrB_free(&edge_indices);
    GrB_free(&ej);
    GrB_free(&extr_vals_ej);
    GrB_free(&pt_indices);
}
