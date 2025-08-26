#include "Polyhedron3D.h"

Polyhedron3D* new_Polyhedron3D(){
    GrB_Info infogrb;
    Polyhedron3D* p = (Polyhedron3D*) malloc(sizeof(Polyhedron3D));

    p->vertices = alloc_empty_vec_pts3D();
    p->edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    p->faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    p->volumes = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    infogrb = GrB_Matrix_new(p->edges, GrB_INT8, 0,0);
    infogrb = GrB_Matrix_new(p->faces, GrB_INT8, 0,0);
    infogrb = GrB_Matrix_new(p->volumes, GrB_INT8, 0,0);
    p->status_face = alloc_empty_vec_int();
    
    return p;
}

Polyhedron3D* new_Polyhedron3D_vefs(Vector_points3D* vertices, GrB_Matrix* edges, GrB_Matrix* faces, Vector_int* status_face){
    GrB_Info infogrb;
    GrB_Index nb_rows;
    Polyhedron3D* p = (Polyhedron3D*) malloc(sizeof(Polyhedron3D));

    p->vertices = alloc_empty_vec_pts3D();
    copy_vec_pts3D(vertices, p->vertices);
    p->edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix_dup(p->edges, *edges);
    p->faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix_dup(p->faces, *faces);

    p->volumes = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    infogrb = GrB_Matrix_ncols(&nb_rows, *(p->faces));
    infogrb = GrB_Matrix_new(p->volumes, GrB_INT8, nb_rows, 1);

    p->status_face = alloc_empty_vec_int();
    copy_vec_int(status_face, p->status_face);
    
    return p;
}

Polyhedron3D* new_Polyhedron3D_vefvs(Vector_points3D* vertices, GrB_Matrix* edges, GrB_Matrix* faces, GrB_Matrix* volumes, Vector_int* status_face){
    Polyhedron3D* p = (Polyhedron3D*) malloc(sizeof(Polyhedron3D));

    p->vertices = alloc_empty_vec_pts3D();
    copy_vec_pts3D(vertices, p->vertices);
    p->edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix_dup(p->edges, *edges);
    p->faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix_dup(p->faces, *faces);
    p->volumes = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix_dup(p->volumes, *volumes);
    p->status_face = alloc_empty_vec_int();
    copy_vec_int(status_face, p->status_face);
    
    return p;
}

Polyhedron3D* Polyhedron3D_from_vertices(const my_real* x_v, unsigned long int n_x, const my_real* y_v, unsigned long int n_y, const my_real* z_v, unsigned long int n_z){
    GrB_Info infogrb;
    GrB_Index nb_pts, nb_edges, nb_faces, nb_volumes;
    uint64_t curr_face, curr_pt, curr_edge, curr_vol, iy, ix, iz;
    uint64_t ind_pt_SW, ind_pt_SE, ind_pt_NW;
    uint64_t ind_pt_SW_up;//, ind_pt_SE_up, ind_pt_NW_up;
    uint64_t ind_e_W, ind_e_S, ind_e_E, ind_e_N;
    uint64_t ind_e_W_up, ind_e_S_up;
    uint64_t ind_e_SW_up, ind_e_SE_up, ind_e_NW_up;
    uint64_t face_W, face_B, face_S, face_E, face_N, face_U;
    my_real yS, xW, zB;
    Point3D ptSWB;
    Vector_points3D* vertices;
    Vector_int* status_face;
    GrB_Matrix* edges = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix* faces = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));
    GrB_Matrix* volumes = (GrB_Matrix*) malloc(sizeof(GrB_Matrix));

    if (n_x<2 || n_y<2 || n_z<2)
        return new_Polyhedron3D();

    nb_pts = n_x*n_y*n_z;
    nb_edges = ((3*(n_x-1) + 2)*(n_y-1) + 2*(n_y-1) + 1)*(n_z-1) + (2*(n_x-1) + 1)*(n_y-1) + (n_y-1);
    nb_faces = ((3*(n_x-1) + 1)*(n_y-1) + (n_y-1))*(n_z-1) + (n_x-1)*(n_y-1);
    nb_volumes = (n_x-1)*(n_y-1)*(n_z-1);

    vertices = alloc_with_capacity_vec_pts3D(nb_pts);
    infogrb = GrB_Matrix_new(edges, GrB_INT8, nb_pts, nb_edges);
    infogrb = GrB_Matrix_new(faces, GrB_INT8, nb_edges, nb_faces);
    infogrb = GrB_Matrix_new(volumes, GrB_INT8, nb_faces, nb_volumes);

    ptSWB = (Point3D){x_v[0], y_v[0], z_v[0]};
    push_back_vec_pts3D(vertices, &ptSWB);
    curr_face = 0;
    curr_pt = 0;
    curr_edge = 0;
    curr_vol = 0;
    for (iz = 0; iz < n_z-2; iz++){
        zB = z_v[iz];
        for (iy = 0; iy < n_y-2; iy++){
            yS = y_v[iy];
            for (ix = 0; ix < n_x-1; ix++){
                xW = x_v[ix];

                ind_pt_SW = curr_pt;
                ind_pt_SE = curr_pt + 1;
                ind_pt_NW = curr_pt + n_x;
                ind_pt_SW_up = curr_pt + n_x*n_y;

                ind_e_W = curr_edge;
                ind_e_SW_up = curr_edge + 1;
                ind_e_S = curr_edge + 2;
                ind_e_E = curr_edge + 3;
                ind_e_SE_up = curr_edge + 4;
                ind_e_NW_up = curr_edge + (3*(n_x-1) + 3);
                ind_e_N = curr_edge + (3*(n_x-1) + 4);
                ind_e_W_up = ind_e_W + (3*(n_x-1)+2)*(n_y-1) +  + 2*(n_y-1) + 1;
                ind_e_S_up = ind_e_S + (3*(n_x-1)+2)*(n_y-1) +  + 2*(n_y-1) + 1;

                face_W = curr_face;
                face_B = curr_face + 1;
                face_S = curr_face + 2;
                face_E = curr_face + 3;
                face_N = face_S + (3*(n_x-1) + 1);
                face_U = face_B + (3*(n_x-1) + 1)*(n_y-1) + (n_y-1);

                ptSWB = (Point3D){xW,yS,zB};

                set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

                infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_W    );
                infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW   , ind_e_W    );
                infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_SW_up);
                infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, ind_e_SW_up);
                infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_S    );
                infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE   , ind_e_S    );

                infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SW_up, face_W);
                infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_W    , face_W);
                infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_NW_up, face_W);
                infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W_up , face_W);

                infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S, face_B);
                infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_E, face_B);
                infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_N, face_B);
                infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W, face_B);

                infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_S    , face_S);
                infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SE_up, face_S);
                infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S_up , face_S);
                infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_SW_up, face_S);

                infogrb = GrB_Matrix_setElement(*volumes,  1, face_E, curr_vol);
                infogrb = GrB_Matrix_setElement(*volumes, -1, face_B, curr_vol);
                infogrb = GrB_Matrix_setElement(*volumes, -1, face_S, curr_vol);
                infogrb = GrB_Matrix_setElement(*volumes, -1, face_W, curr_vol);
                infogrb = GrB_Matrix_setElement(*volumes,  1, face_N, curr_vol);
                infogrb = GrB_Matrix_setElement(*volumes,  1, face_U, curr_vol);

                curr_pt += 1;
                curr_edge += 3;
                curr_face += 3;
                curr_vol += 1;
            }
            //Treat east-most point and edge
            xW = x_v[n_x-1];
            ptSWB = (Point3D){xW,yS,zB};
            ind_pt_SW = curr_pt;
            ind_pt_NW = curr_pt + n_x;
            ind_pt_SW_up = curr_pt + n_x*n_y;

            ind_e_W = curr_edge;
            ind_e_SW_up = curr_edge + 1;
            ind_e_NW_up = curr_edge +  + (3*n_x);
            ind_e_W_up = ind_e_W + (3*(n_x-1)+2)*(n_y-1) +  + 2*(n_y-1) + 1;

            set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_W    );
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW   , ind_e_W    );
            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_SW_up);
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, ind_e_SW_up);

            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_W    , curr_face);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SW_up, curr_face);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_NW_up, curr_face);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W_up , curr_face);

            curr_pt += 1;
            curr_edge += 2;
            curr_face += 1;
        }

        //Treat north-most cells
        iy = n_y-2;
        yS = y_v[iy];
        ind_e_NW_up = curr_edge + (3*(n_x-1));
        face_N = curr_face + 3*(n_x-1);
        for (ix = 0; ix < n_x-1; ix++){
            ind_pt_SW = curr_pt;
            ind_pt_SE = curr_pt + 1;
            ind_pt_NW = curr_pt + n_x;
            ind_pt_SW_up = curr_pt + n_x*n_y;

            ind_e_W = curr_edge;
            ind_e_SW_up = curr_edge + 1;
            ind_e_S = curr_edge + 2;
            ind_e_E = curr_edge + 3;
            ind_e_SE_up = curr_edge + 4;
            ind_e_NW_up += 2;
            ind_e_N = ind_e_NW_up + 1;
            ind_e_W_up = ind_e_W + (3*(n_x-1)+2)*(n_y-1) +  + 2*(n_y-1) + 1;
            ind_e_S_up = ind_e_S + (3*(n_x-1)+2)*(n_y-1) +  + 2*(n_y-1) + 1;

            face_W = curr_face;
            face_B = curr_face + 1;
            face_S = curr_face + 2;
            face_E = curr_face + 3;
            face_N += 1;
            face_U = face_B + (3*(n_x-1) + 1)*(n_y-1) + (n_y-1);
            
            xW = x_v[ix];
            ptSWB = (Point3D){xW,yS,zB};

            set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_W    );
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW   , ind_e_W    );
            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_SW_up);
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, ind_e_SW_up);
            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_S    );
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE   , ind_e_S    );

            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SW_up, face_W);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_W    , face_W);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_NW_up, face_W);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W_up , face_W);

            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S, face_B);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_E, face_B);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_N, face_B);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W, face_B);

            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_S    , face_S);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SE_up, face_S);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S_up , face_S);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_SW_up, face_S);

            infogrb = GrB_Matrix_setElement(*volumes,  1, face_E, curr_vol);
            infogrb = GrB_Matrix_setElement(*volumes, -1, face_B, curr_vol);
            infogrb = GrB_Matrix_setElement(*volumes, -1, face_S, curr_vol);
            infogrb = GrB_Matrix_setElement(*volumes, -1, face_W, curr_vol);
            infogrb = GrB_Matrix_setElement(*volumes,  1, face_N, curr_vol);
            infogrb = GrB_Matrix_setElement(*volumes,  1, face_U, curr_vol);

            curr_pt += 1;
            curr_edge += 3;
            curr_face += 3;
            curr_vol += 1;
        }
        //Treat east-most point and edge
        xW = x_v[n_x-1];
        ptSWB = (Point3D){xW,yS,zB};
        ind_pt_SW = curr_pt;
        ind_pt_NW = curr_pt + n_x;
        ind_pt_SW_up = curr_pt + n_x*n_y;

        ind_e_W = curr_edge;
        ind_e_SW_up = curr_edge + 1;
        ind_e_NW_up += 2;
        ind_e_W_up = ind_e_W + (3*(n_x-1)+2)*(n_y-1) +  + 2*(n_y-1) + 1;

        set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_W    );
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW   , ind_e_W    );
        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_SW_up);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, ind_e_SW_up);

        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_W    , curr_face);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SW_up, curr_face);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_NW_up, curr_face);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W_up , curr_face);

        curr_pt += 1;
        curr_edge += 2;
        curr_face += 1;
    

        //Treat north-most points and edges
        yS = y_v[n_y-1];
        for (ix = 0; ix < n_x-1; ix++){
            ind_pt_SW = curr_pt;
            ind_pt_SE = curr_pt + 1;
            ind_pt_SW_up = curr_pt + n_x*n_y;

            ind_e_SW_up = curr_edge;
            ind_e_S = curr_edge + 1;
            ind_e_SE_up = curr_edge + 2;
            ind_e_S_up = ind_e_S + (3*(n_x-1)+2)*(n_y-1) +  + 2*(n_y-1) + 1;

            xW = x_v[ix];
            ptSWB = (Point3D){xW,yS,zB};

            set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_SW_up);
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, ind_e_SW_up);
            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_S    );
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE   , ind_e_S    );

            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_S    , curr_face);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SE_up, curr_face);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S_up , curr_face);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_SW_up, curr_face);

            curr_pt += 1;
            curr_edge += 2;
            curr_face += 1;
        }
        xW = x_v[n_x-1];
        ptSWB = (Point3D){xW,yS,zB};
        ind_pt_SW = curr_pt;
        ind_pt_SW_up = curr_pt + n_x*n_y;
        
        set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);
        
        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , curr_edge);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, curr_edge);        

        curr_pt += 1;
        curr_edge += 1;
    }

    //Make the top-most cubes of the grid 
    zB = z_v[n_z-2];
    ind_e_W_up = curr_edge + (3*(n_x-1)+2)*(n_y-1) +  + 2*(n_y-1) + 1;
    face_U = curr_face + (3*(n_x-1) + 1)*(n_y-1) + (n_y-1);
    for (iy = 0; iy < n_y-2; iy++){
        yS = y_v[iy];
        for (ix = 0; ix < n_x-1; ix++){
            xW = x_v[ix];

            ind_pt_SW = curr_pt;
            ind_pt_SE = curr_pt + 1;
            ind_pt_NW = curr_pt + n_x;
            ind_pt_SW_up = curr_pt + n_x*n_y;

            ind_e_W = curr_edge;
            ind_e_SW_up = curr_edge + 1;
            ind_e_S = curr_edge + 2;
            ind_e_E = curr_edge + 3;
            ind_e_SE_up = curr_edge + 4;
            ind_e_NW_up = curr_edge + 3*n_x;
            ind_e_N = ind_e_NW_up + 1;
            ind_e_S_up = ind_e_W_up + 1;

            face_W = curr_face;
            face_B = curr_face + 1;
            face_S = curr_face + 2;
            face_E = curr_face + 3;
            face_N = face_S + (3*(n_x-1) + 1);

            ptSWB = (Point3D){xW,yS,zB};

            set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_W    );
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW   , ind_e_W    );
            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_SW_up);
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, ind_e_SW_up);
            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_S    );
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE   , ind_e_S    );

            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SW_up, face_W);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_W    , face_W);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_NW_up, face_W);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W_up , face_W);

            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S, face_B);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_E, face_B);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_N, face_B);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W, face_B);

            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_S    , face_S);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SE_up, face_S);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S_up , face_S);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_SW_up, face_S);

            infogrb = GrB_Matrix_setElement(*volumes,  1, face_E, curr_vol);
            infogrb = GrB_Matrix_setElement(*volumes, -1, face_B, curr_vol);
            infogrb = GrB_Matrix_setElement(*volumes, -1, face_S, curr_vol);
            infogrb = GrB_Matrix_setElement(*volumes, -1, face_W, curr_vol);
            infogrb = GrB_Matrix_setElement(*volumes,  1, face_N, curr_vol);
            infogrb = GrB_Matrix_setElement(*volumes,  1, face_U, curr_vol);

            curr_pt += 1;
            curr_edge += 3;
            curr_face += 3;
            curr_vol += 1;
            ind_e_W_up += 2;
            face_U += 1;
        }
        //Treat east-most point and edge
        xW = x_v[n_x-1];
        ptSWB = (Point3D){xW,yS,zB};
        ind_pt_SW = curr_pt;
        ind_pt_NW = curr_pt + n_x;
        ind_pt_SW_up = curr_pt + n_x*n_y;

        ind_e_W = curr_edge;
        ind_e_SW_up = curr_edge + 1;
        ind_e_NW_up = curr_edge +  + 3*n_x;

        set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_W    );
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW   , ind_e_W    );
        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_SW_up);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, ind_e_SW_up);

        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_W    , curr_face);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SW_up, curr_face);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_NW_up, curr_face);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W_up , curr_face);

        curr_pt += 1;
        curr_edge += 2;
        curr_face += 1;
        ind_e_W_up += 1;
    }

    //Treat north-most cells
    iy = n_y-2;
    yS = y_v[iy];
    ind_e_NW_up = curr_edge + (3*(n_x-1));
    face_N = curr_face + 3*(n_x-1);
    ind_e_S_up = curr_edge + (3*(n_x-1)+2) + (2*(n_x-1)+1)*(n_y-1) - 1;
    for (ix = 0; ix < n_x-1; ix++){
        ind_pt_SW = curr_pt;
        ind_pt_SE = curr_pt + 1;
        ind_pt_NW = curr_pt + n_x;
        ind_pt_SW_up = curr_pt + n_x*n_y;

        ind_e_W = curr_edge;
        ind_e_SW_up = curr_edge + 1;
        ind_e_S = curr_edge + 2;
        ind_e_E = curr_edge + 3;
        ind_e_SE_up = curr_edge + 4;
        ind_e_NW_up += 2;
        ind_e_N = ind_e_NW_up + 1;
        ind_e_S_up += 2;

        face_W = curr_face;
        face_B = curr_face + 1;
        face_S = curr_face + 2;
        face_E = curr_face + 3;
        face_N += 1;

        xW = x_v[ix];
        ptSWB = (Point3D){xW,yS,zB};

        set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_W    );
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW   , ind_e_W    );
        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_SW_up);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, ind_e_SW_up);
        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_S    );
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE   , ind_e_S    );

        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SW_up, face_W);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_W    , face_W);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_NW_up, face_W);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W_up , face_W);

        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S, face_B);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_E, face_B);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_N, face_B);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W, face_B);

        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_S    , face_S);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SE_up, face_S);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S_up , face_S);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_SW_up, face_S);

        infogrb = GrB_Matrix_setElement(*volumes,  1, face_E, curr_vol);
        infogrb = GrB_Matrix_setElement(*volumes, -1, face_B, curr_vol);
        infogrb = GrB_Matrix_setElement(*volumes, -1, face_S, curr_vol);
        infogrb = GrB_Matrix_setElement(*volumes, -1, face_W, curr_vol);
        infogrb = GrB_Matrix_setElement(*volumes,  1, face_N, curr_vol);
        infogrb = GrB_Matrix_setElement(*volumes,  1, face_U, curr_vol);

        curr_pt += 1;
        curr_edge += 3;
        curr_face += 3;
        curr_vol += 1;
        ind_e_W_up += 2;
        face_U += 1;
    }
    //Treat east-most point and edge
    xW = x_v[n_x-1];
    ptSWB = (Point3D){xW,yS,zB};
    ind_pt_SW = curr_pt;
    ind_pt_NW = curr_pt + n_x;
    ind_pt_SW_up = curr_pt + n_x*n_y;

    ind_e_W = curr_edge;
    ind_e_SW_up = curr_edge + 1;
    ind_e_NW_up += 2;

    set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

    infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_W    );
    infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW   , ind_e_W    );
    infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_SW_up);
    infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, ind_e_SW_up);

    infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_W    , curr_face);
    infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SW_up, curr_face);
    infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_NW_up, curr_face);
    infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W_up , curr_face);

    curr_pt += 1;
    curr_edge += 2;
    curr_face += 1;
    ind_e_W_up += 1;
    

    //Treat north-most points and edges
    yS = y_v[n_y-1];
    ind_e_S_up = curr_edge + (2*(n_x-1)+1) + (2*(n_x-1)+1)*(n_y-1) - 1;
    for (ix = 0; ix < n_x-1; ix++){
        ind_pt_SW = curr_pt;
        ind_pt_SE = curr_pt + 1;
        ind_pt_SW_up = curr_pt + n_x*n_y;

        ind_e_SW_up = curr_edge;
        ind_e_S = curr_edge + 1;
        ind_e_SE_up = curr_edge + 2;
        ind_e_S_up += 1;

        xW = x_v[ix];
        ptSWB = (Point3D){xW,yS,zB};
        set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_SW_up);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, ind_e_SW_up);
        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , ind_e_S    );
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE   , ind_e_S    );

        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_S    , curr_face);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_SE_up, curr_face);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S_up , curr_face);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_SW_up, curr_face);

        curr_pt += 1;
        curr_edge += 2;
        curr_face += 1;
    }
    xW = x_v[n_x-1];
    ptSWB = (Point3D){xW,yS,zB};
    ind_pt_SW = curr_pt;
    ind_pt_SW_up = curr_pt + n_x*n_y;

    set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);
    
    infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW   , curr_edge);
    infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SW_up, curr_edge);        

    curr_pt += 1;
    curr_edge += 1;

    //Finally, deal with the top-most faces, edges and points
    zB = z_v[n_z-1];
    for (iy = 0; iy < n_y-2; iy++){
        yS = y_v[iy];
        for (ix = 0; ix < n_x-1; ix++){
            xW = x_v[ix];

            ind_pt_SW = curr_pt;
            ind_pt_SE = curr_pt + 1;
            ind_pt_NW = curr_pt + n_x;

            ind_e_W = curr_edge;
            ind_e_S = curr_edge + 1;
            ind_e_E = curr_edge + 2;
            ind_e_N = curr_edge + (2*n_x);

            face_B = curr_face;

            ptSWB = (Point3D){xW,yS,zB};

            set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_W);
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW, ind_e_W);
            infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_S);
            infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE, ind_e_S);

            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S, face_B);
            infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_E, face_B);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_N, face_B);
            infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W, face_B);

            curr_pt += 1;
            curr_edge += 2;
            curr_face += 1;
        }
        //Treat east-most point and edge
        xW = x_v[n_x-1];
        ptSWB = (Point3D){xW,yS,zB};
        ind_pt_SW = curr_pt;
        ind_pt_NW = curr_pt + n_x;

        ind_e_W = curr_edge;

        set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_W);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW, ind_e_W);        

        curr_pt += 1;
        curr_edge += 1;
    }

    //Treat north-most cells
    yS = y_v[n_y-2];
    ind_e_N = curr_edge + (2*(n_x-1));
    for (ix = 0; ix < n_x-1; ix++){
        ind_pt_SW = curr_pt;
        ind_pt_SE = curr_pt + 1;
        ind_pt_NW = curr_pt + n_x;

        ind_e_W = curr_edge;
        ind_e_S = curr_edge + 1;
        ind_e_E = curr_edge + 2;
        ind_e_N += 1;

        face_B = curr_face ;
            
        xW = x_v[ix];
        ptSWB = (Point3D){xW,yS,zB};

        set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_W);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW, ind_e_W);
        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_S);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE, ind_e_S);

        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_S, face_B);
        infogrb = GrB_Matrix_setElement(*faces,  1, ind_e_E, face_B);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_N, face_B);
        infogrb = GrB_Matrix_setElement(*faces, -1, ind_e_W, face_B);

        curr_pt += 1;
        curr_edge += 2;
        curr_face += 1;
    }
    //Treat east-most point and edge
    xW = x_v[n_x-1];
    ptSWB = (Point3D){xW,yS,zB};
    ind_pt_SW = curr_pt;
    ind_pt_NW = curr_pt + n_x;

    ind_e_W = curr_edge;
    set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

    infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_W);
    infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_NW, ind_e_W);        

    curr_pt += 1;
    curr_edge += 1;
    

    //Treat north-most points and edges
    yS = y_v[n_y-1];
    for (ix = 0; ix < n_x-1; ix++){
        ind_pt_SW = curr_pt;
        ind_pt_SE = curr_pt + 1;

        ind_e_S = curr_edge;

        xW = x_v[ix];
        ptSWB = (Point3D){xW,yS,zB};

        set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);

        infogrb = GrB_Matrix_setElement(*edges, -1, ind_pt_SW, ind_e_S);
        infogrb = GrB_Matrix_setElement(*edges,  1, ind_pt_SE, ind_e_S);        

        curr_pt += 1;
        curr_edge += 1;
    }
    xW = x_v[n_x-1];
    ptSWB = (Point3D){xW,yS,zB};
    ind_pt_SW = curr_pt;
    set_ith_elem_vec_pts3D(vertices, ind_pt_SW, &ptSWB);
            

    status_face = alloc_empty_vec_int();
    status_face->data = (long int*) calloc(nb_faces, sizeof(long int));
    status_face->capacity = nb_faces;
    status_face->size = nb_faces;
    return new_Polyhedron3D_vefvs(vertices, edges, faces, volumes, status_face);
}

void copy_Polyhedron3D(const Polyhedron3D *src, Polyhedron3D *dest){
    if (src == NULL){
        if (dest != NULL) {
            dealloc_vec_pts3D(dest->vertices);
            GrB_free(dest->edges);
            GrB_free(dest->faces);
            GrB_free(dest->volumes);
            dealloc_vec_int(dest->status_face);
            dest->vertices = NULL;
            dest->edges = NULL;
            dest->faces = NULL;
            dest->volumes = NULL;
            dest->status_face = NULL;
        }
    } else {
        if (dest != NULL){
            GrB_free(dest->edges);
            GrB_free(dest->faces);
            GrB_free(dest->volumes);
            copy_vec_pts3D(src->vertices, dest->vertices);
            GrB_Matrix_dup(dest->edges, *(src->edges));
            GrB_Matrix_dup(dest->faces, *(src->faces));
            GrB_Matrix_dup(dest->volumes, *(src->volumes));
            copy_vec_int(src->status_face, dest->status_face);
        }
    }
}

void dealloc_Polyhedron3D(Polyhedron3D* p){
    if(p){
        dealloc_vec_pts3D(p->vertices);
        GrB_free(p->edges);
        GrB_free(p->faces);
        GrB_free(p->volumes);
        dealloc_vec_int(p->status_face);
    }
}