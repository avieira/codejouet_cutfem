#include "array_double.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

Array_double* alloc_empty_arr_double(){
    Array_double* v = (Array_double*)malloc(sizeof(Array_double));
    v->data = NULL;
    v->nrows = 0;
    v->ncols = 0;
    v->capacity = 0;

    return v;
}

//Allocated with all data = 0
Array_double* alloc_with_capacity_arr_double(uint64_t nrows, uint64_t ncols){
    Array_double* v = (Array_double*)malloc(sizeof(Array_double));
    v->data = (double*)calloc(nrows*ncols, sizeof(double));
    v->nrows = nrows;
    v->ncols = ncols;
    v->capacity = nrows*ncols;

    return v;
}

Array_double* alloc_with_init_arr_double(const double* data, uint64_t ncols, uint64_t nrows){
    Array_double* v = (Array_double*)malloc(sizeof(Array_double));
    v->data = (double*)calloc(ncols*nrows, sizeof(double));
    memcpy(v->data, data, ncols*nrows*sizeof(double));
    v->nrows = nrows;
    v->ncols = ncols;
    v->capacity = ncols*nrows;

    return v;
}

void double_capacity_arr_double(Array_double* v){
    v->capacity *= 2;
    v->data = (double*) realloc(v->data, v->capacity*sizeof(double));
}

void dealloc_arr_double(Array_double* v){
    if(v!=NULL){
        free(v->data);
        v->data = NULL;
        v->capacity = 0;
        v->nrows = 0;
        v->ncols = 0;
    }
}

double* get_ijth_elem_arr_double(const Array_double* v, const uint64_t i, const uint64_t j){
    double* p;
    if ((i<v->nrows) && (j<v->ncols)) {
        return v->data + i*v->ncols + j;
    } else {
        p = (double*) malloc(sizeof(double));
        *p = 0.0/0.0;
        return p;
    }
}

void set_ijth_elem_arr_double(Array_double* v, const uint64_t i, const uint64_t j, double* d){
    if (j >= v->ncols){
        printf("ERROR in set set_ijth_elem_arr_double : j bigger than number of columns!");
    }
    while(i*v->ncols+j >= v->capacity) double_capacity_arr_double(v);
    if (i >= v->nrows) v->nrows = i+1;
    
    v->data[i*v->ncols+j] = *d;
}

void copy_arr_double(const Array_double* src, Array_double* dest){
    if (src == NULL){
        if (dest != NULL) {
            free(dest->data);
            dest->nrows = 0;
            dest->ncols = 0;
            dest->capacity = 0;
        }
    } else {
        if (dest != NULL){
            free(dest->data);
            dest->data = (double*)malloc(src->capacity*sizeof(double));
            memcpy(dest->data, src->data, src->nrows*src->ncols*sizeof(double));
            dest->nrows = src->nrows;
            dest->ncols = src->ncols;
            dest->capacity = src->capacity;
        }
    }
}

void print_arr_double(const Array_double* v){
    unsigned long int i, j;
    double* vi;

    printf("v = [");
    for (i = 0; i<v->nrows-1; i++){
        for (j = 0; j<v->ncols-1; j++){
            vi = get_ijth_elem_arr_double(v, i, j);
            printf("%.3e, ", *vi);
        }
        vi = get_ijth_elem_arr_double(v, i, v->ncols-1);
        printf("%.3e\n", *vi);
    }
    for (j = 0; j<v->ncols-1; j++){
        vi = get_ijth_elem_arr_double(v, v->nrows-1, j);
        printf("%.3e, ", *vi);
    }
    vi = get_ijth_elem_arr_double(v, v->nrows-1, v->ncols-1);
    printf("%.3e]\n", *vi);
}