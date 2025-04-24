#include "vector_double.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

Vector_double* alloc_empty_vec_double(){
    Vector_double* v = (Vector_double*)malloc(sizeof(Vector_double));
    v->data = NULL;
    v->size = 0;
    v->capacity = 0;

    return v;
}

Vector_double* alloc_with_capacity_vec_double(uint64_t size){
    Vector_double* v = (Vector_double*)malloc(sizeof(Vector_double));
    v->data = (double*)malloc(size*sizeof(double));
    v->size = 0;
    v->capacity = size;

    return v;
}

Vector_double* alloc_with_init_vec_double(const double* data, uint64_t size){
    Vector_double* v = (Vector_double*)malloc(sizeof(Vector_double));
    v->data = (double*)malloc(size*sizeof(double));
    memcpy(v->data, data, size*sizeof(double));
    v->size = size;
    v->capacity = size;

    return v;
}

void double_capacity_vec_double(Vector_double* v){
    v->capacity *= 2;
    v->data = (double*) realloc(v->data, v->capacity*sizeof(double));
}

void push_back_vec_double(Vector_double* v, const double* point){
    if ((v->size == 0) || (v->data == NULL)){
        free(v->data);
        *v = *alloc_with_init_vec_double(point, 1);
    }
    else if (v->size >= v->capacity)
    {
        double_capacity_vec_double(v);
        v->data[v->size] = *point;
        v->size += 1;
    }
    else {
        v->data[v->size] = *point;
        v->size += 1;
    }
}

void dealloc_vec_double(Vector_double* v){
    if(v!=NULL){
        free(v->data);
        v->data = NULL;
        v->capacity = 0;
        v->size = 0;
    }
}

double* get_ith_elem_vec_double(const Vector_double* v, uint64_t i){
    double* p;
    if (i<v->size) {
        return v->data + i;
    } else {
        p = (double*) malloc(sizeof(double));
        *p = 0.0/0.0;
        return p;
    }
}

void set_ith_elem_vec_double(Vector_double* v, uint64_t i, double* d){
    while(i >= v->capacity) double_capacity_vec_double(v);
    if (i >= v->size) v->size = i+1;
    
    v->data[i] = *d;
}

void copy_vec_double(const Vector_double* src, Vector_double* dest){
    if (src == NULL){
        if (dest != NULL) {
            free(dest->data);
            dest->size = 0;
            dest->capacity = 0;
        }
    } else {
        if (dest != NULL){
            free(dest->data);
            dest->data = (double*)malloc(src->capacity*sizeof(double));
            memcpy(dest->data, src->data, src->size*sizeof(double));
            dest->size = src->size;
            dest->capacity = src->capacity;
        }
    }
}

void print_vec_double(const Vector_double* v){
    unsigned int i;
    double* vi;

    printf("v = [");
    for (i = 0; i<v->size-1; i++){
        vi = get_ith_elem_vec_double(v, i);
        printf("%.3e, ", *vi);
    }
    vi = get_ith_elem_vec_double(v, v->size-1);
    printf("%.3e]\n", *vi);
}