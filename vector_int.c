#include "vector_int.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

Vector_int* alloc_empty_vec_int(){
    Vector_int* v = (Vector_int*)malloc(sizeof(Vector_int));
    v->data = NULL;
    v->size = 0;
    v->capacity = 0;

    return v;
}

Vector_int* alloc_with_capacity_vec_int(uint64_t size){
    Vector_int* v = (Vector_int*)malloc(sizeof(Vector_int));
    v->data = (long int*)malloc(size*sizeof(long int));
    v->size = 0;
    v->capacity = size;

    return v;
}

Vector_int* alloc_with_init_vec_int(const long int* data, uint64_t size){
    Vector_int* v = (Vector_int*)malloc(sizeof(Vector_int));
    v->data = (long int*)malloc(size*sizeof(long int));
    memcpy(v->data, data, size*sizeof(long int));
    v->size = size;
    v->capacity = size;

    return v;
}

void double_capacity_vec_int(Vector_int* v){
    v->capacity *= 2;
    v->data = (long int*) realloc(v->data, v->capacity*sizeof(long int));
}

void push_back_vec_int(Vector_int* v, const long int* point){
    if ((v->size == 0) || (v->data == NULL)){
        free(v->data);
        *v = *alloc_with_init_vec_int(point, 1);
    }
    else if (v->size >= v->capacity)
    {
        double_capacity_vec_int(v);
        v->data[v->size] = *point;
        v->size += 1;
    }
    else {
        v->data[v->size] = *point;
        v->size += 1;
    }
}

void dealloc_vec_int(Vector_int* v){
    if(v!=NULL){
        free(v->data);
        v->data = NULL;
        v->capacity = 0;
        v->size = 0;
    }
}

long int* get_ith_elem_vec_int(const Vector_int* v, uint64_t i){
    long int* p;
    if (i<v->size) {
        return v->data + i;
    } else {
        p = (long int*) malloc(sizeof(long int));
        *p = 0.0/0.0;
        return p;
    }
}

void set_ith_elem_vec_int(Vector_int* v, uint64_t i, long int* d){
    while(i >= v->capacity) double_capacity_vec_int(v);
    if (i >= v->size) v->size = i+1;
    
    v->data[i] = *d;
}

void copy_vec_int(const Vector_int* src, Vector_int* dest){
    if (src == NULL){
        if (dest != NULL) {
            free(dest->data);
            dest->size = 0;
            dest->capacity = 0;
        }
    } else {
        if (dest != NULL){
            free(dest->data);
            dest->data = (long int*)malloc(src->capacity*sizeof(long int));
            memcpy(dest->data, src->data, src->size*sizeof(long int));
            dest->size = src->size;
            dest->capacity = src->capacity;
        }
    }
}

void print_vec_int(const Vector_int* v){
    uint64_t i;
    long int* vi;

    printf("v = [");
    for (i = 0; i<v->size-1; i++){
        vi = get_ith_elem_vec_int(v, i);
        printf("%ld, ", *vi);
    }
    vi = get_ith_elem_vec_int(v, v->size-1);
    printf("%ld]\n", *vi);
}

Vector_uint* alloc_empty_vec_uint(){
    Vector_uint* v = (Vector_uint*)malloc(sizeof(Vector_uint));
    v->data = NULL;
    v->size = 0;
    v->capacity = 0;

    return v;
}

Vector_uint* alloc_with_capacity_vec_uint(uint64_t size){
    Vector_uint* v = (Vector_uint*)malloc(sizeof(Vector_uint));
    v->data = (long unsigned int*)malloc(size*sizeof(long unsigned int));
    v->size = 0;
    v->capacity = size;

    return v;
}

Vector_uint* alloc_with_init_vec_uint(const long unsigned int* data, uint64_t size){
    Vector_uint* v = (Vector_uint*)malloc(sizeof(Vector_uint));
    v->data = (long unsigned int*)malloc(size*sizeof(long unsigned int));
    memcpy(v->data, data, size*sizeof(long unsigned int));
    v->size = size;
    v->capacity = size;

    return v;
}

void double_capacity_vec_uint(Vector_uint* v){
    v->capacity *= 2;
    v->data = (long unsigned int*) realloc(v->data, v->capacity*sizeof(long unsigned int));
}

void push_back_vec_uint(Vector_uint* v, const long unsigned int* point){
    if ((v->size == 0) || (v->data == NULL)){
        free(v->data);
        *v = *alloc_with_init_vec_uint(point, 1);
    }
    else if (v->size >= v->capacity)
    {
        double_capacity_vec_uint(v);
        v->data[v->size] = *point;
        v->size += 1;
    }
    else {
        v->data[v->size] = *point;
        v->size += 1;
    }
}

void dealloc_vec_uint(Vector_uint* v){
    if(v!=NULL){
        free(v->data);
        v->data = NULL;
        v->capacity = 0;
        v->size = 0;
    }
}

long unsigned int* get_ith_elem_vec_uint(const Vector_uint* v, uint64_t i){
    long unsigned int* p;
    if (i<v->size) {
        return v->data + i;
    } else {
        p = (long unsigned int*) malloc(sizeof(long unsigned int));
        *p = 0.0/0.0;
        return p;
    }
}

void set_ith_elem_vec_uint(Vector_uint* v, uint64_t i, long unsigned int* d){
    while(i >= v->capacity) double_capacity_vec_uint(v);
    if (i >= v->size) v->size = i+1;
    
    v->data[i] = *d;
}

void copy_vec_uint(const Vector_uint* src, Vector_uint* dest){
    if (src == NULL){
        if (dest != NULL) {
            free(dest->data);
            dest->size = 0;
            dest->capacity = 0;
        }
    } else {
        if (dest != NULL){
            free(dest->data);
            dest->data = (long unsigned int*)malloc(src->capacity*sizeof(long unsigned int));
            memcpy(dest->data, src->data, src->size*sizeof(long unsigned int));
            dest->size = src->size;
            dest->capacity = src->capacity;
        }
    }
}

void print_vec_uint(const Vector_uint* v){
    uint64_t i;
    long unsigned int* vi;

    printf("v = [");
    for (i = 0; i<v->size-1; i++){
        vi = get_ith_elem_vec_uint(v, i);
        printf("%ld, ", *vi);
    }
    vi = get_ith_elem_vec_uint(v, v->size-1);
    printf("%ld]\n", *vi);
}