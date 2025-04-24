#ifndef VECTOR_INT_H
#define VECTOR_INT_H

#include <stdint.h>

typedef struct {
    long int* data;
    uint64_t size, capacity;
} Vector_int;

void push_back_vec_int(Vector_int* v, const long int* point);
void double_capacity_vec_int(Vector_int* v);
Vector_int* alloc_empty_vec_int();
Vector_int* alloc_with_init_vec_int(const long int* points, uint64_t size);
Vector_int* alloc_with_capacity_vec_int(uint64_t size);
void dealloc_vec_int(Vector_int* v);
long int* get_ith_elem_vec_int(const Vector_int* v, uint64_t i);
void set_ith_elem_vec_int(Vector_int* v, uint64_t i, long int* d);
void copy_vec_int(const Vector_int* src, Vector_int* dest);
void print_vec_int(const Vector_int* v);


typedef struct {
    long unsigned int* data;
    uint64_t size, capacity;
} Vector_uint;

void push_back_vec_uint(Vector_uint* v, const long unsigned int* point);
void double_capacity_vec_uint(Vector_uint* v);
Vector_uint* alloc_empty_vec_uint();
Vector_uint* alloc_with_init_vec_uint(const long unsigned int* points, uint64_t size);
Vector_uint* alloc_with_capacity_vec_uint(uint64_t size);
void dealloc_vec_uint(Vector_uint* v);
long unsigned int* get_ith_elem_vec_uint(const Vector_uint* v, uint64_t i);
void set_ith_elem_vec_uint(Vector_uint* v, uint64_t i, long unsigned int* d);
void copy_vec_uint(const Vector_uint* src, Vector_uint* dest);
void print_vec_uint(const Vector_uint* v);


#endif