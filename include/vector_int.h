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
void set_ith_elem_vec_int(Vector_int* v, uint64_t i, const long int* d);
void copy_vec_int(const Vector_int* src, Vector_int* dest);
void print_vec_int(const Vector_int* v);
Vector_int* cat_vec_int(const Vector_int* v1, const Vector_int* v2);


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
Vector_uint* cat_vec_uint(const Vector_uint* v1, const Vector_uint* v2);

typedef struct {
    int64_t* data;
    uint64_t size, capacity;
} Vector_int64;

void push_back_vec_int64(Vector_int64* v, const int64_t* point);
void push_back_unique_vec_int64(Vector_int64* v, const int64_t* point);
void double_capacity_vec_int64(Vector_int64* v);
Vector_int64* alloc_empty_vec_int64();
Vector_int64* alloc_with_init_vec_int64(const int64_t* points, uint64_t size);
Vector_int64* alloc_with_capacity_vec_int64(uint64_t size);
void dealloc_vec_int64(Vector_int64* v);
int64_t* get_ith_elem_vec_int64(const Vector_int64* v, uint64_t i);
void set_ith_elem_vec_int64(Vector_int64* v, uint64_t i, int64_t* d);
void copy_vec_int64(const Vector_int64* src, Vector_int64* dest);
void print_vec_int64(const Vector_int64* v);
void sort_vec_int64(Vector_int64* v);
Vector_int64* cat_vec_int64(const Vector_int64* v1, const Vector_int64* v2);


#endif