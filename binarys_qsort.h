#ifndef _BINARYS_QSORT_H
#define _BINARYS_QSORT_H

int binsearch_range(uint64_t key, uint32_t *v, int64_t n, int64_t *range, int8_t k_off);

int64_t binsearch_interval_unipath64(uint64_t x, uint64_t v[], uint64_t n);

int64_t binsearch_interval_unipath(uint32_t x, uint32_t v[], uint32_t n);

int compare_uniid(const void * a, const void * b);

int compare_uniseed(const void *a , const void *b);

int compare_uniseed2(const void *a , const void *b);

#endif
