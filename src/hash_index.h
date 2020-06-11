#ifndef _HASH_INDEX_H
#define _HASH_INDEX_H

#include "read_seeding.h"

typedef struct _HashEntry
{
    int value_index; //order of kmer appear
    int *value;//kmer position in reference
}hash_entry;

typedef struct _HashTable
{
    hash_entry *bucket;
}hash_table;

int bucket_num;
int value_num;
int hit_num;

void initHashTable(uint8_t h_k);
void freeHashTable();
void local_hash_process(uint8_t **target, uint32_t target_len, anchored_exons** anchored_exon, uint32_t* anchored_exon_num, uint32_t *uniseed2_length, uint8_t tid);

#endif
