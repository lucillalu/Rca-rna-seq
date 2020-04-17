/*************************************************************************
	> File Name: read_seeding.h
	> Author: 
	> Mail: 
 ************************************************************************/

#ifndef _READ_SEEDING_H
#define _READ_SEEDING_H

#include <stdio.h>
#include <stdint.h>
#include <getopt.h> 
#include <zlib.h>

#define INDEX_KMER 22
#define SEEDING_KMER 15
#define SEED_STEP 5
#define MAX_UNI_POS 50
#define MIN_FRAG_DIS 20
#define BATCH_SIZE 100000
#define TEMP_FILE_PERFIRX "1pass_"
#define MAX_READLEN 1000000 //2^15+1 > readlen_max = 30000
#define LOCAL_HASH_KMER 8
#define MAX_READ_JOIN_GAP 2000

typedef struct READ
{
	char* read_seq;
	char* name;
	uint32_t read_length;
} READ_t;

typedef struct HITSEED
{
	uint64_t uid;
	uint32_t seed_id;
	uint32_t read_pos;
	uint32_t uni_pos_off;
	uint32_t length;
	uint32_t pos_n;
}hit_seed;

typedef struct UNI_SEED
{
	uint32_t read_begin;
	uint32_t read_end;
	uint32_t seed_id; //record the first seed_id of those mems which can be merged
	uint32_t ref_begin;
	uint32_t ref_end;
}uni_seed;

typedef struct ANCHORED_EXON
{
	uint32_t id;
	uint32_t ref_begin;
	uint32_t ref_end;
	uint32_t cov;
	uint32_t length;
}anchored_exons;

typedef struct HCAE
{
	uint32_t id;
	uint32_t ref_begin;
	uint32_t ref_end;
	uint32_t cov;
	uint32_t length;
}hcaes;

typedef struct PATH
{
	float dist;
	int32_t pre_node;
	float mean;
	uint32_t read_end;
} PATH_t;

typedef struct{
    int batch_size;	
	uint8_t read_type;
	uint16_t Eindel;
	uint8_t seed_step;
	uint8_t k_t; //index kmer
	uint8_t seed_k_t; //alignment kmer
	uint8_t hash_kmer;
 	uint16_t pos_n_max;
	char *temp_file_perfix;
	int readlen_max;
	int max_read_join_gap;
}param_map;

hit_seed** hitseed;
// vertex_u** vertexu;
anchored_exons** anchored_exon;
READ_t* query_info;

//variable in this file
uint8_t k_r;
uint8_t re_b;
uint8_t re_bt;
uint8_t re_2bt;
uint8_t top_n;
uint8_t seed_step;
uint8_t hash_kmer;
int8_t seed_offset;
uint16_t pos_n_max;
uint16_t uni_pos_n_max;
int batch_size;
int seed_num;
int max_read_join_gap;
char temp_hit_dir[1024];
char temp_ae_dir[1024];
char temp_uni_dir[1024];
char temp_uni2_dir[1024];
char temp_hcae_dir[1024];
char temp_re_dir[1024];
FILE *fp_tff;
FILE *fp_ae;
FILE *fp_tfu;
FILE *fp_tfu2;
FILE *fp_hcae;
FILE *fp_re;

//global variable
extern uni_seed** uniseed;
extern uni_seed** uniseed3;
extern hcaes** hcae;
extern uint32_t anchored_exon_num[2];
extern uint8_t k_t;
extern uint8_t seed_k_t;
extern uint32_t new_seed_cnt;
extern uint16_t Eindel;
extern int waitingLen;
extern int readlen_max;

uint64_t read_bit1[2][((MAX_READLEN - 1) >> 5) + 1];

int help_usage();
int desalt_aln(int argc, char *argv[], const char *version);

#endif
