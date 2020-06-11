#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>

#include "read_seeding.h"
#include "bit_operation.h"
#include "binarys_qsort.h"
#include "load_unipath_size.h"
#include "bseq.h"
#include "rrs_index.h"
#include "ktime.h"
#include "hash_index.h"
#include "graph.h"
#include "poa.h"

// #define _DEBUG
// #define _OUTRE

//variable extern
uint8_t thread_n;
uni_seed*** uniseed = NULL;
uni_seed*** uniseed3 = NULL;
anchored_exons*** anchored_exon = NULL;
int readlen_max;
uint8_t k_t;
uint8_t seed_k_t;
uint32_t new_seed_cnt;

uint16_t Eindel ; //less than 0.01% of intron length less than 20bp
int waitingLen;

//varibale in this file
uint8_t k_first_level = 14;
int TOTAL_READ_COUNTs;
char command[1024];
uint8_t read_type = 5;
int *read_len = NULL;
int POS_N_MAX = 0;
int THREAD_READ_I;

static void init_map_param(param_map *opt)
{

    opt->batch_size = 655350;	
	opt->Eindel = 30;
	opt->k_t = 22;
	opt->seed_k_t = 15;
	opt->hash_kmer = 8;
	opt->seed_step = 5;
	opt->pos_n_max = 50;
	opt->readlen_max = MAX_READLEN;
	opt->max_read_join_gap = 3000;
	opt->thread_n = 1;
}

static inline void find_anchored_exon(uni_seed **uniseed, uint32_t *uniseed_length ,uint32_t *anchored_exon_num, uint8_t tid)
{
	uint8_t r_c;
	uint32_t su_i, j, k, s1, e1, co_num;
	uint32_t id, cov;
	uint32_t max_intron_length = 200000;
	uint32_t max_len = 0, tmp_id = 0;
	#ifdef _DEBUG
	fprintf(stderr, "id\tref_begin\tref_end\tcov\tlength\n");
	#endif
	for (r_c = 0; r_c < 2; ++r_c)
	{	
		#ifdef _DEBUG
		fprintf(stdout, "strand = %d\n",r_c);
		fprintf(fp_ae, "strand = %d\n",r_c);
		#endif
        if(uniseed_length[r_c] == 0)
        {
            anchored_exon_num[r_c] = 0;
            continue;
        }

        qsort(uniseed[r_c], uniseed_length[r_c], sizeof(uni_seed), compare_uniseed);
        //sort by ref_begin position

        j = 0; 
        id = 0;
        // printf("uniseed_length[%d] = %d\n", r_c, uniseed_length[r_c]);
        while(j < uniseed_length[r_c])
        {
        	s1 = j;
        	cov = uniseed[r_c][j].ref_end - uniseed[r_c][j].ref_begin + 1;
        	j++;
        	while(j < uniseed_length[r_c])
        	{
        		int diff = (int)(uniseed[r_c][j].ref_begin - uniseed[r_c][j-1].ref_end - 1);//distance
        		if (diff > Eindel)//too far
					break;
				else
				{
					//*****************************cov
					cov += uniseed[r_c][j].ref_end - uniseed[r_c][j].ref_begin + 1;
					j++;

				}
			}
			e1 = j - 1;
			if (e1 - s1 + 1 > 2)//condition 1
			{
				anchored_exon[tid][r_c][id].id = id;
				anchored_exon[tid][r_c][id].ref_begin = uniseed[r_c][s1].ref_begin;
				anchored_exon[tid][r_c][id].ref_end = uniseed[r_c][e1].ref_end;
				anchored_exon[tid][r_c][id].cov = cov;
				anchored_exon[tid][r_c][id].length = uniseed[r_c][e1].ref_end - uniseed[r_c][s1].ref_begin + 1;
				// printf("%d\t%d\t%d\t%d\t%d\n", anchored_exon[r_c][id].id,anchored_exon[r_c][id].ref_begin,anchored_exon[r_c][id].ref_end,anchored_exon[r_c][id].cov,anchored_exon[r_c][id].length);
				if (max_len < anchored_exon[tid][r_c][id].length)
				{
					max_len = anchored_exon[tid][r_c][id].length;
					tmp_id = id;
				}
				id++;
			}
			else if (e1 - s1 + 1 <= 2)
			{
				co_num = 0;
				for (k = s1; k <= e1; ++k)
				{
					co_num += uniseed[r_c][k].num;
				}
				if (co_num > 2)//condition 2
				{
					anchored_exon[tid][r_c][id].id = id;
					anchored_exon[tid][r_c][id].ref_begin = uniseed[r_c][s1].ref_begin;
					anchored_exon[tid][r_c][id].ref_end = uniseed[r_c][e1].ref_end;
					anchored_exon[tid][r_c][id].cov = cov;
					anchored_exon[tid][r_c][id].length = uniseed[r_c][e1].ref_end - uniseed[r_c][s1].ref_begin + 1;
					// printf("%d\t%d\t%d\t%d\t%d\n", anchored_exon[r_c][id].id,anchored_exon[r_c][id].ref_begin,anchored_exon[r_c][id].ref_end,anchored_exon[r_c][id].cov,anchored_exon[r_c][id].length);
					if (max_len < anchored_exon[tid][r_c][id].length)
					{
						max_len = anchored_exon[tid][r_c][id].length;
						tmp_id = id;
					}
					id++;
				}
			}
        }
        anchored_exon_num[r_c] = id;
        // printf("max_len = %d\n",max_len);
        id = 0;
        for(j = 0; j < anchored_exon_num[r_c]; ++j)
        {
        	if (j != tmp_id && abs(anchored_exon[tid][r_c][j].ref_begin - anchored_exon[tid][r_c][tmp_id].ref_begin) > max_intron_length)
        	{
        		continue;
        	}
        	else
        	{
        		anchored_exon[tid][r_c][id].id = id;
				anchored_exon[tid][r_c][id].ref_begin = anchored_exon[tid][r_c][j].ref_begin;
				anchored_exon[tid][r_c][id].ref_end = anchored_exon[tid][r_c][j].ref_end;
				anchored_exon[tid][r_c][id].cov = anchored_exon[tid][r_c][j].cov;
				anchored_exon[tid][r_c][id].length = anchored_exon[tid][r_c][j].length;
        		#ifdef _DEBUG
         		fprintf(fp_ae, "%d\t%d\t%d\t%d\t%d\n", anchored_exon[tid][r_c][id].id,anchored_exon[tid][r_c][id].ref_begin,anchored_exon[tid][r_c][id].ref_end,anchored_exon[tid][r_c][id].cov,anchored_exon[tid][r_c][id].length);
        		fprintf(stdout, "%d\t%d\t%d\t%d\t%d\n", anchored_exon[tid][r_c][id].id,anchored_exon[tid][r_c][id].ref_begin,anchored_exon[tid][r_c][id].ref_end,anchored_exon[tid][r_c][id].cov,anchored_exon[r_c][id].length);
        		#endif
        		id++;
        	}
        }
        anchored_exon_num[r_c] = id;
	}
}

static inline void map_to_ref(co_hitseed **hitseed2, uint32_t *memid, uint32_t *uniseed_length,uint8_t tid)
{
	uint8_t r_i;
	uint32_t i, m;
	uint32_t su_i = 0;
	for (r_i = 0; r_i < 2; ++r_i)//r_i=0, 1
	{
        if(memid[r_i] == 0)
        {
            uniseed_length[r_i] = 0;
            continue;
        }
		su_i = 0;
		// printf("memid[%d]=%d\n",r_i,memid[r_i]);
		for (i = 0; i < memid[r_i]; ++i)
		{
			// printf("memid=%d\n",i);
			for (m = 0; m < hitseed2[r_i][i].pos_n; ++m)//record the reference position
			{
				uniseed[tid][r_i][su_i].seed_id = i;
				uniseed[tid][r_i][su_i].read_begin = hitseed2[r_i][i].read_pos;
				uniseed[tid][r_i][su_i].read_end = hitseed2[r_i][i].read_pos + hitseed2[r_i][i].length1 - 1;
				uniseed[tid][r_i][su_i].ref_begin =  buffer_p[m + buffer_pp[hitseed2[r_i][i].uid]] + hitseed2[r_i][i].uni_pos_off - 1;
				uniseed[tid][r_i][su_i].ref_end = uniseed[tid][r_i][su_i].ref_begin + hitseed2[r_i][i].length2 - 1;
                uniseed[tid][r_i][su_i].num = hitseed2[r_i][i].num;
                // printf("uniseed[%d][%d][%d].seed_id = %d\n",r_i,su_i,uniseed[r_i][su_i].seed_id);
                // printf("uniseed[%d][%d][%d].read_begin = %d\n",r_i,su_i,uniseed[r_i][su_i].read_begin);
                // printf("uniseed[%d][%d][%d].read_end = %d\n",r_i,su_i,uniseed[r_i][su_i].read_end);
                // printf("uniseed[%d][%d][%d].ref_begin = %d\n",r_i,su_i,uniseed[r_i][su_i].ref_begin);
                // printf("uniseed[%d][%d][%d].ref_end = %d\n",r_i,su_i,uniseed[r_i][su_i].ref_end);
				su_i++;
			}
		}
		uniseed_length[r_i] = su_i;

		qsort(uniseed[tid][r_i], uniseed_length[r_i], sizeof(uni_seed), compare_uniseed2);
        //sort by read_begin position

        #ifdef _DEBUG
		fprintf(fp_tfu, "strand=%d\n",r_i);
        for (su_i = 0; su_i < uniseed_length[r_i]; ++su_i)
        {
        	fprintf(fp_tfu, "%d\t%d\t%d\t%d\t%d\n", uniseed[tid][r_i][su_i].seed_id,uniseed[tid][r_i][su_i].read_begin,uniseed[tid][r_i][su_i].read_end,uniseed[tid][r_i][su_i].ref_begin,uniseed[tid][r_i][su_i].ref_end);
        }
        #endif
	}
}

int get_split_read(uint32_t dist_max_index, PATH_t *dist_path, uint32_t *read_end)
{
	int32_t j;
	uint32_t cnt = 0;
	float tmp_mean;
	uint32_t tmp_read_end;

	tmp_mean = dist_path[dist_max_index].mean;
	tmp_read_end = dist_path[dist_max_index].read_end;
	if (tmp_mean != 0 && tmp_read_end != 0)
	{
		read_end[cnt] = dist_path[dist_max_index].read_end;
		cnt++;
	}

	j = dist_path[dist_max_index].pre_node;
	while(j != -1)
	{
		if (dist_path[j].mean != tmp_mean || dist_path[j].read_end != tmp_read_end)
		{
			read_end[cnt] = dist_path[j].read_end;
			cnt++;
			tmp_mean = dist_path[j].mean;
			tmp_read_end = dist_path[j].read_end;
		}
		j = dist_path[j].pre_node;
	}

	#ifdef _OUTRE
	for(j = 0; j < cnt; ++j)
	{
		fprintf(fp_re, "%d\t", read_end[j]);
	}
	fprintf(fp_re, "\n");
	#endif

	return cnt;
}

int get_read_seq(uint8_t *read, uint8_t* read_seq, uint32_t start, uint32_t length)
{
	uint32_t i, j = 0;
	for (i = start; i < start + length; ++i)
	{
		read_seq[j] = read[i];
		j++;
	}
	return 0;
}

int find_hit_in_HCAE(uint32_t seqi, uint32_t read_length, uint8_t tid, uint32_t *anchored_exon_num)
{
	int i,j,k;
	uint8_t r_i;
	//hash_kmer = 8
	uint8_t *qseq0[2];
	uint32_t ref_begin[2],ref_end[2];
	uint32_t uniseed2_length[2] = {0,0}; 

	uint32_t read_length_msa = read_length;

	for(r_i = 0;r_i < 2; ++r_i)
	{
		qseq0[r_i] = (uint8_t* )calloc(read_length, sizeof(uint8_t));
	}

	for (i = 0; i < read_length; ++i)
	{
		qseq0[0][i] = nt_table[(uint8_t)(query_info[seqi].read_seq)[i]];
		if (qseq0[0][i] >= 4)
		{
			qseq0[0][i] = rand()%4;
		}
		qseq0[1][read_length - 1 - i] = 3 - qseq0[0][i];
	}
	// #ifdef _DEBUG
	// fprintf(stderr, "[Process-Info] Processing local hash...\n");
	// #endif
	local_hash_process(qseq0, read_length, anchored_exon[tid], anchored_exon_num, uniseed2_length, tid);
	// printf("uniseed2_length[0] = %d\n", uniseed2_length[0]);
	// printf("uniseed2_length[1] = %d\n", uniseed2_length[1]);

	float dist_0 = 0;
	float dist_1 = 0;
	uint32_t max_index0 = 0;
	uint32_t max_index1 = 0;
 	PATH_t* path0;
    PATH_t* path1;
    path0 = (PATH_t* )malloc(uniseed2_length[0]*sizeof(PATH_t));
    path1 = (PATH_t* )malloc(uniseed2_length[1]*sizeof(PATH_t));
	uint32_t **read_end;
	int **seq_len;
	uint8_t ***read_seq;
	int cnt0 = 0,cnt1 = 0;

    uint32_t num = (uniseed2_length[0] < uniseed2_length[1]) ? uniseed2_length[1]: uniseed2_length[0];
	uint8_t **out_degree = (uint8_t** )calloc(2, sizeof(uint8_t* ));
    for ( i = 0; i < 2; ++i)
	{
		out_degree[i] = (uint8_t* )calloc(num, 1);
	}

	for(r_i = 0;r_i < 2; ++r_i)
	{
		if (anchored_exon_num[r_i] == 0)
			continue;
		ref_begin[r_i] = anchored_exon[tid][r_i][0].ref_begin;
		ref_end[r_i] = anchored_exon[tid][r_i][0].ref_end;
		for(i = 1; i < anchored_exon_num[r_i]; ++i)
		{
			if (ref_begin[r_i] > anchored_exon[tid][r_i][i].ref_begin)
			{
				ref_begin[r_i] = anchored_exon[tid][r_i][i].ref_begin;
			}
			if (ref_end[r_i] < anchored_exon[tid][r_i][i].ref_end)
			{
				ref_end[r_i] = anchored_exon[tid][r_i][i].ref_end;
			}
		}
	}
	if (uniseed2_length[0] != 0)
	{
		dist_0 = creatGraph(uniseed3[tid][0], uniseed2_length[0], ref_begin[0], ref_end[0], path0, &max_index0, out_degree[0], tid);
	}
	else
	{
		dist_0 = 0;
	}

	if (uniseed2_length[1] != 0)
	{
		dist_1 = creatGraph(uniseed3[tid][1], uniseed2_length[1], ref_begin[1], ref_end[1], path1, &max_index1, out_degree[1], tid);
	}
	else
	{
		dist_1 = 0;
	}

	float dist_tmp; 
    uint8_t temp_strand; //1 - 0 +  2 +/-(both)
	uint32_t tmp_read_length;
	int thre = 30;
	if (dist_0 >= dist_1 && dist_0 > thre)
	{
		dist_tmp = dist_0 * 0.9;
		if (dist_1 > dist_tmp)
			temp_strand = 2;
		else
			temp_strand = 0;
	}
	else if(dist_0 < dist_1 && dist_1 > thre)
	{
		dist_tmp = dist_1 * 0.9;
		if(dist_0 > dist_tmp)
			temp_strand = 2;
		else
			temp_strand = 1;
	}
	else
	{
        if(path0 != NULL)	free(path0);
        if(path1 != NULL)	free(path1);
        for ( i = 0; i < 2; ++i)
        {
            if (out_degree[i] != NULL)	free(out_degree[i]);
        }
        if(out_degree != NULL)	free(out_degree);
		return 0;
	}
	// printf("temp_strand = %d\n", temp_strand);

	read_end = (uint32_t**)calloc(2, sizeof(uint32_t*));
	for (r_i = 0; r_i < 2; ++r_i)
	{
		read_end[r_i] = (uint32_t*)calloc(uniseed2_length[r_i], sizeof(uint32_t));
	}
	if (read_end == NULL)
	{
		fprintf(stderr, "memory wrong, read_end\n" );
	}

	seq_len = (int**)calloc(2,sizeof(int*));
	for (r_i = 0; r_i < 2; ++r_i)
	{
		seq_len[r_i] = (int*)calloc(uniseed2_length[r_i], sizeof(int));
	}
	if (seq_len == NULL)
	{
		fprintf(stderr, "memory wrong, seq_len\n" );
	}

	read_seq = (uint8_t***)calloc(2,sizeof(uint8_t**));

	if (read_seq == NULL)
	{
		fprintf(stderr, "memory wrong, read_seq\n" );
	}

    if (temp_strand == 1) //-
    {
    	cnt1 = get_split_read(max_index1, path1, read_end[1]);
    	if (cnt1 != 0)
    	{
    		read_seq[1] = (uint8_t** )calloc(cnt1, sizeof(uint8_t*));
    		for (i = 0; i < cnt1; ++i)
    		{
    			tmp_read_length = read_length_msa - read_end[1][i];
    			seq_len[1][i] = tmp_read_length;
    			read_seq[1][i] = (uint8_t*)calloc(tmp_read_length, sizeof(uint8_t));
    			get_read_seq(qseq0[1], read_seq[1][i], read_end[1][i], tmp_read_length);
    			read_length_msa = read_end[1][i];
    		}
    		#ifdef _DEBUG
    		fprintf(stderr, "[Phase-Info] Processing Partial order alignment, generating consensus sequence...\n");
    		#endif
    		fprintf(fp_hqr, ">%s\n", query_info[seqi].name);
    		seq_msa(cnt1, read_seq[1], seq_len[1]);
    	}
    	else
    	{
    		fprintf(fp_hqr, ">%s\n", query_info[seqi].name);
    		fprintf(fp_hqr, "%s\n", query_info[seqi].read_seq);
    	}
    }
	else if (temp_strand == 0)
    {
    	cnt0 = get_split_read(max_index0, path0, read_end[0]);
    	if (cnt0 != 0)
    	{
    		read_seq[0] = (uint8_t** )calloc(cnt0, sizeof(uint8_t*));
    		for (i = 0; i < cnt0; ++i)
    		{
    			tmp_read_length = read_length_msa - read_end[0][i];
    			seq_len[0][i] = tmp_read_length;
    			read_seq[0][i] = (uint8_t*)calloc(tmp_read_length, sizeof(uint8_t));
    			get_read_seq(qseq0[0], read_seq[0][i], read_end[0][i], tmp_read_length);
    			read_length_msa = read_end[0][i];
    		}
    		#ifdef _DEBUG
    		fprintf(stderr, "[Phase-Info] Processing Partial order alignment, generating consensus sequence...\n");
    		#endif
    		fprintf(fp_hqr, ">%s\n", query_info[seqi].name);
    		seq_msa(cnt0, read_seq[0], seq_len[0]);
    	}
    	else
    	{
    		fprintf(fp_hqr, ">%s\n", query_info[seqi].name);
    		fprintf(fp_hqr, "%s\n", query_info[seqi].read_seq);
    	}
    }
	else  //both strand
	{
    	cnt1 = get_split_read(max_index1, path1, read_end[1]);
    	if (cnt1 != 0)
    	{
    		read_seq[1] = (uint8_t** )calloc(cnt1, sizeof(uint8_t*));
    		for (i = 0; i < cnt1; ++i)
    		{
    			tmp_read_length = read_length_msa - read_end[1][i];
    			seq_len[1][i] = tmp_read_length;
    			read_seq[1][i] = (uint8_t*)calloc(tmp_read_length, sizeof(uint8_t));
    			get_read_seq(qseq0[1], read_seq[1][i], read_end[1][i], tmp_read_length);
    			read_length_msa = read_end[1][i];
    		}
    		#ifdef _DEBUG
    		fprintf(stderr, "[Phase-Info] Processing Partial order alignment, generating consensus sequence...\n");
    		#endif
    		fprintf(fp_hqr, ">%s\n", query_info[seqi].name);
    		seq_msa(cnt1, read_seq[1], seq_len[1]);
    	}
    	else
    	{
    		fprintf(fp_hqr, ">%s\n", query_info[seqi].name);
    		fprintf(fp_hqr, "%s\n", query_info[seqi].read_seq);
    	}

    	read_length_msa = read_length;

    	cnt0 = get_split_read(max_index0, path0, read_end[0]);
    	if (cnt0 != 0)
    	{
    		read_seq[0] = (uint8_t** )calloc(cnt0, sizeof(uint8_t*));
    		for (i = 0; i < cnt0; ++i)
    		{
    			tmp_read_length = read_length - read_end[0][i];
    			seq_len[0][i] = tmp_read_length;
    			read_seq[0][i] = (uint8_t*)calloc(tmp_read_length, sizeof(uint8_t));
    			get_read_seq(qseq0[0], read_seq[0][i], read_end[0][i], tmp_read_length);
    			read_length = read_end[0][i];
    		}
    		#ifdef _DEBUG
    		fprintf(stderr, "[Phase-Info] Processing Partial order alignment, generating consensus sequence...\n");
    		#endif
    		fprintf(fp_hqr, ">%s\n", query_info[seqi].name);
    		seq_msa(cnt0, read_seq[0], seq_len[0]);
    	}
    	else
    	{
    		fprintf(fp_hqr, ">%s\n", query_info[seqi].name);
    		fprintf(fp_hqr, "%s\n", query_info[seqi].read_seq);
    	}
	}

	if(path0 != NULL)	free(path0);
	if(path1 != NULL)	free(path1);
    for ( i = 0; i < 2; ++i)
	{
		if (out_degree[i] != NULL)	free(out_degree[i]);
	}
	if(out_degree != NULL)	free(out_degree);

	for(r_i = 0;r_i < 2;++r_i)
	{
		if (qseq0[r_i] != NULL) free(qseq0[r_i]);
	}

	for(r_i = 0;r_i < 2;++r_i)
	{
		if (read_end[r_i] != NULL) free(read_end[r_i]);
	}
	if (read_end != NULL) free(read_end);


	for(r_i = 0;r_i < 2;++r_i)
	{
		if (seq_len[r_i] != NULL) free(seq_len[r_i]);
	}
	if (seq_len != NULL) free(seq_len);

	for (r_i = 0; r_i < cnt1; ++r_i)
    {
    	if (read_seq[1][r_i] != NULL) free(read_seq[1][r_i]);
    }
	for (r_i = 0; r_i < cnt0; ++r_i)
    {
    	if (read_seq[0][r_i] != NULL) free(read_seq[0][r_i]);
    }

    if (read_seq != NULL) free(read_seq);

    return 0;
}

int change_hit_to_HCAE(uint64_t (*read_bit)[((MAX_READLEN - 1) >> 5) + 1], uint32_t read_length, uint32_t seqi, uint8_t tid)
{
	uint32_t read_off = 0;
	uint32_t r_b_v = 0;
	uint8_t re_d = 0;
	uint64_t kmer_bit = 0;
	uint32_t seed_hash = 0;
	uint32_t seed_kmer = 0;
	int64_t seed_id_r = 0;
	uint32_t ref_pos_n = 0;
	uint32_t uni_offset_s_l = 0;
	uint32_t uni_offset_s_r = 0;

	uint32_t left_i = 1;
	uint32_t right_i = 1;
	uint32_t k_way;
	uint64_t hit_i;
	uint64_t hit_binary0, hit_binary1;
	uint32_t mem_i;

	uint32_t max_right_i = 0;
	uint32_t read_pos = 0;
	uint32_t mem_length = 0;
	int64_t seed_binary[2] = {0,0};
	uint32_t memid[2] = {0,0};
	uint32_t uniseed_length[2] = {0,0}; 
	uint32_t anchored_exon_num[2] = {0,0};

	int result = 0;	
	uint8_t r_i = 0;

	uint32_t su_i = 0;
	uint64_t uni_id_temp;
	int32_t j;

	int t = (k_first_level > seed_k_t)? (k_first_level - seed_k_t) : 0;
	uint64_t kmer_pos_uni = 0;

	// fprintf(fp_tff, "uid\tseed_id\tread_pos\tuni_pos_off\tlength\tpos_n\n");
	#ifdef _DEBUG
	fprintf(fp_ae, "id\tref_begin\tref_end\tcov\tlength\n");
	fprintf(fp_tfu, "seed_id\tread_begin\tread_end\tref_begin\tref_end\tcov\n");
	//seed
	fprintf(stderr, "[Process-Info] Using seed-and-extension to find hits between read and reference...\n");
	#endif

	for ( r_i = 0; r_i < 2; ++r_i)//r_i=0, 1
	{
		// fprintf(fp_tff, "strand=%d\n",r_i);
		k_way = 1;
		mem_i = 0;
		r_b_v = 0;
        read_off = 0;
		for(read_off = 0; read_off <= read_length - seed_k_t; read_off += seed_step) //every seed_l we choose a seed
		{
			if(read_off + seed_k_t - 1 <= r_b_v)
			{
				continue;
			}

			re_d = (read_off & 0X1f);/////re_d = read_off & 00011111 = read_off while read_off < 31 , 32 loop
			//printf("read_off=%d, re_d=%d\n",read_off,re_d);
			if(re_d <= re_b)  //re_b = 32 - seed_k_t = 17
			{
				kmer_bit = ((read_bit[r_i][read_off >> 5] & bit_tran_re[re_d]) >> ((re_b - re_d) << 1));
				//printf("kmer_bit1=%lx\n",kmer_bit);
			}
			else
			{
				kmer_bit = (((read_bit[r_i][read_off >> 5] & bit_tran_re[re_d]) << ((re_d - re_b) << 1)) | (read_bit[r_i][(read_off >> 5) + 1] >> ((re_2bt - re_d) << 1)));
				//printf("kmer_bit2=%lx\n",kmer_bit);
			}
			if (seed_k_t == k_first_level) //k = 14
			{
				seed_hash = kmer_bit; //
//				printf("seed_hash1=%d\n",seed_hash);
				seed_binary[0] = 0;
				seed_binary[0] += buffer_hash_g[seed_hash];
				seed_binary[1] = 0;
                seed_binary[1] += buffer_hash_g[seed_hash + 1] - 1;
				if (seed_binary[1] < seed_binary[0])
					continue;
			}
			else
			{
				seed_kmer = (kmer_bit & bit_tran[k_r]); 
//				printf("seed_kmer=%d\n",seed_kmer);
				seed_hash = (kmer_bit >> (k_r << 1));
//				printf("seed_hash2=%d\n",seed_hash);
				result = binsearch_range(seed_kmer, buffer_kmer_g + buffer_hash_g[seed_hash], buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], seed_binary, seed_offset<<1);

				if (result == -1)
				{
					continue;
				}
				seed_binary[0] += buffer_hash_g[seed_hash];
				seed_binary[1] += buffer_hash_g[seed_hash];
			}
			// printf("seed_binary[0]=%ld\n",seed_binary[0]);
			// printf("seed_binary[1]=%ld\n",seed_binary[1]);

			max_right_i = 1;
			hit_binary0 = seed_binary[0];
			hit_binary1 = seed_binary[1];

			if ((hit_binary1 - hit_binary0 + 1) > uni_pos_n_max)
			{
				continue;
			}
			if (mem_i + hit_binary1 - hit_binary0 + 1 >= new_seed_cnt)
			{
				break;
			}
            //extension
            int ref_cnt = 0;
			for (hit_i = hit_binary0; hit_i <= hit_binary1; ++hit_i)
			{
				kmer_pos_uni = buffer_off_g[hit_i];//this kmer's offset on unipath seq
				//find the UID of this kmer
				seed_id_r = binsearch_interval_unipath64(kmer_pos_uni, buffer_seqf, result_seqf);

				ref_pos_n = buffer_pp[seed_id_r + 1] - buffer_pp[seed_id_r];
                ref_cnt += ref_pos_n;

				//if (ref_cnt > POS_N_MAX)
                if(ref_pos_n > pos_n_max)
				{
                    continue;
					//break;
				}
				uni_offset_s_l = kmer_pos_uni - buffer_seqf[seed_id_r];
				uni_offset_s_r = buffer_seqf[seed_id_r + 1] - (kmer_pos_uni + seed_k_t);

				// printf("uni_offset_s_l=%d\n",uni_offset_s_l);
				// printf("uni_offset_s_r=%d\n",uni_offset_s_r);
				//extend the kmer to a exact match 
 				for(left_i = 1; (left_i <= uni_offset_s_l) && (left_i <= read_off); left_i++)
 				{
					if(((buffer_seq[(kmer_pos_uni - left_i) >> 5] >> ((31 - ((kmer_pos_uni - left_i) & 0X1f)) << 1)) & 0X3)
					        != ((read_bit[r_i][(read_off - left_i) >> 5] >> ((31 - ((read_off - left_i) & 0X1f)) << 1)) & 0X3)
					        ) break;
				}

				for(right_i = 1; (right_i <= uni_offset_s_r) && (right_i <= read_length - read_off - seed_k_t); right_i++)
				{
					if(((buffer_seq[(kmer_pos_uni + seed_k_t - 1 + right_i) >> 5] >> ((31 - ((kmer_pos_uni + seed_k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
					        != ((read_bit[r_i][(read_off + seed_k_t - 1 + right_i) >> 5] >> ((31 - ((read_off + seed_k_t - 1 + right_i) & 0X1f)) << 1)) & 0X3)
					        ) break;
 				}

				read_pos = read_off + 1 - left_i;
				mem_length = seed_k_t + left_i + right_i - 2;

				hitseed[tid][r_i][mem_i].uid = seed_id_r;
				hitseed[tid][r_i][mem_i].seed_id = k_way;
				hitseed[tid][r_i][mem_i].read_pos = read_pos;
				hitseed[tid][r_i][mem_i].uni_pos_off = uni_offset_s_l + 1 - left_i;
				hitseed[tid][r_i][mem_i].length = mem_length;
				hitseed[tid][r_i][mem_i].pos_n = ref_pos_n;
				// fprintf(fp_tff, "%ld\t%d\t%d\t%d\t%d\t%d\n", hitseed[r_i][mem_i].uid,hitseed[r_i][mem_i].seed_id,hitseed[r_i][mem_i].read_pos,hitseed[r_i][mem_i].uni_pos_off,hitseed[r_i][mem_i].length,hitseed[r_i][mem_i].pos_n);
				if (right_i > max_right_i)
				{
					max_right_i = right_i;
				}
				
				++mem_i;
			}

			k_way++;
			r_b_v = read_off + seed_k_t + max_right_i - 1;
		}

		//merge colinear seeds in the same unipath
		su_i = 0;
		uint32_t s1, e1;
		if (mem_i == 0)
		{
			memid[r_i] = 0;
			continue;
		}
		if (mem_i == 1)
		{
			memid[r_i] = mem_i;
			hitseed2[tid][r_i][0].uid = hitseed[tid][r_i][0].uid;
			hitseed2[tid][r_i][0].read_pos = hitseed[tid][r_i][0].read_pos;
			hitseed2[tid][r_i][0].uni_pos_off = hitseed[tid][r_i][0].uni_pos_off; //whether to set uni_pos_off to the leftest position
			hitseed2[tid][r_i][0].pos_n = hitseed[tid][r_i][0].pos_n;
			hitseed2[tid][r_i][0].length1 = hitseed[tid][r_i][0].length;
			hitseed2[tid][r_i][0].length2 = hitseed[tid][r_i][0].length;
			hitseed2[tid][r_i][0].cov = hitseed[tid][r_i][0].length;
			hitseed2[tid][r_i][0].num = 1;
			continue;
		}
		qsort(hitseed[tid][r_i], mem_i, sizeof(hit_seed), compare_uniid);

		uni_id_temp = hitseed[tid][r_i][0].uid;
		j = 0;
		uint32_t cov = 0;
		while (j < mem_i)
		{
			s1 = j;
			cov = hitseed[tid][r_i][s1].length;
			j++;
			while ((uni_id_temp == hitseed[tid][r_i][j].uid) && (hitseed[tid][r_i][j].uni_pos_off > hitseed[tid][r_i][j-1].uni_pos_off) && (j < mem_i))
			{
				int diff = (int)(hitseed[tid][r_i][j].read_pos - hitseed[tid][r_i][j-1].read_pos - hitseed[tid][r_i][j-1].length);
				if (diff > waitingLen)
					break;
				if (abs((hitseed[tid][r_i][j].uni_pos_off - hitseed[tid][r_i][j-1].uni_pos_off) - (hitseed[tid][r_i][j].read_pos - hitseed[tid][r_i][j-1].read_pos)) < Eindel)
				{
					cov += (diff > 0)? hitseed[tid][r_i][j].length : (diff + hitseed[tid][r_i][j].length);
					++j;
				}
				else
					break;
			}
			e1 = j - 1;
			hitseed2[tid][r_i][su_i].uid = hitseed[tid][r_i][s1].uid;
			hitseed2[tid][r_i][su_i].read_pos = hitseed[tid][r_i][s1].read_pos;
			hitseed2[tid][r_i][su_i].uni_pos_off = hitseed[tid][r_i][s1].uni_pos_off; //whether to set uni_pos_off to the leftest position
			hitseed2[tid][r_i][su_i].pos_n = hitseed[tid][r_i][s1].pos_n;
			hitseed2[tid][r_i][su_i].cov = cov;
			hitseed2[tid][r_i][su_i].num = e1 - s1 + 1;
			cov = 0;
			
			//cal length
			if (s1 == e1)
			{
				hitseed2[tid][r_i][su_i].length1 = hitseed[tid][r_i][s1].length;
				hitseed2[tid][r_i][su_i].length2 = hitseed[tid][r_i][s1].length;
			}else{
				hitseed2[tid][r_i][su_i].length1 = hitseed[tid][r_i][e1].read_pos + hitseed[tid][r_i][e1].length - hitseed[tid][r_i][s1].read_pos;
				hitseed2[tid][r_i][su_i].length2 = hitseed[tid][r_i][e1].uni_pos_off + hitseed[tid][r_i][e1].length - hitseed[tid][r_i][s1].uni_pos_off;
			}
			
			uni_id_temp = hitseed[tid][r_i][j].uid;
			++su_i;
		}
		memid[r_i] = su_i;
	}

    if ((memid[0] == 0) && (memid[1] == 0))
    {
		return 0;
    }

    #ifdef _DEBUG
	fprintf(stderr, "[Process-Info] Mapping hits to reference genome...\n");
	#endif
	map_to_ref(hitseed2[tid], memid, uniseed_length, tid);

	#ifdef _DEBUG
	fprintf(stderr, "[Process-Info] Finding anchored exon region...\n");
	#endif
	find_anchored_exon(uniseed[tid], uniseed_length, anchored_exon_num, tid);

	find_hit_in_HCAE(seqi, read_length, tid, anchored_exon_num);

	return 0;
}

int seeding_core(int read_seq_core, uint8_t tid)
{
	
	uint32_t read_length = 0;
	uint32_t read_length_a = 0;
	uint16_t read_bit_char = 0;

	uint32_t r_i = 0;
	uint32_t seqi = 0;
	uint8_t rc_i;
	uint8_t c_tmp = 0;
	char tmp_char;
	
	seqi = read_seq_core;
	read_length = query_info[seqi].read_length;
	#ifdef _OUTRE
		fprintf(fp_re, "%d\t", read_length);
	#endif
	#ifdef _DEBUG
	fprintf(stderr,"[Process-Info] processing the %d-th read\n", seqi);
	fprintf(stderr,"[Process-Info] The %d-th read, name: %s, length: %d\n",seqi,query_info[seqi].name,read_length); 
	#endif
    readlen_max = (read_length > readlen_max)? read_length : readlen_max;
    
    read_len[tid] += read_length;

	if (read_length < 30)
	{
		return 1;
	}

	read_length_a = read_length - 1;
	read_bit_char = (((uint16_t )((read_length_a >> 5) + 1)) << 3);

	memset(read_bit1[tid][0], 0, read_bit_char);
	memset(read_bit1[tid][1], 0, read_bit_char);

	r_i = 0;
	while ((query_info[seqi].read_seq)[r_i])
	{
		tmp_char = (query_info[seqi].read_seq)[r_i];
		if (tmp_char == 'N')
		{
			//random
			tmp_char = "ACGT"[rand()%4];
		}
		c_tmp = charToDna5n[(uint8_t)tmp_char];
		read_bit1[tid][0][r_i >> 5] |= (((uint64_t )c_tmp) << ((31 - (r_i & 0X1f)) << 1));
		read_bit1[tid][1][(read_length_a - r_i) >> 5] |= (((uint64_t )(c_tmp ^ 0X3)) << ((31 - ((read_length_a - r_i) & 0X1f)) << 1));
		r_i++;
	}

	change_hit_to_HCAE(read_bit1[tid], read_length, seqi, tid);

	return 0;
}

static void *seeding_core_thread(void *aux)
{
	thread_aln_t *d = (thread_aln_t* )aux;
	int _read_lines;

	while(1)
	{
		pthread_rwlock_wrlock(&rwlock);
		_read_lines = THREAD_READ_I++;
		pthread_rwlock_unlock(&rwlock);

		if (_read_lines < d->seqn)
		{
			// fprintf(stderr, "the %d-th read\n", _read_lines);
			#ifdef _OUTRE
				fprintf(fp_re, "%d\t", _read_lines);
			#endif
			seeding_core(_read_lines, d->tid);
		}
		else break;
	}
	return 0;
}

int load_fasta_1pass(bseq_file_t *bf)
{
	
	uint32_t read_in = batch_size;
	uint32_t seqii = read_in;
	uint32_t r_i = 0;
	uint32_t r_ii = 0;
	int32_t r_iii = 0;
	uint32_t primary;

	k_r = seed_k_t - k_first_level;//k_first_level = 14 k_r = 1
    re_b = 32 - seed_k_t;//re_b = 17
    re_bt = (re_b << 1);//re_bt = 34
	//printf("re_bt=%d\n",re_bt);
	re_2bt = 64 - seed_k_t;//re_2bt = 49

	#ifdef _DEBUG
	fp_ae = fopen(temp_ae_dir, "w");
	if (fp_ae == NULL)
	{
		fprintf(stderr, "[Wrong] Failed to open file %s!!!\n", temp_ae_dir);
		exit(0);
	}

	fp_tfu = fopen(temp_uni_dir, "w");
	if (fp_tfu == NULL)
	{
		fprintf(stderr, "[Wrong] Failed to open file %s!!!\n", temp_uni_dir);
		exit(0);
	}

	fp_tfu2 = fopen(temp_uni2_dir, "w");
	if (fp_tfu2 == NULL)
	{
		fprintf(stderr, "[Wrong] Failed to open file %s!!!\n", temp_uni2_dir);
		exit(0);
	}
	fprintf(fp_tfu2, "id\tread_begin\tread_end\tref_begin\tref_end\n");
	#endif
	#ifdef _OUTRE
	fp_re = fopen(temp_re_dir, "w");
	if (fp_re == NULL)
	{
		fprintf(stderr, "[Wrong] Failed to open file %s!!!\n", temp_re_dir);
		exit(0);
	}
	#endif

	fp_hqr = fopen(temp_hqr_dir, "w");
	if (fp_hqr == NULL)
	{
		fprintf(stderr, "[Wrong] Failed to open file %s!!!\n", temp_hqr_dir);
		exit(0);
	}

	query_info = (READ_t* )calloc(read_in, sizeof(READ_t));

	hitseed = (hit_seed*** )calloc(thread_n, sizeof(hit_seed** ));
	for ( r_i = 0; r_i < thread_n; ++r_i)
	{
		hitseed[r_i] = (hit_seed** )calloc(2, sizeof(hit_seed* ));
		for(r_ii = 0; r_ii < 2; ++r_ii)
		{
			hitseed[r_i][r_ii] = (hit_seed* )calloc(new_seed_cnt, sizeof(hit_seed));
		}
	}
	if (hitseed == NULL)
	{
		fprintf(stderr, "memory wrong, hitseed\n" );
	}

	hitseed2 = (co_hitseed*** )calloc(thread_n, sizeof(co_hitseed** ));
	for ( r_i = 0; r_i < thread_n; ++r_i)
	{
		hitseed2[r_i] = (co_hitseed** )calloc(2, sizeof(co_hitseed* ));
		for(r_ii = 0; r_ii < 2; ++r_ii)
		{
			hitseed2[r_i][r_ii] = (co_hitseed* )calloc(new_seed_cnt, sizeof(co_hitseed));
		}
	}
	if (hitseed2 == NULL)
	{
		fprintf(stderr, "memory wrong, hitseed2\n" );
	}

	anchored_exon = (anchored_exons*** )calloc(thread_n, sizeof(anchored_exons** ));
	for ( r_i = 0; r_i < thread_n; ++r_i)
	{
		anchored_exon[r_i] = (anchored_exons** )calloc(2, sizeof(anchored_exons* ));
		for(r_ii = 0; r_ii < 2; ++r_ii)
		{
			anchored_exon[r_i][r_ii] = (anchored_exons* )calloc(new_seed_cnt, sizeof(anchored_exons));
		}
	}
	if (anchored_exon == NULL)
	{
		fprintf(stderr, "memory wrong, anchored_exon\n" );
	}

	uniseed = (uni_seed*** )calloc(thread_n, sizeof(uni_seed** ));
	for ( r_i = 0; r_i < thread_n; ++r_i)
	{
		uniseed[r_i] = (uni_seed** )calloc(2, sizeof(uni_seed* ));
		for(r_ii = 0; r_ii < 2; ++r_ii)
		{
			uniseed[r_i][r_ii] = (uni_seed* )calloc(new_seed_cnt*pos_n_max, sizeof(uni_seed));
		}
	}
	if (uniseed == NULL)
	{
		fprintf(stderr, "memory wrong, uniseed\n" );
	}

	uniseed3 = (uni_seed***)calloc(thread_n, sizeof(uni_seed**));
	for ( r_i = 0; r_i < thread_n; ++r_i)
	{
    	uniseed3[r_i] = (uni_seed**)calloc(2, sizeof(uni_seed*));
    	for(r_ii = 0;r_ii < 2; ++r_ii)
    	{
        	uniseed3[r_i][r_ii] = (uni_seed* )calloc(10000 , sizeof(uni_seed)); //the threshold should be change
    	}
    }
    if (uniseed3 == NULL)
	{
		fprintf(stderr, "memory wrong, uniseed3\n" );
	}

	read_len = (int* )calloc(thread_n, sizeof(int));

   	pthread_rwlock_init(&rwlock, NULL);
	thread_aln_t* aux;
	aux = (thread_aln_t* )calloc(thread_n, sizeof(thread_aln_t));
    double t_s;

	while(seqii == read_in)
    {
    	t_s = realtime();
    	seqii = bseq_read(bf, read_in, query_info);

		fprintf(stderr,"[Process-Info] total read number = %d\n",seqii);

        TOTAL_READ_COUNTs += seqii;

		THREAD_READ_I = 0;
		if (thread_n <= 1)
		{
			uint32_t seqii_i;
			for(seqii_i = 0; seqii_i < seqii; seqii_i++)
			// for(seqii_i = 0; seqii_i < 5000; seqii_i++)
			// seqii_i = 1725;
			{
				#ifdef _OUTRE
				fprintf(fp_re, "%d\t", seqii_i);
				#endif
				#ifdef _DEBUG
				fprintf(stderr, "the %d-th read\n", seqii_i);
				#endif
				seeding_core(seqii_i,0);
			}
		}
		else
		{
			pthread_t* tid;
			pthread_attr_t attr;

			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

			tid = (pthread_t* )calloc(thread_n, sizeof(pthread_t));

			for(r_i = 0; r_i < thread_n; ++r_i)
			{
				aux[r_i].tid = r_i;
				aux[r_i].seqn = seqii;
				// aux[r_i].seqn = 5000;
				int res = pthread_create(&tid[r_i], &attr, seeding_core_thread, aux + r_i);
				if(res != 0)
				{
					fprintf(stderr, "create pthread error");
					exit(1);
				}
			}

			for(r_i = 0; r_i < thread_n; ++r_i)	pthread_join(tid[r_i], 0);
			free(tid);
		}

		//free
		for (r_i = 0; r_i < read_in; ++r_i)
		{
			if (query_info[r_i].name != NULL)
			{
				free(query_info[r_i].name);
				query_info[r_i].name = NULL;
			}
			if (query_info[r_i].read_seq != NULL)
			{
				free(query_info[r_i].read_seq);
				query_info[r_i].read_seq = NULL;
			}
		}

		int t_len = 0;
        for (r_i = 0; r_i < thread_n; r_i++)
        {
            t_len += read_len[r_i];
            read_len[r_i] = 0;
        }
        
        fprintf(stderr, "[Phase-Info] Processing %d reads, total %d bases.\n", seqii, t_len); 

		//seqii = 0;
    }
   	pthread_rwlock_destroy(&rwlock);
	free(aux);
	free(read_len);
	// free memory

	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		for (r_ii = 0; r_ii < 2; ++r_ii)
		{
			if (hitseed[r_i][r_ii] != NULL)	free (hitseed[r_i][r_ii]);
		}
		if (hitseed[r_i] != NULL)	free (hitseed[r_i]);
	}
	if (hitseed != NULL) free(hitseed);

	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		for (r_ii = 0; r_ii < 2; ++r_ii)
		{
			if (hitseed2[r_i][r_ii] != NULL)	free (hitseed2[r_i][r_ii]);
		}
		if (hitseed2[r_i] != NULL)	free (hitseed2[r_i]);
	}
	if (hitseed2 != NULL) free(hitseed2);

	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		for (r_ii = 0; r_ii < 2; ++r_ii)
		{
			if (anchored_exon[r_i][r_ii] != NULL)	free (anchored_exon[r_i][r_ii]);
		}
		if (anchored_exon[r_i] != NULL)	free (anchored_exon[r_i]);
	}
	if (anchored_exon != NULL) free(anchored_exon);

	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		for (r_ii = 0; r_ii < 2; ++r_ii)
		{
			if (uniseed[r_i][r_ii] != NULL)	free (uniseed[r_i][r_ii]);
		}
		if (uniseed[r_i] != NULL)	free (uniseed[r_i]);
	}
	if (uniseed != NULL) free(uniseed);

	for (r_i = 0; r_i < thread_n; ++r_i)
	{
		for(r_ii = 0;r_ii < 2; ++r_ii)
    	{
        	if (uniseed3[r_i][r_ii] != NULL) free(uniseed3[r_i][r_ii]);
    	}
    	if (uniseed3[r_i] != NULL) free(uniseed3[r_i]);
    }
    if (uniseed3 != NULL)    free(uniseed3);

	for (r_i = 0; r_i < read_in; ++r_i)
	{
		if (query_info[r_i].read_seq != NULL)	free(query_info[r_i].read_seq);
		if (query_info[r_i].name != NULL)	free(query_info[r_i].name);
	}
	if(query_info != NULL)	free(query_info);

	#ifdef _DEBUG
	fclose(fp_ae);
	fclose(fp_tfu);
	fclose(fp_tfu2);
	#endif
	#ifdef _OUTRE
	fclose(fp_re);
	#endif
	fclose(fp_hqr);

	return 0;
}

void init_memory(param_map *opt, char *index_dir)
{
	load_index_file(index_dir);
	initGraph();
	initHashTable(hash_kmer);
}

void del_deBGAmemory()
{
	//if (buffer_ref_seq)	free(buffer_ref_seq);   //free after the program finish
	free(buffer_seqf);
	free(buffer_seq);
	free(buffer_pp);
	free(buffer_p);
	free(buffer_hash_g);
	free(buffer_kmer_g);
	free(buffer_off_g);
}

void del_finalmemory()
{
	free(buffer_ref_seq);   //free after the program finish
}

static int aln_usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:\tde Brijn Graph-based 3rd RNA sequence alignment\n");
	fprintf(stderr, "Usage:\t\trrs aln [options] <index_route> <read.fa/fq>\n\n");

	fprintf(stderr, "Algorithm options:\n\n");
	fprintf(stderr, "    -k --index-kmer       [INT]    K-mer length of RdBG-index. [%u]\n", INDEX_KMER);
	fprintf(stderr, "    -l --seeding-lmer     [INT]    K-mer length of seeding process (no long than RdBG-index). [%u]\n", SEEDING_KMER);
	fprintf(stderr, "    -s --seed-step        [INT]    The interval of seeding. [%u]\n", SEED_STEP);
    fprintf(stderr, "    -B --batch-size       [INT]    The number of reads to be processed in one loop. [%u]\n", BATCH_SIZE);
    fprintf(stderr, "    -n --max-uni-pos      [INT]    Maximum allowed number of hits per seed. [%u]\n", MAX_UNI_POS);
	fprintf(stderr, "    -i --min-frag-dis     [INT]    Maximum allowed distance of two fragment can be connect. [%u]\n", MIN_FRAG_DIS);
    fprintf(stderr, "    -x --read-type        [STR]    Specifiy the type of reads and set multiple paramters unless overriden.\n");
	fprintf(stderr, "                                   [null] default parameters.\n");
	fprintf(stderr, "                                   ccs (PacBio SMRT CCS reads): error rate 1%%\n");
	fprintf(stderr, "                                   clr (PacBio SMRT CLR reads): error rate 15%%\n");
	fprintf(stderr, "                                   ont1d (Oxford Nanopore 1D reads): error rate > 20%%\n");
	fprintf(stderr, "                                   ont2d (Oxford Nanopore 2D reads): error rate > 12%%\n\n");
	fprintf(stderr, "    -a --local-hash-kmer  [INT]    K-mer length of local hash process. [%u]\n", LOCAL_HASH_KMER);
	fprintf(stderr, "    -L --max-readlen      [INT]    Maximum allowed read length. [%u]\n", MAX_READLEN);
	fprintf(stderr, "    -g --max-read-gap     [INT]    Maximum allowed gap in read when chaining. [%u]\n", MAX_READ_JOIN_GAP);
	fprintf(stderr, "    -t --thread           [INT]    Number of threads. [4]\n");
	
	fprintf(stderr, "Output options\n\n");
	fprintf(stderr, "    -f --temp-file-perfix [STR]    Route of temporary files after the first-pass alignment. [%s]\n", TEMP_FILE_PERFIRX);
	return 1;
}

int help_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	rrs (Third generation RNA sequence alignment)\n");

	fprintf(stderr, "Usage:		rrs <command> [options]\n\n");
	fprintf(stderr, "Command: \n");
	fprintf(stderr, "		index		index reference sequence\n");
	fprintf(stderr, "		aln		align long RNA sequence to reference\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:	rrs index <ref.fa> <index_route>\n");
	fprintf(stderr, "		build deBGA index file using default k-mer length of deBGA. You can get more deBGA information from https://github.com/HongzheGuo/deBGA");
	fprintf(stderr, "\n\n");
	fprintf(stderr, "Usage:	rrs aln [options] <index_route> <read.fa/fq>\n\n");

	fprintf(stderr, "Algorithm options:\n\n");
	fprintf(stderr, "    -k --index-kmer       [INT]    K-mer length of RdBG-index. [%u]\n", INDEX_KMER);
	fprintf(stderr, "    -l --seeding-lmer     [INT]    K-mer length of seeding process (no long than RdBG-index). [%u]\n", SEEDING_KMER);
	fprintf(stderr, "    -s --seed-step        [INT]    The interval of seeding. [%u]\n", SEED_STEP);
    fprintf(stderr, "    -B --batch-size       [INT]    The number of reads to be processed in one loop. [%u]\n", BATCH_SIZE);
	fprintf(stderr, "    -n --max-uni-pos      [INT]    Maximum allowed number of hits per seed. [%u]\n", MAX_UNI_POS);
	fprintf(stderr, "    -i --min-frag-dis     [INT]    Maximum allowed distance of two fragment can be connected. [%u]\n", MIN_FRAG_DIS);
    fprintf(stderr, "    -x --read-type        [STR]    Specify the type of reads and set multiple paramters unless overriden.\n");
	fprintf(stderr, "                                   [null] default parameters.\n");
	fprintf(stderr, "                                   ccs (PacBio SMRT CCS reads): error rate 1%%\n");
	fprintf(stderr, "                                   clr (PacBio SMRT CLR reads): error rate 15%%\n");
	fprintf(stderr, "                                   ont1d (Oxford Nanopore 1D reads): error rate > 20%%\n");
	fprintf(stderr, "                                   ont2d (Oxford Nanopore 2D reads): error rate > 12%%\n\n");
	fprintf(stderr, "    -a --local-hash-kmer  [INT]    K-mer length of local hash process. [%u]\n", LOCAL_HASH_KMER);
	fprintf(stderr, "    -L --max-readlen      [INT]    Maximum allowed read length. [%u]\n", MAX_READLEN);
	fprintf(stderr, "    -g --max-read-gap     [INT]    Maximum allowed gap in read when chaining. [%u]\n", MAX_READ_JOIN_GAP);
	fprintf(stderr, "    -t --thread           [INT]    Number of threads. [4]\n");
	
	fprintf(stderr, "Output options\n\n");
	fprintf(stderr, "    -f --temp-file-perfix [STR]    Route of temporary files after the first-pass alignment. [%s]\n", TEMP_FILE_PERFIRX);
	return 1;
}

static const char *short_option = "k:l:B:n:f:s:i:x:a:L:g:t:";

static struct option long_option[] = {
	{"index-kmer", required_argument, NULL, 'k'},
	{"seeding-lmer", required_argument, NULL, 'l'},
	{"seed-step", required_argument, NULL, 's'},
	{"batch-size", required_argument, NULL, 'B'},
    {"max-uni-pos", required_argument, NULL, 'n'},
	{"min-frag-dis", required_argument, NULL, 'i'},
	{"temp-file-perfix", required_argument, NULL, 'f'},
	{"read-type", required_argument, NULL, 'x'},
	{"local-hash-kmer", required_argument, NULL, 'a'},
	{"max-readlen", required_argument, NULL, 'L'},
	{"max-read-gap", required_argument, NULL, 'g'},
	{"thread", required_argument, NULL, 't'},
	{0,0,0,0}
};

int rrs_par(int argc, char *argv[], const char *version)
{ 
	fprintf(stderr, "[Main] rrs - Partitioning Long Transcriptome Reads by Rolling Circle Amplyfication Sequencing\n");
	param_map *opt = (param_map* )calloc(1, sizeof(param_map));
	init_map_param(opt);
	int c;
	char *p;

    sprintf(command, "@PGrrs\tPN:rrs\tVN:%s\tCL:%s", version, argv[0]);
//      printf("1%s\n",command);
    for (c = 1; c < argc; ++c) sprintf(command+strlen(command), " %s", argv[c]);
//printf("2%ld\n",strlen(command));
	while((c = getopt_long(argc, argv, short_option, long_option, NULL)) != -1)
	{
		switch(c)
		{
			case 'k': opt->k_t = atoi(optarg); break;
			case 'l': opt->seed_k_t = atoi(optarg); break;
			case 's': opt->seed_step = atoi(optarg); break;
            case 'B': opt->batch_size = atoi(optarg); break;
			case 'n': opt->pos_n_max = atoi(optarg); break;
			case 'i': opt->Eindel = atoi(optarg); break;
            case 'f': opt->temp_file_perfix = strdup(optarg); break;
            case 'a': opt->hash_kmer = atoi(optarg); break;
            case 'L': opt->readlen_max = atoi(optarg); break;
            case 'g': opt->max_read_join_gap = atoi(optarg); break;
			case 't': opt->thread_n = atoi(optarg); break;
			case 'x': if (strcmp(optarg, "ccs") == 0) opt->read_type = 1;
					  else if (strcmp(optarg, "clr") == 0) opt->read_type = 2;
					  else if (strcmp(optarg, "ont1d") == 0) opt->read_type = 3;
					  else if (strcmp(optarg, "ont2d") == 0) opt->read_type = 4;
					  else {
						  fprintf(stderr, "[main:] Unkown parameter: %s\n", optarg);
						  return aln_usage();
					  } 
			default: return aln_usage(); break;
		}
	} 

    if (argc - optind < 3)
	return aln_usage();
    if ((opt->seed_k_t < 14) || (opt->seed_k_t > opt->k_t))
    {
        fprintf(stderr, "Input error: -k cannot be less than 14 or more than %d\n", opt->k_t);
        exit(1);
    }
    if ((opt->seed_step < 1) || (opt->seed_step > 10))
    {
        fprintf(stderr, "Input error: -s cannot be less than 1 or more than 10\n");
        exit(1);
    }
    if ((opt->pos_n_max < 30) || (opt->pos_n_max > 5000))
    {
        fprintf(stderr, "Input error: -n cannot be less than 30 or more than 5000\n");
        exit(1);
    }
    if ((opt->thread_n < 1) || (opt->thread_n > 48))
    {
        fprintf(stderr, "Input error: -t cannot be less than 1 or more than 48\n");
        exit(1);
    }
    if ((opt->Eindel < 15) || (opt->Eindel > 115))
    {
        fprintf(stderr, "Input error: -i cannot be less than 15 or more than 30\n");
        exit(1);
    }
    if ((opt->hash_kmer < 6) || (opt->hash_kmer > 10))
    {
        fprintf(stderr, "Input error: -a cannot be less than 6 or more than 10\n");
        exit(1);
    }
	char *index_dir;
	char *read_fastq;
	index_dir = strdup(argv[optind + 1]);
	if (index_dir[strlen(index_dir) - 1] != '/') strcat(index_dir, "/");
	read_fastq = strdup(argv[optind + 2]);

	#ifdef _DEBUG
	memset(temp_uni_dir, 0, 1024);
	memset(temp_uni2_dir, 0, 1024);
	memset(temp_ae_dir, 0, 1024);
    if (opt->temp_file_perfix == NULL)
    {
        strcpy(temp_uni_dir, "./uni.lines");
        strcpy(temp_ae_dir, "./ae.lines");
        strcpy(temp_uni2_dir, "./mb.lines");
    }
    else
    {
        strcpy(temp_uni_dir, opt->temp_file_perfix);
        strcat(temp_uni_dir, "uni.lines");
        strcpy(temp_ae_dir, opt->temp_file_perfix);
        strcat(temp_ae_dir, "ae.lines");
        strcpy(temp_uni2_dir, opt->temp_file_perfix);
        strcat(temp_uni2_dir, "mb.lines");
    }
	#endif
	#ifdef _OUTRE
	memset(temp_re_dir, 0, 1024);
    if (opt->temp_file_perfix == NULL)
    {
        strcpy(temp_re_dir, "./re.lines");
    }
    else
    {
        strcpy(temp_re_dir, opt->temp_file_perfix);
        strcat(temp_re_dir, "re.lines");
    }
	#endif
	memset(temp_hqr_dir, 0, 1024);
    if (opt->temp_file_perfix == NULL)
    {
        strcpy(temp_hqr_dir, "./hq_reads.fasta");
    }
    else
	{
        strcpy(temp_hqr_dir, opt->temp_file_perfix);
        strcat(temp_hqr_dir, "hq_reads.fasta");
    }

    //variable in this file
    thread_n = opt->thread_n;
	k_t = opt->k_t;
	seed_k_t = opt->seed_k_t;
	Eindel = opt->Eindel;
  	batch_size = opt->batch_size;
    read_type = opt->read_type;
	seed_step = opt->seed_step;
	pos_n_max = opt->pos_n_max;
	hash_kmer = opt->hash_kmer;
	readlen_max = opt->readlen_max;
	seed_offset = k_t - seed_k_t; 
	max_read_join_gap = opt->max_read_join_gap;
    seed_num = 10000; // extract at most 10000 seed every read

    // variable uni_pos_n_max and POS_N_MAX got from experience, which considering the speed and accuracy.
    //uni_pos_n_max
    if (seed_offset < 3)
        uni_pos_n_max = pow(4, seed_offset);
    else if (seed_offset < 6)
        uni_pos_n_max = pow(2, seed_offset);
    else
        uni_pos_n_max = 64; //64

    POS_N_MAX = 25;
    if ((seed_k_t == 15) || (seed_k_t == 16))
        POS_N_MAX = pos_n_max;
    else if (seed_k_t == 14)
        POS_N_MAX = 35;

    new_seed_cnt = seed_num * (uni_pos_n_max + 1);

    //waitlength
	float els = 0.05;
	float error = 0.2;

	if (read_type == 1) //ccs
		error = 0.02;
	else if (read_type == 2) //clr
		error = 0.15; 
	else if (read_type == 3) //ont1d
		error = 0.25;
	else if (read_type == 4)
		error = 0.13;

	float t = log10(els)/log10(1 - pow(1 - error, seed_k_t));
	float q = 1/error - seed_k_t * pow(1 - error, seed_k_t)/(1 - pow(1 - error, seed_k_t));
	waitingLen = (int)(t * q);
	// printf("waitingLen = %d \n", waitingLen);

    printf("read_fastq=%s\n",read_fastq);
    printf("index_dir=%s\n",index_dir);    
    bseq_file_t *bf;
    bf = bseq_open(read_fastq);
    if(bf == 0)
    {
        fprintf(stderr, "[Waring] Wrong input file route or name: %s \n", read_fastq);
        exit(1);
    }

	fprintf(stderr, "[Phase-Info] Loading Index and Reads\n");
	init_memory(opt, index_dir);

	fprintf(stderr, "[Phase-Info] Seeding and Chaining Phase\n");
    double tt1 = realtime();

    load_fasta_1pass(bf);
    fprintf(stderr, "[Phase-Info] Total %d reads were processed in %.3f seconds.\n", TOTAL_READ_COUNTs, realtime() - tt1);
   
	del_deBGAmemory();
	delGraph();
	freeHashTable();
	bseq_close(bf);
	del_finalmemory();

	free(opt);

	return 0;
}
