#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "graph.h"
#include "read_seeding.h"
#include "binarys_qsort.h"

// #define _DEBUG

double DAG_time = 0;
double findpath_time = 0;
double qsort_time = 0;

//global variant
Graph Dp_graph;

void initGraph()
{
	uint8_t r_i;

	Dp_graph.vnode = (VNode* )calloc(new_seed_cnt*pos_n_max,sizeof(VNode));
	if(Dp_graph.vnode == NULL) 
	{
		printf("wrong for allocate graph memory!\n");
	}
}

void delGraph()
{	
	uint8_t r_i;
	uint32_t r_ii;
	for (r_ii = 0; r_ii < new_seed_cnt*pos_n_max; ++r_ii)
	{
		if (Dp_graph.vnode[r_ii].preedge != NULL)	free(Dp_graph.vnode[r_ii].preedge);
	}
	if (Dp_graph.vnode != NULL)	free(Dp_graph.vnode);
}

static float min(float a, float b)
{
    return a<b?a:b;
}

void Find_longest_path(Graph *graph, uint32_t target, PATH_t *dist_path, uni_seed *vertexArr)
{
	float tmp_dist = 0;
    float penalty;
    float tmp_mean = 0, tmp = 0, diff = 0, mean = 0;
	float current_dist = 0;
	float mean_penalty;
	uint8_t reverse_edge;
	int32_t pre_node = 0;
    int32_t ver = 0;
    int32_t weight;
	uint32_t adjvex;
	uint32_t tmp_read_end, read_end = 0;
    //travel the processer of vertex target
    pANode arcnode;
	int i, param = 6;
	for (i = 0; i < graph->vnode[target].adjacent_node; ++i)
	{
		arcnode = &graph->vnode[target].preedge[i];
		weight = arcnode->weight;
		penalty = arcnode->penalty;
		adjvex = arcnode->adjvex;
		reverse_edge = arcnode->reverse_edge;

		mean_penalty = 0;

		if (reverse_edge == 0)
		{
			tmp_mean = dist_path[adjvex].mean;
			tmp_read_end = dist_path[adjvex].read_end;
		}
		else if (reverse_edge == 1)
		{
			if (dist_path[adjvex].mean == 0)//the first reverse edge
			{
				tmp_mean = vertexArr[adjvex].read_end;
				tmp_read_end = vertexArr[adjvex].read_end;
			}
			else
			{
				tmp_read_end = vertexArr[adjvex].read_end;//- the last reverse edge read end
				tmp = tmp_read_end - dist_path[adjvex].read_end;
				tmp_mean = (dist_path[adjvex].mean + tmp)/2;
				diff = abs(dist_path[adjvex].mean - tmp_mean);
				mean_penalty = (diff*param)/dist_path[adjvex].mean;//*
				// printf("mean_penalty = %f\n", mean_penalty);
			}
			//if there are two reverse edge, choose the longer one
			// if (tmp > tmp_mean)//record which reverse edge is chosen
			// {
				// tmp_mean = tmp;
			// }
		}


		tmp_dist = dist_path[adjvex].dist + weight - penalty - mean_penalty;//*

		if (current_dist <= tmp_dist)
        {
            current_dist = tmp_dist;
            pre_node = adjvex;
            mean = tmp_mean;
            read_end = tmp_read_end;
        }
	}

	dist_path[target].dist = current_dist;   // change the distance because there are overlaps between mems
	dist_path[target].pre_node = pre_node; //the front node
	dist_path[target].mean = mean;
	dist_path[target].read_end = read_end;
}

void dynamic_programming_path(Graph *graph, uint32_t vertexNum, PATH_t *dist_path, uint8_t *out_degree, uni_seed *vertexArr)
{
	uint32_t i;

	// for all For each v∈V\S in Topological Sorting order do dilg(v)=max(u,v)∈E{dilg(u)+w(u, v)}
	for (i = 0; i < vertexNum; ++i)
	{
		if (graph->vnode[i].adjacent_node > 0) //calculate dist[i]
		{
			Find_longest_path(graph, i, dist_path, vertexArr);
		}
	}

	// for (i = 0; i < vertexNum; ++i)
	// {
	// 	printf("vertex%d: %d, %d, %d, %d\n", i, vertexArr[i].read_begin,vertexArr[i].read_end,vertexArr[i].ref_begin,vertexArr[i].ref_end);
	// }
	#ifdef _DEBUG
	printf("show the path\n");
	int j;
	for (i = 0; i < vertexNum; ++i)
	{
		j = i;
		if ((out_degree[i] == 0))
		{
			printf("vertex: %d\tdist = %f\tmean = %f\tread_end = %d\n", i, dist_path[i].dist, dist_path[j].mean, dist_path[j].read_end);
			// printf("path: ");
			// printf("%d->", i);

			j = dist_path[j].pre_node;
			while(j != -1)
			{
				// printf("%d->", j);
				printf("vertex: %d\tdist = %f\tmean = %f\tread_end = %d\n", j, dist_path[j].dist, dist_path[j].mean, dist_path[j].read_end);
				j = dist_path[j].pre_node;
			}
			printf("\n");
		}
	}
	#endif

}

float creatGraph(uni_seed *vertexArr, uint32_t vertexNum, uint32_t ref_begin, uint32_t ref_end, PATH_t *dist_path, uint32_t *max_index, uint8_t *out_degree)
{
	int32_t i,j;

	int32_t weight = 0;
	uint32_t ref_pos2 = 0;
	uint32_t read_pos2 = 0;

	// float max_distance = 0;
    int32_t gap;
    int32_t ove1;
    int32_t ove2;
    int32_t ove3;
    int32_t dis1, dis2;
    int32_t ref_range;
    int8_t param = 2;
	int8_t intron_penalty = Eindel*param/seed_k_t;
	uint32_t non_ioslated_point= 0;
	uint32_t thre_num;
	int search_step = (vertexNum < 50)? vertexNum : 50;
	float max_distance = 0;
	
	Graph *graph = &Dp_graph;

	double time1 = clock();
	qsort(vertexArr, vertexNum, sizeof(uni_seed), compare_uniseed2);//by read start position
	qsort_time += (clock() - time1)/CLOCKS_PER_SEC;

	// re inital
	for (i = 0; i < vertexNum; ++i)
	{
        // vertexArr[i].id = i;
		graph->vnode[i].head_vertex = i;
		graph->vnode[i].adjacent_node = 0;
		graph->vnode[i].preedge = (ANode *)calloc(search_step, sizeof(ANode));

		dist_path[i].dist = weight;
		dist_path[i].pre_node = -1;
		dist_path[i].mean = 0;
		dist_path[i].read_end = 0;
	}


	// add for print dist
	memset(out_degree, 0, vertexNum);

	ref_range = ref_end - ref_begin;
	#ifdef _DEBUG
	printf("ref_range = %d\n",ref_range);
	#endif
	time1 = clock();
	for ( i = 0; i < vertexNum-1; ++i)
	{
		read_pos2 = vertexArr[i].read_end;
		ref_pos2 = vertexArr[i].ref_end;

		thre_num = (vertexNum < (i + search_step))? vertexNum : (i + search_step);

		for (j = i + 1; j < thre_num; ++j)
		{
			// if (vertexArr[j].ref_begin > ref_pos2 + max_intron_length)
			// {
			// 	break;
			// }
			//if two MB from the same seed, they can not be connected
			if (vertexArr[i].read_end == vertexArr[j].read_end || vertexArr[i].read_begin == vertexArr[j].read_begin)
				continue;
        
            ove1 = (int32_t)(vertexArr[j].read_begin - read_pos2);
            ove2 = (int32_t)(vertexArr[j].ref_begin - ref_pos2);
            ove3 = (int32_t)(vertexArr[i].ref_begin - vertexArr[j].ref_end);

            // if (ove1 > max_read_join_gap)
            // 	continue;
            

			if ((ove1 > 0 && ove2 >= 0 && ( ove2 > ove1/2 || abs(ove2 - ove1)==1 )) || (ove1 >= -8 && ove1 <= 0  && ove2 >= -8))  //5  1/error rate
            {
				//the first part is normal, if (ove1 - ove2) < Eindel, there is an edge
				//the last part, beaucse the begin of exon and the begin of intron have the same bases 
            	graph->vnode[j].preedge[graph->vnode[j].adjacent_node].adjvex = i;
            	graph->vnode[j].preedge[graph->vnode[j].adjacent_node].reverse_edge = 0;
                dis1 = vertexArr[j].read_end - read_pos2;
                dis2 = vertexArr[j].ref_end - ref_pos2;
                gap = (int32_t)(dis1 - dis2);
			
				int diff = (read_pos2 >= vertexArr[j].read_begin)? (read_pos2 + 1 - vertexArr[j].read_begin) : 0;
				weight = vertexArr[j].read_end - vertexArr[j].read_begin + 1 - diff;
				
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].weight = min(weight, min(dis1, dis2));
                if (gap >= 0) //in the same exon
					graph->vnode[j].preedge[graph->vnode[j].adjacent_node].penalty = gap*param/(float)weight;
                else
					graph->vnode[j].preedge[graph->vnode[j].adjacent_node].penalty = min(abs(gap)*param/(float)weight, intron_penalty);
                
				// if ((j == 35 || j == 36 || j == 37)&&(i == 35 || i == 36 || i == 37))
				// {
				// 	printf("(%d,%d)\n", i, j);
				// 	printf("penalty = %f\n", graph->vnode[j].preedge[graph->vnode[j].adjacent_node].penalty);
				// 	printf("weight = %d\n", graph->vnode[j].preedge[graph->vnode[j].adjacent_node].weight);
				// }
                graph->vnode[j].adjacent_node++;
				out_degree[i] = 1;

        		non_ioslated_point++;
            }
            else if (ove2 == ove1) //linear
            {
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].adjvex = i;
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].weight = vertexArr[j].read_end - vertexArr[j].read_begin + 1 - (read_pos2 + 1 - vertexArr[j].read_begin);
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].penalty = 0.0;
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].reverse_edge = 0;
				graph->vnode[j].adjacent_node++;
				out_degree[i] = 1;

				non_ioslated_point++;
            }
            else if (ove1 >= -8 && (ove3 > (1*ref_range)/2 || ove3 > ove1))
            // else if (ove1 >= -8 && ove3 > (1*ref_range)/2)
            // else if (ove1 >= -8 && ove3 > ove1)
            {
            	// printf("(%d,%d)\n", i, j);
            	// printf("ove1 = %d\n", ove1);
            	// printf("ove3 = %d\n", ove3);
            	// printf("vertexArr[%d].read_begin = %d\n",j,vertexArr[j].read_begin);
            	// printf("vertexArr[%d].read_begin = %d\n",i,vertexArr[i].read_begin);
            	graph->vnode[j].preedge[graph->vnode[j].adjacent_node].adjvex = i;
            	
            	dis1 = vertexArr[j].read_end - read_pos2;
				int diff = (read_pos2 >= vertexArr[j].read_begin)? (read_pos2 + 1 - vertexArr[j].read_begin) : 0;
				weight = vertexArr[j].read_end - vertexArr[j].read_begin + 1 - diff;
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].weight = min(weight, dis1);
				
				graph->vnode[j].preedge[graph->vnode[j].adjacent_node].penalty = (float)ref_range/(ove3*0.14);
				// printf("penalty = %f\n", graph->vnode[j].preedge[graph->vnode[j].adjacent_node].penalty);
            	graph->vnode[j].preedge[graph->vnode[j].adjacent_node].reverse_edge = 1;

    //         	if ((j == 35 || j == 36 || j == 37)&&(i == 35 || i == 36 || i == 37))
				// {
				// 	printf("(%d,%d)\n", i, j);
				// 	printf("penalty = %f\n", graph->vnode[j].preedge[graph->vnode[j].adjacent_node].penalty);
				// 	printf("weight = %d\n", graph->vnode[j].preedge[graph->vnode[j].adjacent_node].weight);
				// }
            	graph->vnode[j].adjacent_node++;

            	out_degree[i] = 1;
            	non_ioslated_point++;
            }
		}
	}
	DAG_time += (clock() - time1)/CLOCKS_PER_SEC;
	#ifdef _DEBUG
	fprintf(stderr, "[Process-Info] Generating Directed Acyclic Graph, total time %f seconds..\n",DAG_time);
	#endif
	if (non_ioslated_point != 0)
	{
		time1 = clock();
		#ifdef _DEBUG
		fprintf(stderr, "[Process-Info] Processing dynamic programming...");
		#endif
		dynamic_programming_path(graph, vertexNum, dist_path, out_degree, vertexArr);
		findpath_time += (clock() - time1)/CLOCKS_PER_SEC;
		#ifdef _DEBUG
		fprintf(stderr, "[Process-Info] Dynamic programming, total time %f seconds..\n",findpath_time);
		#endif
	}

	for ( i = 0; i < vertexNum; ++i)
	{
		if (max_distance < dist_path[i].dist)
		{
			max_distance = dist_path[i].dist;
			*max_index = i;
		}
	}

	//free
	for (i = 0; i < vertexNum; ++i)
	{
		free(graph->vnode[i].preedge);
		graph->vnode[i].preedge = NULL;
	}

	return max_distance;

}