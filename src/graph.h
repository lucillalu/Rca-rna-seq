#ifndef _GRAPH_H
#define _GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "read_seeding.h"

typedef struct ArcNode
{
	uint32_t adjvex;
    int weight; 
    float penalty;
    uint8_t reverse_edge;
}ANode, *pANode;

typedef struct VertexNode
{
	uint32_t head_vertex;
	uint32_t adjacent_node;
    ANode* preedge;
}VNode;

typedef struct dpGraph
{
	VNode* vnode;
}Graph;

extern Graph* graph;

float creatGraph(uni_seed* vertexArr, uint32_t vertexNum, uint32_t ref_begin, uint32_t ref_end, PATH_t* dist_path, uint32_t *max_index, uint8_t *out_degree, uint8_t tid);
void initGraph();
void delGraph();
// void show_vertexm(vertex_m *vertexArr, uint32_t vertexNum);
// void show_vertexu(vertex_u *vertexArr, uint32_t vertexNum);
#endif
