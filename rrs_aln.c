#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<getopt.h>
#include<unistd.h>
#include<zlib.h>

int rrs_aln_usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:	rrs index <ref.fa> <index_route>\n");
	fprintf(stderr, "		build deBGA index file with default 22-kmer. You can get more deBGA information from https://github.com/HongzheGuo/deBGA");
	fprintf(stderr, "\n");
	return 1;
}

int get_bin_dir2(char *bin, char *dir)
{
	char *end = strrchr(bin, '/');
	if (end == NULL)
		return 1;

	bin[end-bin] = '\0';
	strcpy(dir, bin);
	return 0;
}

int deSALT_aln(char *dir, char *index_route, char *read_route)
{
	char cmd[1024];
	sprintf(cmd, "%sdeSALT aln %s %s", dir, index_route, read_route);
	fprintf(stderr, "[rrs_aln] Ececuting deSALT aln ...\n");
	if (system(cmd) != 0)
	{
		fprintf(stderr, "\n[rrs_aln] Aligning undoing, deSALT aln exit abnormally. \n"); 
		exit(1);
	}
	fprintf(stderr, "[rrs_aln] Done!\n");
}

int rrs_aln(int argc, char *argv[])
{
	// clock_t t = clock();
	char *index_route = 0;
	char *read_route = 0;
    
    char dir[1024];
    char rrs_path[1024];
    int r;

	if(optind + 3 > argc)	return rrs_aln_usage();

	index_route = strdup(argv[optind+1]);
	read_route = strdup(argv[optind+2]);
    
    r = readlink("/proc/self/exe", dir, 2048);
    if (r < 0 || r >= 2048)
    {
        fprintf(stderr, "Failed, could not find the program of rrs!\n");
    }
    dir[r] = '\0';

    if (!get_bin_dir2(dir, rrs_path))
    {
        strcat(rrs_path, "/");
    }

    char path_deSALT[1024];
    strcpy(path_deSALT, rrs_path);
    strcat(path_deSALT, "deSALT");

    if((access(path_deSALT, F_OK)) == -1)
    {
        fprintf(stderr, "[Wrong!] %s not exist, please check!\n", path_deSALT);
        exit(1);
    }

	deSALT_aln(rrs_path, index_route, read_route);

	return 0;
}
