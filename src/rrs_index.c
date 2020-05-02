#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<getopt.h>
#include<unistd.h>
#include<zlib.h>

int rrs_index_usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:	rrs index <ref.fa> <index_route>\n");
	fprintf(stderr, "		build deBGA index file with default 22-kmer. You can get more deBGA information from https://github.com/HongzheGuo/deBGA");
	fprintf(stderr, "\n");
	return 1;
}

int get_bin_dir(char *bin, char *dir)
{
	char *end = strrchr(bin, '/');
	if (end == NULL)
		return 1;

	bin[end-bin] = '\0';
	strcpy(dir, bin);
	return 0;
}

int deBGA_index(char *dir, char *ref_fa, char *index_route)
{
	char cmd[1024];
	sprintf(cmd, "%sdeBGA index %s %s", dir, ref_fa, index_route);
	fprintf(stderr, "[rrs_index] Ececuting deBGA index ...\n");
	if (system(cmd) != 0)
	{
		fprintf(stderr, "\n[rrs_index] Indexing undoing, deBGA index exit abnormally. \n"); 
		exit(1);
	}
	fprintf(stderr, "[rrs_index] Done!\n");
}

int rrs_index(int argc, char *argv[])
{
	// clock_t t = clock();
	char *ref_fa = 0;
	char *index_route = 0;
    
    char dir[1024];
    char rrs_path[1024];
    int r;

	if(optind + 3 > argc)	return rrs_index_usage();

	ref_fa = strdup(argv[optind+1]);
	index_route = strdup(argv[optind+2]);
    
    r = readlink("/proc/self/exe", dir, 2048);
    if (r < 0 || r >= 2048)
    {
        fprintf(stderr, "Failed, could not find the program of rrs!\n");
    }
    dir[r] = '\0';

    if (!get_bin_dir(dir, rrs_path))
    {
        strcat(rrs_path, "/");
    }

    char path_deBGA[1024];
    strcpy(path_deBGA, rrs_path);
    strcat(path_deBGA, "deBGA");

    if((access(path_deBGA, F_OK)) == -1)
    {
        fprintf(stderr, "[Wrong!] %s not exist, please check!\n", path_deBGA);
        exit(1);
    }

	deBGA_index(rrs_path, ref_fa, index_route);

	return 0;
}
