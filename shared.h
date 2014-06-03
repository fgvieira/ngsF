#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <zlib.h>
#include <semaphore.h>
#include <sys/mman.h>
#include "zlib.h"
#include "bgzf/bgzf.h"

extern bool SIG_COND;


// Struct to store all input arguments //GZIP
typedef struct {
#ifdef _USE_BGZF
	BGZF* in_glf_fh;
#else
	FILE* in_glf_fh;
	double **data;
#endif
	char *in_glf;
	char *in_glf_type;
	char *init_values;
	char *out_file;
	uint16_t n_ind;
	uint64_t n_sites;
	uint64_t n_chunks;
	int64_t *chunks_voffset;
	uint64_t max_chunk_size;
	bool fast_lkl;
	bool approx_EM;
	bool call_geno;
	int max_iters;
	double min_epsilon;
	unsigned int n_threads;
        uint64_t seed;
	sem_t launch_thread_semaph;
	sem_t running_thread_semaph;
	pthread_mutex_t F_lock;
	bool quick;
	bool version;
	unsigned int verbose;
} params;


// Output data
typedef struct {
	double *site_freq;
	double *site_freq_num;
	double *site_freq_den;
	double *site_prob_var;
	double *site_tmpprob_var;
	double *indF;
	double *indF_num;
	double *indF_den;
	double *ind_lkl;
	double global_lkl;
} out_data;


// Prepare structure
typedef struct {
	params *pars;
	uint64_t chunk_size;
	uint64_t chunk_abs_start_pos;
	double **chunk_data;
	int iter;
	unsigned int *n_run_threads;
	out_data *output;
} pth_params;


// parse_args.cpp
void init_pars(params *);
int parse_cmd_args(int, char **, params *);

// read_data.cpp
int init_output(params *, out_data *);
uint64_t read_chunk(double **, params *, uint64_t);
void call_geno(double *, int, int);
int array_max_pos(double *, int);

// EM.cpp
int do_EM(params *, out_data *);
void *run_chunk(void *);
void EM_iter(params *, double **, uint64_t, uint64_t, out_data *, int);

// shared.cpp
void error(const char *, const char *);
void handler(int);
void catch_SIG();
double check_interv(double, bool);
double addProtect2(double, double);
double addProtect3(double, double, double);
double HWE_like(double *, double, double);
double site_HWE_like(double **, double *, double *, uint16_t, uint16_t, uint64_t, uint64_t);
double full_HWE_like(params *, double *, double *, uint16_t, uint16_t);
