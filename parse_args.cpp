#include <getopt.h>
#include "shared.h"


void init_pars(params *pars) {
	pars->in_glf = NULL;
	pars->in_glf_type = NULL;
	pars->init_values = NULL;
	pars->out_file = NULL;
	pars->n_ind = 0;
	pars->n_sites = 0;
	pars->n_chunks = 0;
	pars->chunks_voffset = NULL;
	pars->max_chunk_size = 10000;
	pars->fast_lkl = false;
	pars->approx_EM = false;
	pars->call_geno = false;
	pars->max_iters = 1500;
	pars->min_epsilon = 1e-5;
	pars->n_threads = 1;
	pars->quick = false;
	pars->version = false;
	pars->verbose = 1;
}



// Parses command line args and stores them into struct params
int parse_cmd_args(int argc, char **argv, params *pars) {

	static struct option long_options[] =
	{
			{"glf", required_argument, NULL, 'g'},
			{"init_values", required_argument, NULL, 'x'},
			{"out", required_argument, NULL, 'o'},
			{"n_ind", required_argument, NULL, 'i'},
			{"n_sites", required_argument, NULL, 's'},
			{"chunk_size", required_argument, NULL, 'c'},
			{"fast_lkl", no_argument, NULL, 'l'},
			{"approx_EM", no_argument, NULL, 'H'},
			{"call_geno", no_argument, NULL, 'G'},
			{"max_iters", required_argument, NULL, 'n'},
			{"min_epsilon", required_argument, NULL, 'e'},
			{"n_threads", required_argument, NULL, 'p'},
			{"quick", no_argument, NULL, 'q'},
			{"version", no_argument, NULL, 'v'},
			{"verbose", required_argument, NULL, 'd'},
			{0, 0, 0, 0}
	};

	int c = 0;
	while ( (c = getopt_long_only(argc, argv, "g:x:o:i:s:c:lHGn:e:p:qvd:", long_options, NULL)) != -1 )
		switch (c) {
		case 'g':
			pars->in_glf = optarg;
			break;
		case 'x':
			pars->init_values = optarg;
			break;
		case 'o':
			pars->out_file = optarg;
			break;
		case 'i':
			pars->n_ind = atoi(optarg);
			break;
		case 's':
			pars->n_sites = atoi(optarg);
			break;
		case 'c':
			pars->max_chunk_size = atoi(optarg);
			break;
		case 'l':
			pars->fast_lkl = true;
			break;
		case 'H':
			pars->approx_EM = true;
			break;
		case 'G':
			pars->call_geno = true;
			break;
		case 'n':
			pars->max_iters = atoi(optarg);
			break;
		case 'e':
			pars->min_epsilon = atof(optarg);
			break;
		case 'p':
			pars->n_threads = atoi(optarg);
			break;
		case 'q':
			pars->quick = true;
			break;
		case 'v':
			pars->version = true;
			break;
		case 'd':
			pars->verbose = atoi(optarg);
			break;
		default:
			exit(-1);
		}

	// Default value for initial values
	if(pars->init_values == NULL) {
		pars->init_values = new char[2];
		pars->init_values[0] = 'e';
		pars->init_values[1] = '\0';
	}

	return 0;
}
