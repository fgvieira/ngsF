#include <gsl/gsl_rng.h>
#include "shared.h"



int init_output(params *pars, out_data *output) {
	double gsl_rng_max, gsl_rng_min;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set (r, time(NULL));

    gsl_rng_max = 0.99;
    gsl_rng_min = 0.01;
	for(uint16_t i = 0; i < pars->n_ind; i++) {
		output->indF[i] = ( strcmp("r", pars->init_values) == 0 ? gsl_rng_min + gsl_rng_uniform(r) * (gsl_rng_max - gsl_rng_min) : 0.01 );
		output->indF_num[i] = 0;
		output->indF_den[i] = 0;
		output->ind_lkl[i] = 0;
	}

	gsl_rng_max = 0.49;
	gsl_rng_min = 0.01;
	for(uint64_t s = 0; s < pars->n_sites; s++) {
		output->site_freq[s] = ( strcmp("r", pars->init_values) == 0 ? gsl_rng_min + gsl_rng_uniform(r) * (gsl_rng_max - gsl_rng_min) : 0.01 );
		output->site_freq_num[s] = 0;
		output->site_freq_den[s] = 0;
		output->site_prob_var[s] = 1;
	}

	output->global_lkl = 0;
	gsl_rng_free (r);

	// If initial values provided...
	if( strcmp("e", pars->init_values) != 0 &&
		strcmp("r", pars->init_values) != 0 &&
		strcmp("u", pars->init_values) != 0) {

		gzFile init_values_fh = gzopen(pars->init_values, "rb");
		if( init_values_fh == NULL )
			error("cannot open init_f GLF file!");

		// Read ind F...
		if( gzread (init_values_fh, output->indF, sizeof(double) * pars->n_ind) < 0 )
			error("cannot read initial indF from file!");

		// Read site freqs...
		if( gzread (init_values_fh, output->site_freq, sizeof(double) * pars->n_sites) < 0 )
			error("cannot read initial freqs from file!");

		gzclose(init_values_fh);
	}

	return 0;
}



uint64_t read_chunk(double **chunk_data, params *pars, uint64_t chunk) {
	uint64_t total_elems_read = 0;

	if(chunk >= pars->n_chunks)
		error("invalid chunk number!");

	// Define chunk start and end positions
	uint64_t start_pos = chunk * pars->max_chunk_size;
	uint64_t end_pos = start_pos + pars->max_chunk_size - 1;
	if(end_pos >= pars->n_sites)	end_pos = pars->n_sites - 1;
	uint64_t chunk_size = end_pos - start_pos + 1;
	if( pars->verbose >= 6 ) printf("\tReading chunk %ld from position %ld to %ld (%ld)\n", chunk+1, start_pos, end_pos, chunk_size);

	// Search start position
#ifdef _USE_BGZF
	if( bgzf_seek(pars->in_glf_fh, pars->chunks_voffset[chunk], SEEK_SET) < 0 )
		error("cannot seek GLF file (BGZF)!");
#endif

	// Read data from file
	for(uint64_t c = 0; c < chunk_size; c++) {
#ifdef _USE_BGZF
		int bytes_read = bgzf_read(pars->in_glf_fh, chunk_data[c], (int) pars->n_ind * 3 * sizeof(double));
		if(pars->call_geno)
			call_geno(chunk_data[c], pars->n_ind, 3);
		uint64_t elems_read = (uint64_t) bytes_read / sizeof(double);
#else
		chunk_data[c] = pars->data[start_pos+c];
		uint64_t elems_read = pars->n_ind * 3;
#endif
		if( elems_read != pars->n_ind * 3 )
			error("cannot read GLF file!");
		total_elems_read += elems_read;
	}

#ifdef _USE_BGZF
	// Update index for next chunk
	if( chunk+1 != pars->n_chunks && pars->chunks_voffset[chunk+1] == 0 )
		pars->chunks_voffset[chunk+1] = bgzf_tell(pars->in_glf_fh);
#endif

	return( total_elems_read/(pars->n_ind * 3) );
}



void call_geno(double *site_gl, int n_ind, int n_geno) {

	for (int i = 0; i < n_ind; i++) {
		int max_pos = array_max_pos(&site_gl[i*n_geno], n_geno);

        for (int j=0; j < n_geno; j++) site_gl[i*n_geno+j] = -1e4;
        site_gl[i*n_geno+max_pos] = 0;
	}
}


int array_max_pos(double *geno, int size) {
	int max = -1e6;
	int res = 0;

	for (int i = 0; i < size; i++) {
		if (geno[i] > max) {
			res = i;
			max = geno[i];
		}
	}
	return res;
}
