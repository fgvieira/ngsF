#include <gsl/gsl_rng.h>
#include "shared.h"



int init_output(params *pars, out_data *output) {
  if( pars->verbose >= 1 )
    printf("==> Setting initial values\n");

        gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set (r, pars->seed);

        double indF_rng_max = 0.99;
        double indF_rng_min = 0.01;
	for(uint16_t i = 0; i < pars->n_ind; i++) {
		output->indF[i] = ( strcmp("r", pars->init_values) == 0 ? indF_rng_min + gsl_rng_uniform(r) * (indF_rng_max - indF_rng_min) : 0.01 );
		output->indF_num[i] = 0;
		output->indF_den[i] = 0;
		output->ind_lkl[i] = 0;
	}

	double freq_rng_max = 0.49;
	double freq_rng_min = 0.01;
	for(uint64_t s = 0; s < pars->n_sites; s++) {
		output->site_freq[s] = ( strcmp("r", pars->init_values) == 0 ? freq_rng_min + gsl_rng_uniform(r) * (freq_rng_max - freq_rng_min) : 0.01 );
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

	  if( pars->verbose >= 1 )
	    printf("\tReading initial values from file: %s\n", pars->init_values);

	  // Check pars file size
	  struct stat st;
	  stat(pars->init_values, &st);
	  if( (uint64_t) st.st_size != sizeof(double) * (1+pars->n_ind+pars->n_sites) )
	    error(__FUNCTION__, "initial parameters file corrupted!");

		gzFile init_values_fh = gzopen(pars->init_values, "rb");
		if( init_values_fh == NULL )
		  error(__FUNCTION__, "cannot open initial parameters file!");

		// Skip Lkl...
		if( gzread (init_values_fh, output->indF, sizeof(double)) < 0 )
		  error(__FUNCTION__, "cannot read Lkl from file!");

		// Read ind F...
		if( gzread (init_values_fh, output->indF, sizeof(double) * pars->n_ind) < 0 )
		  error(__FUNCTION__, "cannot read initial indF from file!");
		for(uint16_t i = 0; i < pars->n_ind; i++) {
		  if(output->indF[i] > indF_rng_max)
		    output->indF[i] = indF_rng_max;
		  if(output->indF[i] < indF_rng_min)
		    output->indF[i] = indF_rng_min;
		}

		// Read site freqs...
		if( gzread (init_values_fh, output->site_freq, sizeof(double) * pars->n_sites) < 0 )
		  error(__FUNCTION__, "cannot read initial freqs from file!");
		for(uint64_t s = 0; s < pars->n_sites; s++) {
		  if(output->site_freq[s] > freq_rng_max)
		    output->site_freq[s] = freq_rng_max;
		  if(output->site_freq[s] < freq_rng_min)
		    output->site_freq[s] = freq_rng_min;
		}

		gzclose(init_values_fh);
	}

	return 0;
}



uint64_t read_chunk(double **chunk_data, params *pars, uint64_t chunk) {
	uint64_t total_elems_read = 0;

	if(chunk >= pars->n_chunks)
		error(__FUNCTION__,"invalid chunk number!");

	// Define chunk start and end positions
	uint64_t start_pos = chunk * pars->max_chunk_size;
	uint64_t end_pos = start_pos + pars->max_chunk_size - 1;
	if(end_pos >= pars->n_sites)	end_pos = pars->n_sites - 1;
	uint64_t chunk_size = end_pos - start_pos + 1;
	if( pars->verbose >= 6 )
	  printf("\tReading chunk %lu from position %lu to %lu (%lu)\n", chunk+1, start_pos, end_pos, chunk_size);

	// Search start position
#ifdef _USE_BGZF
	if( bgzf_seek(pars->in_glf_fh, pars->chunks_voffset[chunk], SEEK_SET) < 0 )
		error(__FUNCTION__,"cannot seek GLF file (BGZF)!");
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
			error(__FUNCTION__,"cannot read GLF file!");
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
