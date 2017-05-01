#include <pthread.h>
#include "shared.h"




int do_EM (params *pars, out_data *output) {
	SIG_COND = true;
	catch_SIG();

	int iter = ( strcmp("e", pars->init_values) == 0 ? 0 : 1);
	double new_global_lkl = 0;
	double lkl_epsilon = 0;
	double est_epsilon = 0;
	sem_init(&pars->launch_thread_semaph, 0, pars->n_threads);
	sem_init(&pars->running_thread_semaph, 0, 0); //To avoid warnings from Valgrind
	sem_init(&pars->running_thread_semaph, 0, -pars->n_chunks);
	pthread_mutex_init(&pars->F_lock, NULL);

	while( (est_epsilon > pars->min_epsilon || lkl_epsilon > pars->min_epsilon || iter <= 10) && iter <= pars->max_iters && SIG_COND ) {

		time_t iter_start = time(NULL);

		// Next Iteration...
		if( pars->verbose >= 1 ) printf("\nIteration %d:\n", iter);



		////////////////////////////////
		// Loop through all chunks... //
		////////////////////////////////
		for(uint64_t c = 0; c < pars->n_chunks; c++) {
			// Wait for room to launch more threads
			while( sem_wait(&pars->launch_thread_semaph) );

			if( pars->verbose >= 5 ) printf("\tChunk %lu of %lu\n", c+1, pars->n_chunks);

			// Declare structure
			pth_params *pth_struct = new pth_params;
			// Reserve memory for chunk data
			pth_struct->chunk_data = new double* [pars->max_chunk_size];
#ifdef _USE_BGZF
			for(uint64_t s = 0; s < pars->max_chunk_size; s++)
				pth_struct->chunk_data[s] = new double[pars->n_ind * 3];
#endif
			// Fill in PThread structure
			pth_struct->pars = pars;
			pth_struct->chunk_size = read_chunk(pth_struct->chunk_data, pth_struct->pars, c);
			pth_struct->chunk_abs_start_pos = c * pars->max_chunk_size;
			pth_struct->iter = iter;
			pth_struct->output = output;

			// Initialize and set thread detached attribute
			pthread_t thread_id;
			pthread_attr_t pt_attr;
			pthread_attr_init(&pt_attr);
			pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_DETACHED);

			// Launch thread
			int rc = pthread_create(&thread_id, &pt_attr, run_chunk, (void*) pth_struct);
			if(rc) error(__FUNCTION__,"pthread_create() failed!");

			if( pars->verbose >= 6 ) {
				int n_free_threads = 0;
				sem_getvalue(&pars->launch_thread_semaph, &n_free_threads);
				printf("Thread launched! Available slots: %d\n", n_free_threads);
			}
			fflush(stdout);
		}



		////////////////////////////////////
		// Wait for all threads to finish //
		////////////////////////////////////
		int n_free_threads = 0;
		do {
			while( sem_wait(&pars->running_thread_semaph) );
			sem_getvalue(&pars->launch_thread_semaph, &n_free_threads);
			if( pars->verbose >= 6 ) printf("Waiting for all threads to finish: %d\n", pars->n_threads - n_free_threads);
		}while(n_free_threads < (int) pars->n_threads);



		est_epsilon = new_global_lkl = 0;
		/////////////////////////////////////
		// Indiv post-iteration processing //
		/////////////////////////////////////
		if( pars->verbose >= 2 ) printf("\tInd F:\t");
		for(uint16_t i = 0; i < pars->n_ind; i++) {
			// Get new indF and check for interval...
			double new_indF = check_interv(output->indF_num[i] / output->indF_den[i], false);
			// Calculate iter epsilon
			est_epsilon += pow(new_indF - output->indF[i], 2);
			// Store new indF
			new_indF = ( new_indF == 1 ? 0.9999 : new_indF );
			output->indF[i] = new_indF;

			// Calculate overall lkl
			output->ind_lkl[i] = full_HWE_like(pars, output->site_freq, output->indF, i, 1);
			new_global_lkl += output->ind_lkl[i];

			// Reset variables...
			output->indF_num[i] = 0;
			output->indF_den[i] = 0;

			// Debug
			if( pars->verbose >= 2 ) printf("\t%.9f", output->indF[i]);
		}
		if( pars->verbose >= 2 ) printf("\n");



		////////////////////////////////////
                // Site post-iteration processing //
                ////////////////////////////////////
                if( pars->verbose >= 4 ) printf("\tFreq:\t");
                for(uint64_t s = 0; s < pars->n_sites; s++) {
                  if(output->site_freq[s] == 0) continue;
                  if(!pars->freq_fixed){
                    double new_site_freq = check_interv(output->site_freq_num[s] / output->site_freq_den[s], true);
                    est_epsilon += pow(new_site_freq - output->site_freq[s], 2);
                    output->site_freq[s] = (new_site_freq > 0.99 ? 0.99 : new_site_freq);
                  }

                  // Reset variables...
                  output->site_freq_num[s] = 0;
                  output->site_freq_den[s] = 0;
                  output->site_prob_var[s] = output->site_tmpprob_var[s];
                  output->site_tmpprob_var[s] = 0;

                  // Debug
                  if( pars->verbose >= 4 ) printf("\t%.9f", output->site_freq[s]);
                }
                if( pars->verbose >= 4 ) printf("\n");



		///////////////////////
		// Calculate epsilon //
		///////////////////////
		// Parameter epsilon
		est_epsilon = sqrt(est_epsilon)/(pars->n_ind + pars->n_sites);
		// Lkl epsilon calculation - On first iteration, since there is no global_lkl, calculate Lkl epsilon based on current lkl
		output->global_lkl = (output->global_lkl == 0 ? new_global_lkl : output->global_lkl);
		lkl_epsilon = (new_global_lkl - output->global_lkl)/fabs(output->global_lkl);
		output->global_lkl = new_global_lkl;


		if( pars->verbose >= 1 ) {
			time_t iter_end = time(NULL);
			printf("\tLogLkl: %.15f\t epsilon: %.15f %.15f\ttime: %.0f (s)\n", output->global_lkl, lkl_epsilon, est_epsilon, difftime(iter_end, iter_start) );
		}
		iter++;
		fflush(stdout);



		///////////////////////////////
		// Dump iteration parameters //
		///////////////////////////////
		char* pars_file = (char*) malloc( (strlen(pars->out_file)+5+1)*sizeof(char) );
		memset(pars_file, '\0', (strlen(pars->out_file)+5+1)*sizeof(char));
		strcat(pars_file, pars->out_file); strcat(pars_file, ".pars");
		// Write the last iteration to disk
		FILE* last_est_pars = fopen(pars_file, "w");
		if(last_est_pars == NULL)
		  error(__FUNCTION__, "Cannot open PARS file!");
		fwrite(&output->global_lkl, sizeof(double), 1, last_est_pars);
		fwrite(output->ind_lkl, sizeof(double), pars->n_ind, last_est_pars);
		fwrite(output->indF, sizeof(double), pars->n_ind, last_est_pars);
		fwrite(output->site_freq, sizeof(double), pars->n_sites, last_est_pars);
		fclose(last_est_pars);
		free(pars_file);


		if( pars->quick ) break;
	}


	if( iter > pars->max_iters )
		printf("WARN: Maximum number of iterations reached! Check if analysis converged... \n");

	return 0;
}



void *run_chunk(void *pth_struct) {
	pth_params *p = (pth_params *) pth_struct;

	// Get freq estimate per site
	EM_iter(p->pars, p->chunk_data, p->chunk_abs_start_pos, p->chunk_size, p->output, p->iter);

	// Signal semaphores
	if( sem_post(&p->pars->launch_thread_semaph) )
		printf("WARN: launch thread semaphore post failed!\n");
	if( sem_post(&p->pars->running_thread_semaph) )
		printf("WARN: running thread semaphore post failed!\n");

	// Debug
	if( p->pars->verbose >= 6 ) {
		int n_free_threads = 0;
		sem_getvalue(&p->pars->launch_thread_semaph, &n_free_threads);
		printf("Thread finished! Still running: %d\n", p->pars->n_threads - n_free_threads);
	}

	// Free pthread structure memory
#ifdef _USE_BGZF
	for(uint64_t s = 0; s < p->pars->max_chunk_size; s++)
		delete [] p->chunk_data[s];
#endif
	delete [] p->chunk_data;
	delete p;

	pthread_exit(NULL);
}




void EM_iter(params *pars, double **chunk_data, uint64_t chunk_abs_start_pos, uint64_t chunk_size, out_data *output, int iter) {

	// Loop over all sites
	for(uint64_t s = 0; s < chunk_size; s++) {
		uint64_t abs_s = s + chunk_abs_start_pos;
		double p = output->site_freq[abs_s];

		// Skip site if freq == 0
		if (p == 0) continue;

		// Loop over all individuals
		for(uint64_t i = 0; i < pars->n_ind; i++) {
			double F = output->indF[i];
			double p0 = pow(1-p,2) + p*(1-p)*F;
			double p1 = 2*p*(1-p) * (1 - F);
			double p2 = pow(p,2) + p*(1-p)*F;

			// If initial guess assumes uniform priors
			if(iter == 0) p0 = p1 = p2 = 1;

			double norm = addProtect3(log(p0)+chunk_data[s][i*3+0], log(p1)+chunk_data[s][i*3+1], log(p2)+chunk_data[s][i*3+2]);
			double pp0 = p0 * exp(chunk_data[s][i*3+0]-norm);
			double pp1 = p1 * exp(chunk_data[s][i*3+1]-norm);
			double pp2 = p2 * exp(chunk_data[s][i*3+2]-norm);

			output->site_tmpprob_var[abs_s] = fmax(output->site_tmpprob_var[abs_s], 1 - pp0);

			double IBD = 0;
			double indF_num = 0;
			double indF_den = 0;
			// P(IBD)
			if(iter == 0) { //if initial guess do not use any prior
				IBD = check_interv(1 - (pp1/(2*(1-p)*p)), false);
				indF_num = IBD;
				indF_den = 1;
			} else {
				if(pars->approx_EM) { // Vieira et al. algorithm
					double a0 = pp0*(1-p)*p;
					double b0 = pow(1-p,2);
					double c0 = (1-p)*p;
					double a1 = pp1*2*(1-p)*p;
					double b1 = 2*(1-p)*p;
					double c1 = -2*(1-p)*p;
					double a2 = pp2*(1-p)*p;
					double b2 = pow(p,2);
					double c2 = (1-p)*p;

					double s_num0 = a0/(b0+F*c0) + (a0*c0*F)/pow(b0+F*c0,2);
					double s_num1 = a1/(b1+F*c1) + (a1*c1*F)/pow(b1+F*c1,2);
					double s_num2 = a2/(b2+F*c2) + (a2*c2*F)/pow(b2+F*c2,2);
					double s_den0 = (a0*c0)/pow(b0+F*c0,2);
					double s_den1 = (a1*c1)/pow(b1+F*c1,2);
					double s_den2 = (a2*c2)/pow(b2+F*c2,2);

					indF_num = s_num0 - s_num1 + s_num2;
					indF_den = s_den0 - s_den1 + s_den2;
					if( indF_num != indF_num || indF_den != indF_den || indF_num/indF_den != indF_num/indF_den || pars->verbose >= 7 )
					  printf("site IBD: %lu %lu %f %f / %f %f %f %f %f %f %f %f %f / %f %f %f %f %f %f %f %f %f / %f %f %f %f %f %f / %f %f\n", 
						 abs_s, i, p, F, 
						 chunk_data[s][i*3+0], chunk_data[s][i*3+1], chunk_data[s][i*3+2], p0, p1, p2, pp0, pp1, pp2, 
						 a0, b0, c0, a1, b1, c1, a2, b2, c2, 
						 s_num0, s_num1, s_num2, s_den0, s_den1, s_den2, 
						 indF_num, indF_den);
					IBD = check_interv(indF_num/indF_den, false);

				} else { // Hall et al. algorithm
					IBD = pp0 * (1-p) * F / ((1-p) * F + pow(1-p,2) * (1-F)) + pp2 * p * F / (p * F + pow(p,2) * (1-F));
					indF_num = IBD;
					indF_den = 1;
				}
			}

			// Update site freq
			output->site_freq_num[abs_s] += pp1 + pp2 * (2-IBD);
			output->site_freq_den[abs_s] += pp1 + pp2 * (2-IBD) + pp1 + pp0*(2-IBD);

			// Update indiv F
			pthread_mutex_lock(&pars->F_lock);
			output->indF_num[i] += indF_num;// * pow(output->site_prob_var[abs_s], 100);
			output->indF_den[i] += indF_den;// * pow(output->site_prob_var[abs_s], 100);
			pthread_mutex_unlock(&pars->F_lock);

			if( pars->verbose >= 7 ) printf("Ind: %lu\t%.10f %.10f %.10f\tfa: %f\tindF: %f\tp: %f %f %f\tpp: %f %f %f\tCum_freq: %f (%f/%f)\tCumF: %f (%f/%f)\n",
					i+1, chunk_data[s][i*3+0], chunk_data[s][i*3+1], chunk_data[s][i*3+2], \
					p, F, p0, p1, p2, pp0, pp1, pp2, \
					output->site_freq_num[abs_s]/output->site_freq_den[abs_s], output->site_freq_num[abs_s], output->site_freq_den[abs_s], \
					output->indF_num[i]/output->indF_den[i], output->indF_num[i], output->indF_den[i]);
		}

		if( pars->verbose >= 6 ) printf("\t\t%lu\t%f (%f / %f) %f\n", abs_s+1, output->site_freq_num[abs_s]/output->site_freq_den[abs_s], output->site_freq_num[abs_s], output->site_freq_den[abs_s], output->site_prob_var[abs_s]);
	}
}
