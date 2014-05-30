#include <signal.h>
#include "shared.h"


bool SIG_COND;


void error(const char *func, const char *msg) {
        printf("\n[%s] ERROR: %s\n", func, msg);
        fprintf(stderr, "\n[%s] ERROR: %s\n", func, msg);
	perror("\t");
	exit(-1);
}


void handler(int s) {
	if(SIG_COND)
		fprintf(stderr,"\n\"%s\" signal caught! Will try to exit nicely (no more threads are created, we will wait for the current threads to finish)\n", strsignal(s));
	SIG_COND = false;
}


//we are threading so we want make a nice signal handler for ctrl+c
void catch_SIG(){
	struct sigaction sa;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = 0;
	sa.sa_handler = handler;

	//sigaction(SIGKILL, &sa, 0);
	sigaction(SIGTERM, &sa, 0);
	sigaction(SIGINT, &sa, 0);
	sigaction(SIGQUIT, &sa, 0);
	sigaction(SIGABRT, &sa, 0);
	sigaction(SIGPIPE, &sa, 0);
}


double check_interv(double value, bool verbose) {
	double errTol = 1e-5;

	if (value != value) {
	        error(__FUNCTION__,"value is NaN!\n");
	} else if(value < errTol) {
		value = 0;
		if(verbose && value < 0) printf("value %f < 0!\n", value);
	} else if(value > 1 - errTol) {
		value = 1;
		if(verbose && value > 1) printf("value %f > 1!\n", value);
	}

	return value;
}


double addProtect2(double a,double b) {
	//function does: log(exp(a)+exp(b)) while protecting for underflow
	double maxVal, sumVal = 0;// = std::max(a,b));

	if(a>b)
		maxVal=a;
	else
		maxVal=b;

	a = (a-maxVal); b = (b-maxVal);
	if  (a > -100) sumVal = sumVal + exp(a);
	if  (b > -100) sumVal = sumVal + exp(b);

	return log(sumVal) + maxVal;
}


double addProtect3(double a, double b, double c) {
	//function does: log(exp(a)+exp(b)+exp(c)) while protecting for underflow
	double maxVal = 0;// = std::max(a,std::max(b,c));
	double sumVal = 0;

	if(a > b && a > c)
		maxVal = a;
	else if(b > c)
		maxVal = b;
	else
		maxVal = c;

	a = (a-maxVal); b = (b-maxVal); c = (c-maxVal);
	if  (a > -100) sumVal = sumVal + exp(a);
	if  (b > -100) sumVal = sumVal + exp(b);
	if  (c > -100) sumVal = sumVal + exp(c);

	return log(sumVal) + maxVal;
}


double HWE_like(double *GL, double p, double F) {
	double p0 = pow(1-p,2) + p*(1-p)*F;
	double p1 = 2*p*(1-p) * (1 - F);
	double p2 = pow(p,2) + p*(1-p)*F;

	return addProtect3( log(p0)+*(GL+0), log(p1)+*(GL+1), log(p2)+*(GL+2) );
}

double site_HWE_like(double **data, double *freq, double *indF, uint16_t ind, uint16_t n_ind, uint64_t site, uint64_t n_sites) {
	double logLike = 0;

	for(uint64_t s = site; s < site+n_sites; s++) {
		for(uint16_t i = ind; i < ind+n_ind; i++)
			logLike += HWE_like( (data[s]+i*3), freq[s], indF[i]);
	}

	return logLike;
}


//TODO: use pthreads for ML calculation
double full_HWE_like(params *pars, double *site_freq, double *indF, uint16_t ind, uint16_t n_ind) {
	double totalLogLike = 0;

	// Initialize data variable
	double **data = new double*[pars->max_chunk_size];
#ifdef _USE_BGZF
	for(uint64_t cnt = 0; cnt < pars->max_chunk_size; cnt++)
		data[cnt] = new double [pars->n_ind*3];
#endif

	for(uint64_t c = 0; c < pars->n_chunks; c++) {
		uint64_t chunk_size = read_chunk(data, pars, c);
		uint64_t ref_s = c * pars->max_chunk_size;
		totalLogLike += site_HWE_like(data, site_freq+ref_s, indF, ind, n_ind, 0, chunk_size);
	}

	// Free data variable
#ifdef _USE_BGZF
	for(uint64_t cnt = 0; cnt < pars->max_chunk_size; cnt++)
		delete [] data[cnt];
#endif
	delete [] data;

	return totalLogLike;
}
