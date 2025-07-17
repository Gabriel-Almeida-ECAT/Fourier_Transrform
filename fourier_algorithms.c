#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif


#define TRUE 			1
#define FALSE 			0
#define RESULTS_SHOW 	4

#define q				6		/* for 2^16 samples */
#define N_SAMPLES 	(1<<q)
#define SIGNAL_3

/*
NS => Num Samples
TW => Time Window
FS => Sampling Frequency
*/

typedef enum sampling_opt{
	FORCE_NS_IN_TW, // FS is variable
	FORCE_FS_IN_TW, // Number of sample is variable
	FORCE_NS_AND_FS // TW is variable
} sampling_opt_t;


typedef struct{
	uintmax_t num_samples;
	float init_time;
	float end_time;
	float fs;
	double Ts;
	double* val;
	double* ft_freq;
} signal_t;


void generate_signal(signal_t* signal, 
					uint16_t noise, 
					double sample_freq, 
					float init_time, 
					float end_time, 
					sampling_opt_t samp_opt);

//fazer uma func para salvar uma tbl generia como CSV no futuro
//double* write_tbl_csv(signal_t* signal, char* file_name);
void saveSig2file(signal_t* signal, char* file_name, int is_rslt_sig);
void read_sig_file(signal_t* signal, char* file_name);

void init_sig(signal_t* signal);
void deinit_sig(signal_t* signal);

void dft(signal_t* signal, signal_t* dft_rslt);
void fft(signal_t* signal, signal_t* fft_rslt);
static void recursive_fft(double complex *vals, uintmax_t N, double complex *temp);
double time_ft(void (*ft_func)(signal_t* signal, signal_t* fft_rslt), signal_t* signal, signal_t* fft_rslt);

void plot_results(double tim_dft, double tim_fft);


int main(){
	printf("\n\n=================================================");
	printf("\n# Init Fourier Analizys");
	printf("\n# N_SAMPLES %ld", (uintmax_t)N_SAMPLES);

	signal_t sig, dft_rslt, fft_rslt;
	init_sig(&sig);
	init_sig(&dft_rslt);
	init_sig(&fft_rslt);


	/*=========================================
		SIGNAL
	===========================================*/
#if defined(SIGNAL_1)
	generate_signal(&sig, 
					0, 			//noise
					4, 			//sample freq
					0, 			//init time
					3.0/4.0,
					FORCE_NS_IN_TW);	//end time
#elif defined(SIGNAL_2)
	generate_signal(&sig,
					0, 			//noise
					50, 		//sample freq
					0, 			//init time
					10.0,
					FORCE_NS_IN_TW);	//end time
#elif defined(SIGNAL_3)
	generate_signal(&sig,
					0, 			//noise
					100e6, 		//sample freq
					0, 			//init time
					1e-6,
					FORCE_NS_IN_TW);	//end time
#endif

	printf("\n# Generated signal with %ld samples", sig.num_samples);

	saveSig2file(&sig, "signal.csv", FALSE);
	printf("\n# Saved signal to csv file");

	read_sig_file(&sig, "signal.csv");
	printf("\n# Signal read from csv file");

	/*=========================================
		DFT
	===========================================*/
	printf("\n\n# Calculating DFT");
	//dft(&sig, &dft_rslt);
	double dft_tim = time_ft(dft, &sig, &dft_rslt);
	printf("\n# DFT result: ");
	for(int k=0; (k<sig.num_samples && k<RESULTS_SHOW); k++) printf("%.3lf, ", dft_rslt.val[k]);
	printf("\n# DFT freqs: ");
	for(int k=0; (k<sig.num_samples && k<RESULTS_SHOW); k++) printf("%.3lf, ", dft_rslt.ft_freq[k]);

	saveSig2file(&dft_rslt, "dft_result.csv", TRUE);
	printf("\n# Saved DFT result to csv file");

	/*=========================================
		FFT
	===========================================*/
	printf("\n\n# Calculating FFT");
	//fft(&sig, &fft_rslt);
	double fft_tim = time_ft(fft, &sig, &fft_rslt);
	printf("\n# FFT result: ");
	for(int k=0; (k<sig.num_samples && k<RESULTS_SHOW); k++) printf("%.3lf, ", fft_rslt.val[k]);
	printf("\n# FFT freqs: ");
	for(int k=0; (k<sig.num_samples && k<RESULTS_SHOW); k++) printf("%.3lf, ", fft_rslt.ft_freq[k]);

	saveSig2file(&fft_rslt, "fft_result.csv", TRUE);
	printf("\n# Saved DFT result to csv file");


	deinit_sig(&sig);
	deinit_sig(&dft_rslt);
	deinit_sig(&fft_rslt);

	/*=========================================
		PLOT RESUTLS
	===========================================*/
	printf("\n\n# Plotting graphs ");
	plot_results(dft_tim, fft_tim);


	printf("\n# END Fourier Analizys");
	printf("\n=================================================\n\n");
	return 0;
}


void init_sig(signal_t* signal){
	signal->val = calloc((uintmax_t)N_SAMPLES, sizeof(double));
	signal->ft_freq = calloc((uintmax_t)N_SAMPLES, sizeof(double));
}


void deinit_sig(signal_t* signal){
	free(signal->val);
	free(signal->ft_freq);
}


void generate_signal(signal_t* signal, 
					uint16_t noise, 
					double sample_freq, 
					float init_time, 
					float end_time, 
					sampling_opt_t samp_opt){
	
	double time_window;
	double bin_width;

	signal->init_time = init_time;
	signal->end_time = end_time;

	switch(samp_opt){
		case FORCE_NS_IN_TW: // FS is variable
			signal->num_samples = N_SAMPLES;

			time_window = end_time - init_time;
			signal->Ts = time_window/(double)N_SAMPLES;
			signal->fs = 1.0/signal->Ts;
		break;

		case FORCE_FS_IN_TW: // NS is variable
			signal->fs = sample_freq;

			time_window = end_time - init_time;
			signal->Ts = 1.0/signal->fs;
			signal->num_samples = (uintmax_t)(time_window/signal->Ts) + 1;
		break;

		case FORCE_NS_AND_FS: // TW is variable
			signal->fs = sample_freq;
			signal->num_samples = N_SAMPLES;

			signal->Ts = 1.0/signal->fs;
			time_window = (signal->num_samples*signal->Ts);
		break;
	}
	
	bin_width = (double)signal->fs/(double)signal->num_samples;

	printf("\n# generate_signal():\n\tsamp_opt: %d\n\ttime_window: %.3e,\n\tcalc samples: %ld,\n\tsample_freq: %.3e,\n\tTs: %.3e,\n\tbin_width: %.3e", 
		samp_opt, time_window, signal->num_samples, signal->fs, signal->Ts, bin_width);
	
	double t;
	for(int k=0; (k<signal->num_samples && k<(uintmax_t)N_SAMPLES); k++){
		t = init_time + (k*signal->Ts);

#if defined(SIGNAL_1)
		signal->val[k] = 5 + 2*cos((2*M_PI)*t - (M_PI/2)) + 3*cos(4*M_PI*t);
#elif defined(SIGNAL_2)
		signal->val[k] = sin(2*M_PI*15*t) + sin(2*M_PI*20*t);
#elif defined(SIGNAL_3)
		signal->val[k] = 10*exp(- pow(t/(1/(2*5e6)), 2) )*sin(2*M_PI*5e6*t);
#endif

		if(noise) signal->val[k] += ((double)rand()/RAND_MAX)*noise;
		signal->ft_freq[k] = (double)k*bin_width;
	}

	return;
}


void saveSig2file(signal_t* signal, char* file_name, int is_rslt_sig){
	FILE *file_ptr = fopen(file_name, "w");

	if(file_ptr == NULL){ 
		printf("\n\nError creating file!\n\n"); 
		exit(1); 
	}

	fprintf(file_ptr, "sample,val,x_axis\n");

	for(int k=0; (k<signal->num_samples && k<(uintmax_t)N_SAMPLES); k++){
		fprintf(file_ptr, 
			"%ld,%lf,%e\n",
			k,
			signal->val[k],
			is_rslt_sig ? signal->ft_freq[k] : (double)(signal->init_time + k * signal->Ts)
		);
	}

	fclose(file_ptr);
}


void read_sig_file(signal_t* signal, char* file_name){
	FILE *file_ptr = fopen(file_name, "r");

	if(file_ptr == NULL){ 
		printf("\n\n<!error!> Error openning file!\n\n"); 
		exit(1); 
	}

	int first_line_size = 12;
	char header_line[first_line_size];
	fgets(header_line, first_line_size, file_ptr); // to make the file pointer skip the headers


	printf("\n# read_sig_file(): Read values: ");
	int sample;
	double val, x_axis;
	for(int k=0; (k<signal->num_samples && k<(uintmax_t)N_SAMPLES); k++){
		if(fscanf(file_ptr, "%ld,%lf,%e\n",
			&sample, &val, &x_axis) == 3){
			signal->val[sample] = val;
			if(k <= RESULTS_SHOW) printf("%.3lf, ", signal->val[sample]);
		}
	}

	fclose(file_ptr);
	return NULL;
}


void dft(signal_t* signal, signal_t* dft_rslt){
	uintmax_t N = signal->num_samples;
	dft_rslt->num_samples = N;

	for(int k = 0; k < N; k++) {
	    dft_rslt->ft_freq[k] = signal->ft_freq[k];
	}

	double complex acm;
	double f_k;

	for(int n=0; n < N; n++){
		
		acm = 0.0 + 0.0*I;
		for(int k=0; (k < (N-1) && k<(uintmax_t)N_SAMPLES); k++){
			f_k = signal->val[k];
			acm += f_k*cexp(-I * 2 * M_PI * n * k / (double)N);
		}
 
		dft_rslt->val[n] = cabs(acm);
	}


	return NULL;
}


/*
=> Cooley-Tukey fft algorithm
compute the fourier transform using the cooley-tukey akgorithm
The number of samples *MUST* be 2^n, where n is a natural number
*/
void fft(signal_t* signal, signal_t* fft_rslt){
	uintmax_t N = signal->num_samples;

	if(N == 1){
		fft_rslt->val[0] = signal->val[0];
		return;
	}

	fft_rslt->num_samples = N;
	for(int k = 0; k < N; k++) {
	    fft_rslt->ft_freq[k] = signal->ft_freq[k];
	}

	
	// dont want create suport temp var outside of fft function
	// so i`ll use a suport recursive func
	double complex c_vals[N], temp_val[N];
	for(int n=0; n<N; n++) c_vals[n] = signal->val[n];

	recursive_fft(c_vals, N, temp_val);

	for(int n=0; n<N; n++) fft_rslt->val[n] = cabs(c_vals[n]);

	return;
}


static void recursive_fft(double complex *vals, uintmax_t N, double complex *temp){
	if(N <= 1) return;

	double complex *even, *odd;
	double complex w, z;
	even = temp; odd = temp + (N/2);

	// separate even and odd terms
	for(int k=0; k<(N/2); k++){
		even[k] = vals[2*k];
		odd[k] = vals[(2*k)+1];
	}

	recursive_fft(even, N/2, vals);
	recursive_fft(odd, N/2, vals);

	for(int k=0; k<N/2; k++){
		w = cexp(-I * 2 * M_PI * k / (double)N);
		z = ( creal(w)*creal(odd[k]) - cimag(w)*cimag(odd[k]) ) 
				+ I*( creal(w)*cimag(odd[k]) + cimag(w)*creal(odd[k]) );
		
		vals[k]			= even[k] + z;
		vals[k+(N/2)]	= even[k] - z;
	}

	return;
}

// return the time it takes in milliseconds to exeucte the DFT or FFT
double time_ft(void (*ft_func)(signal_t* signal, signal_t* ft_rslt), 
					  signal_t* signal, signal_t* ft_rslt){
#ifdef _WIN32
    LARGE_INTEGER frequency;
    LARGE_INTEGER start, end;
    
    QueryPerformanceFrequency(&frequency);

    QueryPerformanceCounter(&start);

    ft_func(signal, ft_rslt);

    QueryPerformanceCounter(&end);
    
    double ft_time = (double)(end.QuadPart - start.QuadPart) * 1000000.0 / frequency.QuadPart;
#else 
	struct timespec start, end;

	clock_gettime(CLOCK_MONOTONIC, &start);
    
    ft_func(signal, ft_rslt);
    
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    double ft_time = ((end.tv_sec - start.tv_sec) * 1000.0) + ((end.tv_nsec - start.tv_nsec) / 1e6);
#endif
    
    printf("\n# Calculated ft in %.6lf ms", ft_time);

    return ft_time;
}


void plot_results(double tim_dft, double tim_fft){
	char test_cmd[256];
	char script_call_cmd[256];
	    
#ifdef _WIN32
    strcpy(test_cmd, "python --version >nul 2>&1");
	snprintf(script_call_cmd, sizeof(script_call_cmd), "python get_graphs.py %.6lf %.6lf", tim_dft, tim_fft);
#else
    strcpy(test_cmd, "python3 --version >/dev/null 2>&1");
	snprintf(script_call_cmd, sizeof(script_call_cmd), "python3 get_graphs.py %.6lf %.6lf", tim_dft, tim_fft);
#endif

	if (system(test_cmd) != 0) {
        printf("# Python not available on machine -> exting the program.");
        exit(1);
	}

    if (system(script_call_cmd) == 0) {
        printf("=> Python graph script executed successfully\n");
    } else {
        printf("=> Error executing Python script\n");
    }

    return;
}