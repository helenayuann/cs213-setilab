#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <pthread.h>
#include <unistd.h>
#include <sched.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0

// make a struct for each thread
typedef struct {
  signal* sig;
  double* band_power;
  int filter_order;
  int num_bands;
  int id;
  int num_threads;
  int num_processors;
} thread_data;

// make a worker function that given a thread, performs
// the band power calculations and adds it to the band
// power array
//
// pass in a long i as the parameter that will just be
// where we are in the thread array. i * num_bands per thread
// will tell us our start and end....right....?
//
// if there is left over apparently im supposed to
// divide the work evenly between the threads... idk
// how to do that
//
// maybe just do the suboptimal solution first :D
//
// uhhh bands.
//
// divide the num_samples into the bands and the bands
// into the threads i think
void* worker(thread_data* arg) {
  thread_data* curr_thread = (thread_data*)arg;
  int bands_per_thread = curr_thread->num_bands / curr_thread->num_threads;
  int remainder = curr_thread->num_bands % curr_thread->num_threads;

  // put ourselves on the desired processor
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(curr_thread->id % curr_thread->num_processors, &set);
  if (sched_setaffinity(0, sizeof(set), &set) < 0) { // do it
    perror("Can't setaffinity"); // hopefully doesn't fail
    exit(-1);
  }

  // This figures out the chunk of the vector I should
  // work on based on my id
  int mystart = curr_thread->id * bands_per_thread;
  int myend   = 0;
  double filter_coeffs[curr_thread->filter_order + 1];
  if (curr_thread->id < remainder) {
    mystart += curr_thread->id;
    myend += mystart + bands_per_thread + 1;
  } else {
    mystart += remainder;
    myend = mystart + bands_per_thread;
  }

  double Fc = (curr_thread->sig->Fs) / 2;
  double bandwidth = Fc / curr_thread->num_bands;

  // analyze signal
  for (int band = mystart; band < myend; band++) {
    // Make the filter
    generate_band_pass(curr_thread->sig->Fs,
                       band * bandwidth + 0.0001, // keep within limits
                       (band + 1) * bandwidth - 0.0001,
                       curr_thread->filter_order,
                       filter_coeffs);
    hamming_window(curr_thread->filter_order,filter_coeffs);

    // Convolve
    convolve_and_compute_power(curr_thread->sig->num_samples,
                               curr_thread->sig->data,
                               curr_thread->filter_order,
                               filter_coeffs,
                               &(curr_thread->band_power[band]));

  }

  pthread_exit(NULL);
}

void usage() {
  printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands num_threads num_processors\n");
}

double avg_power(double* data, int num) {

  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double* data, int num) {

  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double* data, int num) {

  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double* data, int num) {

  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (int i = 0; i < num; i++) {
    data[i] -= dc;
  }
}


int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub, int num_threads, int num_processors) {
  pthread_t* tid;
  tid = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
  thread_data* threads;
  threads = (thread_data*)malloc(num_threads * sizeof(thread_data));
  double* band_power;
  band_power = (double*)malloc(num_bands * sizeof(double));

  double Fc = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;

  remove_dc(sig->data,sig->num_samples);

  double signal_power = avg_power(sig->data,sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  resources rstart;
  get_resources(&rstart,THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();

  // all this is in the worker function now
  //
  // double filter_coeffs[filter_order + 1];
  // double band_power[num_bands];
  // for (int band = 0; band < num_bands; band++) {
  //   // Make the filter
  //   generate_band_pass(sig->Fs,
  //                      band * bandwidth + 0.0001, // keep within limits
  //                      (band + 1) * bandwidth - 0.0001,
  //                      filter_order,
  //                      filter_coeffs);
  //   hamming_window(filter_order,filter_coeffs);

  //   // Convolve
  //   convolve_and_compute_power(sig->num_samples,
  //                              sig->data,
  //                              filter_order,
  //                              filter_coeffs,
  //                              &(band_power[band]));

  // }

  // need to initialize all the structs and put them in an array
  for (int i = 0; i <  num_threads; i++) {
    threads[i].sig            = sig;
    threads[i].band_power     = band_power;
    threads[i].filter_order   = filter_order;
    threads[i].id             = i;
    threads[i].num_bands      = num_bands;
    threads[i].num_processors = num_processors;
    threads[i].num_threads    = num_threads;
  }

  // launch the threads!!!!!
  for (long i = 0; i < num_threads; i++) {
    int returncode = pthread_create(&(tid[i]),  // thread id gets put here
                                    NULL, // use default attributes
                                    worker, // thread will begin in this function
                                    (void*)&threads[i] // we'll give it threads[i] as the argument
                                    );
    if (returncode != 0) {
      perror("Failed to start thread");
      exit(-1);
    }
  }

  // now we will join all the threads
  for (int i = 0; i < num_threads; i++) {
    int returncode = pthread_join(tid[i], NULL);
    if (returncode != 0) {
      perror("join failed");
      exit(-1);
    }
  }

  // ok and then analyze...

  unsigned long long tend = get_cycle_count();
  double end = get_seconds();

  resources rend;
  get_resources(&rend,THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power,num_bands);
  double avg_band_power = avg_of(band_power,num_bands);
  int wow = 0;
  *lb = -1;
  *ub = -1;

  for (int band = 0; band < num_bands; band++) {
    double band_low  = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ",
           band, band_low, band_high, band_power[band]);

    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

      // band of interest
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0) {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
      } else {
        printf("(meh)");
      }
    } else {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n\
User time        %lf seconds\n\
System time      %lf seconds\n\
Page faults      %ld\n\
Page swaps       %ld\n\
Blocks of I/O    %ld\n\
Signals caught   %ld\n\
Context switches %ld\n",
         rdiff.usertime,
         rdiff.systime,
         rdiff.pagefaults,
         rdiff.pageswaps,
         rdiff.ioblocks,
         rdiff.sigs,
         rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
         "Note that cycle count only makes sense if the thread stayed on one core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  return wow;
}

int main(int argc, char* argv[]) {

  if (argc != 8) {
    usage();
    return -1;
  }

  char sig_type    = toupper(argv[1][0]);
  char* sig_file   = argv[2];
  double Fs        = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  int num_bands    = atoi(argv[5]);
  int num_threads  = atoi(argv[6]);
  int num_processors = atoi(argv[7]);

  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);

  printf("type:     %s\n\
file:     %s\n\
Fs:       %lf Hz\n\
order:    %d\n\
bands:    %d\n\
threads:  %d\n\
processors: %d\n",
         sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
         sig_file,
         Fs,
         filter_order,
         num_bands,
         num_threads,
         num_processors);

  printf("Load or map file\n");

  signal* sig;
  switch (sig_type) {
    case 'T':
      sig = load_text_format_signal(sig_file);
      break;

    case 'B':
      sig = load_binary_format_signal(sig_file);
      break;

    case 'M':
      sig = map_binary_format_signal(sig_file);
      break;

    default:
      printf("Unknown signal type\n");
      return -1;
  }

  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  double start = 0;
  double end   = 0;
  if (analyze_signal(sig, filter_order, num_bands, &start, &end, num_threads, num_processors)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}

