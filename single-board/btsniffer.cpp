/* 
 * Copyright 2019-2020 Marco Cominelli
 * Copyright 2019-2020 Francesco Gringoli
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/types/tune_request.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <iostream>
#include <cstdlib>
#include <pthread.h>
#include <csignal>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <unordered_map>
#include "lapnode.hpp"

namespace po = boost::program_options;

typedef std::complex<double> iqsamp_t;

/*** Select number of channels ***/
#define CHAN_8
#include "btsniffer.hpp"

const unsigned int decfactor = (int) raw_srate/srate;

unsigned long total_num_samps = 0;

volatile unsigned int bufselect;
pthread_spinlock_t lock[2];

volatile static bool stopsig = false;
void sigint_handler(int) {stopsig = true;}

typedef struct ProcPars {
	iqsamp_t *bankA;
	iqsamp_t *bankB;
	size_t bufsize;
} proc_pars_t;

inline bool is_valid_preamble(uint8_t *binbuf, unsigned int k)
{
	uint8_t preamble1 = 0, preamble2 = 0;
	preamble1 = binbuf[k] + binbuf[k+2*srate];
	preamble2 = binbuf[k+1*srate] + binbuf[k+3*srate];
	
	if ((preamble1 == 2 && preamble2 == 0) or (preamble1 == 0 && preamble2 == 2))
		return true;
	else
		return false;
}

inline uint8_t extract_byte(uint8_t *binbuf, unsigned int start)
{
	uint8_t result = 0x00;
	for (int b = 0; b < 8; b++) result |= binbuf[start + b*srate] << b;
	return result;
}

struct timeval t_start;

void compute_tsdiff(struct timeval *x, struct timeval *y, struct timeval *diff)
{
	// x must be bigger than y
	int x_usec = x->tv_usec;
	int x_sec = x->tv_sec;
	int y_usec = y->tv_usec;
	int y_sec = y->tv_sec;
	
	if (x_usec < y_usec) {
		x_sec --;
		x_usec += 1000000;
	}
	
	diff->tv_sec = x_sec - y_sec;
	diff->tv_usec = x_usec - y_usec;
	
  return;
}


int _length (uint64_t word, int left, int right)
{
  if (left == right)
    return left;
  
  int mid = (left + right) / 2;
  if (right == mid)
    return mid;
  
  if (word >= (1LLU << mid))
    return (_length (word, left, mid));
  
  return (_length (word, mid, right));
}


int length (uint64_t word)
{
  return (_length (word, 64, 0));
}

std::unordered_map<uint32_t, lap_node> lap_map;

uint64_t compute_remainder (uint64_t input, uint64_t divisor)
{
  int divisor_length = length (divisor);
  int input_length = length (input);
  
  if (divisor_length + input_length > 63)
    return input;
  
  input = input << divisor_length;
  
  while (length (input) >= divisor_length) {
    uint64_t tmp = divisor << (length (input) - divisor_length);
    input = input ^ tmp;
  }
  
  return input;
}

/* Extracts the header if FEC is perfect, then it tries (bruteforce)
 * all possible clk values to dewhiten the header and finally
 * tries all possible UAP until a valid HEC code is found.
 */
int extract_header_bf (uint8_t *buf, uint32_t* head, uint32_t *clks, int *clk_found, uint32_t uap)
{
  int perfect_rx = 0;
  uint32_t header = 0;
  
  for (int i = 0; i < 54; i += 3) {
    int s0 = 0, s1 = 0;
    for (int j = 0; j < 3; j++) {
      if ( *(buf+(i+j)*srate) ) s1++;
      else s0++;
    }
    header >>= 1;
    if (s1 == 0 || s0 == 0) perfect_rx++;
    if (s1 > s0) header |= 0x20000;
  }

  *head = header;

  if (perfect_rx != 18) return -1;

  // first bit received is the LSB, so header is composed by
  //     +-----+---------+
  // MSB | HEC | LT_ADDR | LSB
  //     +-----+---------+
  //      8 bit   10 bit

  // Brute force the clock
  *clk_found = 0;
  for (uint32_t clk = 0; clk < 64; clk++) {
    uint32_t header_dewhiten = header;
    uint32_t whitener;
    whitener = (clk & 0x3f) | 0x40;

    // Dewhiten header using this clk value
    for (int i = 0; i < 18; i++) {
      uint32_t whitener_out = (whitener >> 6) & 0x1;
      uint32_t whitener_shifted = (whitener << 1) & 0x7f;
      whitener = whitener_shifted ^ (whitener_out | (whitener_out << 4));
      header_dewhiten = header_dewhiten ^ (whitener_out << i);
    }

    // Re-compute the HEC over the dewhitened header
    uint32_t lfsr = uap;
    for (int i = 0; i < 10; i++) {
      uint32_t lfsr_out = (lfsr >> 7) & 0x1;
      uint32_t data_in = (header_dewhiten >> i) & 0x1;
      uint32_t lfsr_in = (lfsr_out ^ data_in);
      uint32_t lfsr_adder =
	(lfsr_in << 7) |
	(lfsr_in << 5) |
	(lfsr_in << 2) |
	(lfsr_in << 1) |
	(lfsr_in << 0);
      lfsr = (lfsr << 1) & 0xff;
      lfsr = lfsr ^ lfsr_adder;
    }

    // Compare HEC computed and HEC received.
    // First bit received is in header_dewhiten[10], last is in header_dewhiten[17].
    // First bit to be transmitted is in position 7, last one is in position 0.
    int kk = 0;
    while (kk < 8) {
      uint32_t bit_rx = (header_dewhiten >> (10 + kk)) & 0x1;
      uint32_t bit_tx = (lfsr >> (7 - kk)) & 0x1;
      if (bit_rx != bit_tx) break;
      kk++;
    }

    if (kk == 8) {
      clks[*clk_found] = clk;
      (*clk_found)++;
    }
  }

  return *clk_found;
}


void* proc_routine(void *routine_params)
{
  // Set priority on current thread
  uhd::set_thread_priority_safe(1, true);
  
  // Read parameters.
  proc_pars_t *pars = (proc_pars_t *) routine_params;
  size_t bufsize = pars->bufsize;
  iqsamp_t *bankA = pars->bankA;
  iqsamp_t *bankB = pars->bankB;

  iqsamp_t *curbuf;
  size_t samples_processed = 0;
  FILE *fptrout = fopen("results.txt","w");
  int local_bufselect;

  // Allocate buffers and auxiliary pointers.
  uint8_t *binbuffer = (uint8_t*) malloc(decfactor * bufsize * sizeof(uint8_t));
  iqsamp_t *sigbuf = (iqsamp_t*) malloc(decfactor * bufsize * sizeof(iqsamp_t));
  iqsamp_t *chanbuf = (iqsamp_t*) malloc(decfactor * bufsize * sizeof(iqsamp_t));

  // Make filter polyphase!
  std::vector<std::vector<double>> poly(decfactor);
  size_t components_length = (size_t) ceil(FILTER_TAP_NUM/decfactor);
  for (size_t i = 0; i < components_length * decfactor; i++) {
	  if (i < FILTER_TAP_NUM) poly[i%decfactor].push_back(filter_taps[i]);
	  else poly[i%decfactor].push_back(0.0);
  }

  // Create twiddle matrix.
  std::complex<double> i_unit(0.0, 1.0);
  std::vector<std::vector<iqsamp_t>> twiddle(decfactor);
  for (size_t row = 0; row < decfactor; row++) {
	  for (size_t col = 0; col < decfactor; col++) {
		  std::complex<double> tmp =
			  exp(-i_unit * 2.0*M_PI * ((double) col*row) / (double)decfactor);
		  twiddle[row].push_back(tmp);
	  }
  }

  while(stopsig == false) {
    local_bufselect = bufselect;
    
    if (local_bufselect) {
      curbuf = bankB;
      pthread_spin_lock(&lock[1]);
    } else {
      curbuf = bankA;
      pthread_spin_lock(&lock[0]);
    }

    // Filtering
    memset(sigbuf, 0, decfactor * bufsize * sizeof(iqsamp_t));

    for (size_t i = FILTER_TAP_NUM-1; i < bufsize*decfactor; i += decfactor) {
	    for (size_t k = 0; k < components_length; k++) {
		    size_t idx = i/decfactor;
		    sigbuf[idx-1] += (curbuf[i-k*decfactor+0] * poly[0].at(k))/((double) components_length);
		    for (size_t ch = 1; ch < decfactor; ch++) {
			    sigbuf[ch*bufsize + idx] += (curbuf[i-k*decfactor+(decfactor-ch)] * poly[ch].at(k)/((double) components_length-1));
		    }
	    }
    }

    const size_t blocksize = 1000; // us
    const size_t blocknsamps = blocksize*srate; // samples per block
    const size_t nblocks = floor(bufsize/blocknsamps);

    // Channelize using DFT.
    memset(chanbuf, 0, decfactor * bufsize * sizeof(iqsamp_t));
    for (unsigned int i = 0; i < bufsize; i++) {
	    for (size_t ch = 0; ch < decfactor; ch++)
		    for (size_t k = 0; k < decfactor; k++)
			    chanbuf[ch*bufsize + i] += sigbuf[k*bufsize+i] * twiddle[ch].at(k);
    }
    
    for (unsigned int ch=0; ch<decfactor; ch++) {
	    iqsamp_t *chan = chanbuf + ch * bufsize;
	    uint8_t *tmpbinbuf = binbuffer + ch * bufsize;

	    // Discriminate bits without using atan2.
	    for (size_t i = 1; i < bufsize; i++) {
		    double tmp = chan[i-1].real() * chan[i].imag() - chan[i-1].imag() * chan[i].real();
		    tmpbinbuf[i] = (tmp > 0) ? 1 : 0;
	    }
    }

    for (size_t block = 1; block < nblocks-1; block++) {
	    for (unsigned int ch = 0; ch < decfactor; ch++) {
		    iqsamp_t *chan = chanbuf + ch * bufsize;
		    uint8_t *binbuf = binbuffer + ch * bufsize;

		    for (unsigned int i = block*blocknsamps; i < (block+1)*blocknsamps; i++) {
			    if (is_valid_preamble(binbuf, i) == false) continue;

			    uint64_t barker = extract_byte(binbuf, i + 62*srate);
			    barker = barker & 0x3f;
			    if (barker != 0x13 && barker != 0x2c) continue;
	    
			    uint64_t lap =
				    (uint64_t) extract_byte(binbuf, i + 54*srate) << 16 |
				    (uint64_t) extract_byte(binbuf, i + 46*srate) << 8  |
				    (uint64_t) extract_byte(binbuf, i + 38*srate);
			    
			    uint64_t code =
				    ((uint64_t) extract_byte(binbuf, i +  4*srate) <<  0) |
				    ((uint64_t) extract_byte(binbuf, i + 12*srate) <<  8) |
				    ((uint64_t) extract_byte(binbuf, i + 20*srate) << 16) |
				    ((uint64_t) extract_byte(binbuf, i + 28*srate) << 24) |
				    ((uint64_t) extract_byte(binbuf, i + 36*srate) << 32);
			    code = code & 0x3FFFFFFFFLLU;
			    
			    uint64_t aw = ((uint64_t) barker << 58) | (lap << 34) | code;
			    
			    // use lap to rebuild access word from scratch, do not use barker
			    // set barker accordingly to extracted lap.
			    uint64_t barker_true =  ((lap & 0x800000) != 0) ? 0x13 : 0x2c;
			    
			    uint64_t x = (barker_true << 24) | lap;
			    uint64_t p = 0x83848D96BBCC54FC;
			    uint64_t xtilde = (p >> 34) ^ x;
			    uint64_t gp = 0157464165547;
			    uint64_t g = (gp << 1) ^ gp;
			    uint64_t ctilde = compute_remainder (xtilde, g);
			    uint64_t stilde = ctilde | (xtilde << 34);
			    uint64_t awfinal = stilde ^ p;
			    
			    uint32_t _lap = (uint32_t) lap;
			    if (aw == awfinal) {
				    if (lap_map.find(_lap) == lap_map.end()) {
					    lap_map[_lap] = lap_node(_lap);
				    }
			    } else {
				    continue;
			    }
	    
			    if (stopsig) break;
	    
#define INVALID_CLK_INDEX -1
#define DELTA_TS_SAME_THRESHOLD 40 // this should depend on frame length!
#define DELTA_TS_SLOT_THRESHOLD 620 // should be 625, give margin
#define SLOT_DURATION 625.0
#define ERROR_THRESHOLD 0.05
			    
			    uint32_t header = 0;
			    uint32_t clk_table[64];
			    long long timenow_sec_us = (samples_processed + i)/srate;
			    
			    if (stopsig == false) {
				    std::cout << boost::format("[%2d] %12lld us -- %06X -- ")
					    % ch % (timenow_sec_us) % (_lap);
			    }
			    
			    lap_map[_lap].increase_processed_packets();
			    switch (lap_map[_lap].get_status()) {
			    case LAP_STATE_NEW:
				    {
					    int a;
					    uint32_t uap;
					    lap_map[_lap].set_ts(timenow_sec_us);
					    lap_map[_lap].set_tstart(timenow_sec_us);
                        //std::cout << boost::format("started %lld -- ") % timenow_sec_us;
					    int valid_uaps = 0;
					    for (uap = 0; uap < 256; uap ++) {
						    int clk_found = extract_header_bf(binbuf+i+72*srate,
						                                      &header, clk_table, &a, uap);
						    if (clk_found <= 0) {
							    lap_map[_lap].bf_failed(); // don't log but keep trace of LAP
							    continue;
						    }
						    else if (clk_found != 2) {
							    // This should never happen.
							    std::cerr << "Invalid number of clk values " << clk_found << std::endl;
							    stopsig = true;
							    continue;
						    }
						    lap_map[_lap].set_uap_data(valid_uaps, true, uap,
						                               clk_table[0], clk_table[1], INVALID_CLK_INDEX);
						    valid_uaps++;
					    }
					    if (valid_uaps != 32) {
						    std::cout << boost::format("Init failed") << std::endl;
						    lap_map[_lap].bf_cannot_init();
					    } else {
						    std::cout << boost::format("Initialized") << std::endl;
						    lap_map[_lap].set_status(LAP_STATE_BRUTE_FORCING);
					    }
				    }
				    break;
			    case LAP_STATE_BRUTE_FORCING:
				    {
					    long long tmpdeltats = timenow_sec_us - lap_map[_lap].get_ts();
					    if (tmpdeltats < 0) {
						    // Skip packet
						    std::cout << "Skip packet in the past" << std::endl;
						    lap_map[_lap].count_packet_inthepast ();
						    continue;
					    }
					    
					    if (llabs(lap_map[_lap].get_ts() - timenow_sec_us) < DELTA_TS_SAME_THRESHOLD) {
						    // This might happen with packet received on side channels.
						    std::cout << "Skip packet (too close to another one) ";
						    int a, clk_index;
						    int confirmed_uap = 0;
						    int valid_uap = 0;
						    for (int jj = 0; jj < 32; jj ++) {
							    uint32_t clks[2];
							    uint32_t uap;
							    bool uap_valid;
							    lap_map[_lap].get_uap_data(jj, &uap_valid, &uap,
							                               &clks[0], &clks[1], &clk_index);
							    if(uap_valid == false) continue;
							    valid_uap ++;
							    int clk_found = extract_header_bf(binbuf+i+72*srate,
							                                      &header, clk_table, &a, uap);
							    if (clk_found != 2) continue;
							    if (clk_table[0] != clks[0] || clk_table[1] != clks[1]) continue;
							    confirmed_uap ++;
						    }
						    
						    if (confirmed_uap == valid_uap) lap_map[_lap].new_packet_too_close (true);
						    else lap_map[_lap].new_packet_too_close (false);
						    
						    std::cout <<
							    boost::format("Confirmed %d out of %d UAPs for LAP %0X. Skipping")
							    % confirmed_uap % valid_uap % _lap << std::endl;
						    lap_map[_lap].set_ts(timenow_sec_us);
						    continue;
					    } else if (llabs(lap_map[_lap].get_ts() - timenow_sec_us) < DELTA_TS_SLOT_THRESHOLD) {
						    // we cannot handle such situation at the moment, so simply exit
						    std::cerr << "Cannot handle packet type, removing LAP" << std::endl;
						    lap_map.erase(_lap);
						    continue;
					    } else {
						    // measure how far away we are from being a multiple of 625
						    long long prevts = lap_map[_lap].get_ts();
						    long long deltats = timenow_sec_us - prevts;
						    float deltats_float = (float) deltats;
						    float periods = deltats_float / 625.0;
						    float periods_round = roundf(periods);
						    float error = fabsf(periods - periods_round);
						    if (error > ERROR_THRESHOLD) {
							    std::cerr << boost::format("Error too big (%f), LAP removed") % error
							              << std::endl;
							    lap_map.erase(_lap);
							    continue;
						    }
						    
						    int a, clk_index;
						    int count_valid_uap = 0;
						    int count_broken_uap = 0;
						    for (int jj = 0; jj < 32; jj ++) {
							    uint32_t clks[2];
							    uint32_t uap;
							    bool uap_valid;
							    lap_map[_lap].get_uap_data(jj, &uap_valid, &uap,
							                               &clks[0], &clks[1], &clk_index);
							    if(uap_valid == false)
								    continue;
							    int clk_found = extract_header_bf(binbuf+i+72*srate,
							                                      &header, clk_table, &a, uap);
							    if (clk_found != 2) {
								    // if we end up here it means either
								    // 1. the connection changed, we should remove the old one and restart
								    // 2. this frame is corrupt
								    // At the moment we handle case 2 only, so we simply ignore this uap
								    // and we make sure at the end that all were broken
								    count_broken_uap++;
								    continue;
							    }
							    if (clk_table[0] == clks[0] && clk_table[1] == clks[1]) {
								    count_valid_uap++;
								    continue;
							    }
							    int slot = ((int) periods_round) % 64;
							    // test only the valid one
							    int old_to_check = 2;
							    if (clk_index != INVALID_CLK_INDEX) {
								    clks[0] = clks[clk_index];
								    old_to_check = 1;
							    }
							    int qnew, qold;
							    int qnew_valid = INVALID_CLK_INDEX;
							    int matched = 0;
							    for (qold = 0; qold < old_to_check; qold ++) {
								    for (qnew = 0; qnew < 2; qnew ++) {
									    int slot_guessed = (clk_table[qnew] - clks[qold]) % 64;
									    if (slot == slot_guessed) {
										    qnew_valid = qnew;
										    matched ++;
									    }
								    }
							    }
							    if (matched > 1) {
								    lap_map[_lap].set_uap_data(jj, true, uap, clk_table[0], clk_table[1],
								                               INVALID_CLK_INDEX);
								    count_valid_uap ++;
							    } else if (matched > 0) {
								    lap_map[_lap].set_uap_data(jj, true, uap, clk_table[0], clk_table[1],
								                               qnew_valid);
								    count_valid_uap ++;
							    } else {
								    lap_map[_lap].set_uap_data(jj, false, 0, 0, 0, INVALID_CLK_INDEX);
							    }
						    }
						    if (count_valid_uap == 0 && count_broken_uap > 0) {
							    lap_map[_lap].count_broken_uap ();
							    std::cerr << "Frame likely broken, skipped" << std::endl;
							    continue;
						    } else if (count_valid_uap == 0) {
							    std::cout << "No valid UAP remaining, LAP removed" << std::endl;
							    lap_map.erase(_lap);
							    continue;
						    } else if (count_valid_uap <= 2) {
							    uint32_t clks[2];
							    uint32_t uap;
							    bool isvalid;
							    int clock_index;
							    uint32_t uap_found[2];
							    int uap_idx = 0;
							    
 							    struct timeval ts;
 							    gettimeofday(&ts, NULL);
 							    
							    for (int i = 0; i < 32; i++) {
								    lap_map[_lap].get_uap_data(i, &isvalid, &uap,
								                               &clks[0], &clks[1], &clock_index);
								    if (isvalid) {
									    uap_found[uap_idx++] = uap;
								    }
							    }
							    
							    std::cout << boost::format("Only two UAP left (%d and %d) - ")
								    % uap_found[0] % uap_found[1];
							    
							    // Check that only 2 UAPs have been found.
							    if (uap_idx > 2)
								    fprintf(fptrout, "ERROR\n");
							    
							    long long tsfirst = lap_map[_lap].get_tstart();
							    long long tsresolv = (samples_processed + i)/srate;
							    long long tsdiff = tsresolv - tsfirst;

                                /* Compute energy of last pkt */
                                double energy = 0;
                                for (size_t k = i; k < i + 64*srate; k++) {
                                    energy += (chan[k].real() * chan[k].real()) + (chan[k].imag() * chan[k].imag());
                                }

                                std::cout << boost::format("first: %lld -- ") % tsfirst;
							    
                                std::cout << boost::format("solved in %lld us")  % tsdiff << std::endl;

                                fprintf(fptrout,
                                    "%ld.%06ld %06X -- %3d %3d -- ts %lld -- tdiff %lld -- energy %lf\n",
 								    ts.tv_sec, ts.tv_usec, _lap, uap_found[0], uap_found[1], tsresolv,
                                    tsdiff, energy);
 							    fflush (fptrout);

                                // LAP resolved, remove LAP to restart.
                                lap_map.erase(_lap);
							    
						    } else {
							    std::cout << boost::format("%d possible UAPs remaining") % count_valid_uap
							              << std::endl;
						    }
						    
						    lap_map[_lap].set_ts(timenow_sec_us);
					    }
				    }
				    break;
			    default:
				    {
				    }
				    break;
			    }
			    
			    i += 100;
		    }
	    }
		    
    }

    // Increment sample count.
    samples_processed += bufsize;
    
    
    // Release buffer.
    if (local_bufselect) pthread_spin_unlock(&lock[1]);
    else pthread_spin_unlock(&lock[0]); 
  }
  
  free(chanbuf);
  free(sigbuf);
  free(binbuffer); 
  fclose(fptrout);
  
  return NULL;
}


// Entry point.
int UHD_SAFE_MAIN(int argc, char *argv[])
{
  // Set highest priority on current thread.
  uhd::set_thread_priority_safe(1, true);

  std::cout << "Realtime Bluetooth Sniffer";
  std::cout << std::endl << std::endl;

  std::string args, subdev, antname, channel_list;
  double rate, freq, gain;
  
  // Initialise program options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
      "print this help message")
    ("args", po::value<std::string>(&args)->default_value(""),
     "USRP device arguments")
    ("channels", po::value<std::string>(&channel_list)->default_value("0"),
     "select channel to use (e.g. \"0\", \"0,1\")")
    ("subdev", po::value<std::string>(&subdev)->default_value("A:A"),
     "set frontend specification (e.g. \"A:A\", \"A:A A:B\")")
    ("rate", po::value<double>(&rate)->default_value(8e6),
     "set sampling rate")
    ("freq", po::value<double>(&freq)->default_value(2420e6),
     "set central frequency")
    ("gain", po::value<double>(&gain)->default_value(40),
     "set RF gain of the receiving chains")
    ("ant", po::value<std::string>(&antname)->default_value("TX/RX"),
     "select antenna on the frontend");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  // Print help message.
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return ~0;
  }

  // Create a multi-USRP device.
  std::cout << "Creating USRP device..." << std::endl;
  uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);
  std::cout << "USRP device created." << std::endl;

  // Select the RX subdevice first; this mapping affects all the settings.
  if (vm.count("subdev")) usrp->set_rx_subdev_spec(subdev);

  // Set the RX sample rate on all channels.
  std::cout << boost::format("-- Asking for RX rate %f MHz...") % (rate/1e6);
  std::cout << std::endl;
  usrp->set_rx_rate(rate);
  std::cout << boost::format("-- Actually got RX rate %f MHz") % (usrp->get_rx_rate()/1e6);
  std::cout << std::endl;

  // Lock motherboard clock to internal reference source and reset time register.
  std::cout << "-- Setting device timestamp to 0.0 s... ";
  usrp->set_clock_source("internal");
  usrp->set_time_now(uhd::time_spec_t(0.0));
  std::cout << "done" << std::endl;
  
  // Set the RX center frequency.
  std::cout << boost::format("-- Asking for freq %f MHz...") % (freq/1e6) << std::endl;
  uhd::tune_request_t tunereq(freq);
  usrp->set_rx_freq(tunereq);
  std::cout << boost::format("-- Actually got freq %f MHz") % (usrp->get_rx_freq()/1e6) << std::endl;

  // Set the RX gain.
  std::cout << boost::format("-- Asking for gain %f dB...") % gain << std::endl;
  usrp->set_rx_gain(gain);
  std::cout << boost::format("-- Actually got gain %f dB") % usrp->get_rx_gain() << std::endl;

  // Set the RX antennas.
  std::cout << boost::format("-- Asking for antenna \"%s\"...") % antname << std::endl;
  usrp->set_rx_antenna(antname);
  std::cout << boost::format("-- Actually got antenna \"%s\"...") % antname << std::endl;

  // Select the RX channels.
  std::vector<std::string> channel_strings;
  std::vector<size_t> channel_nums;
  boost::split(channel_strings, channel_list, boost::is_any_of("\"',"));
  for (size_t c = 0; c < channel_strings.size(); c++) {
    size_t chan = boost::lexical_cast<int>(channel_strings[c]);
    if (chan < usrp->get_rx_num_channels()) {
      channel_nums.push_back(boost::lexical_cast<int>(channel_strings[c]));
    } else {
      throw std::runtime_error("Invalid channel(s) specified.");
    }
  }

  // Create a receive streamer.
  uhd::stream_args_t stream_args("fc64","sc16");
  stream_args.channels = channel_nums;
  uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

  // Allocate data buffers.
  const double nseconds = 2.0;
  const double bufsize_us = nseconds * rate;
  const size_t usrp_bufsize = rx_stream->get_max_num_samps();
  const size_t usrp_nbuffers = (size_t) ceil(bufsize_us / (float) usrp_bufsize);
  std::vector<iqsamp_t> bankA(usrp_bufsize * usrp_nbuffers);
  std::vector<iqsamp_t> bankB(usrp_bufsize * usrp_nbuffers);
  iqsamp_t *bankptr;

  // Auxiliary data for recv().
  uhd::rx_metadata_t md;
  double timeout = 10.0;
  unsigned long num_acc_samps = 0;
  unsigned long num_rx_samps = 0;

  // Setup streaming options.
  double stream_delay = 1.0;
  uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
  stream_cmd.num_samps = total_num_samps;
  stream_cmd.stream_now = false;
  stream_cmd.time_spec = uhd::time_spec_t(stream_delay);

  // Initialise frontend spinlocks.
  pthread_spin_init(&lock[0], PTHREAD_PROCESS_PRIVATE);
  pthread_spin_init(&lock[1], PTHREAD_PROCESS_PRIVATE);

  // Create thread for processing data.
  pthread_t proc_thread;
  proc_pars_t proc_pars = {&bankA.front(), &bankB.front(), usrp_bufsize/decfactor * usrp_nbuffers};
  pthread_create(&proc_thread, NULL, proc_routine, &proc_pars);

  std::signal(SIGINT, &sigint_handler);

  // Start streaming.
  rx_stream->issue_stream_cmd(stream_cmd);
  std::cout << std::endl << "Streaming... Press CTRL+C to stop." << std::endl;
  std::cout << std::endl << "      Timestamp       LAP    Info " << std::endl;
  
  while (stopsig == false) {

    bufselect = (bufselect + 1) % 2;
    
    // Lock buffer.
    pthread_spin_lock(&lock[bufselect]);
    
    // Claim buffer for writing.
    if (bufselect) {
      bankptr = &bankA.front();
    } else {
      bankptr = &bankB.front();
    }
      
    // Receive data.
    for (unsigned int n = 0; n < usrp_nbuffers; n++) {
      num_rx_samps = rx_stream->recv(bankptr + usrp_bufsize*n, usrp_bufsize, md, timeout, false);
      
      // Throw an exception on error.
      if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW) {
        std::cerr << std::endl << "OVERFLOW DETECTED" << std::endl;
      } else if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
        std::cout << boost::format("\nReceived %u samples.") % num_acc_samps;
        std::cout << std::endl;
        throw std::runtime_error(md.strerror());
      }
      
      // Increment number of received samples.
      num_acc_samps += num_rx_samps;

    }

    // Unlock buffer.
    pthread_spin_unlock(&lock[bufselect]);  
  }

  pthread_join(proc_thread, NULL);
  
  pthread_spin_destroy(&lock[0]);
  pthread_spin_destroy(&lock[1]);
  std::cout << "Done." << std::endl;

  return 0;
}
