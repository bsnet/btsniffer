/*
 * Copyright 2019-2020 Marco Cominelli
 * Copyright 2017-2020 Francesco Gringoli
 * Copyright 2014      Omri Iluz
 * Copyright 2012      Jiang Wei
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

#include "btdecoder.hpp"
#include <pcap.h>

#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include <math.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>

#include <netinet/ether.h>
#include <net/ethernet.h>

#include <unordered_map>

#define CARRIER_LO 2421.5
#define CARRIER_HI 2460.5

#define PLAYBUFSIZE 31
#define PI 3.141592653589793
#define NYQFREQ 44.0
#define NYQFREQ_INT 44
#define DECIMA 11
#define SAMPLES_PER_BLOCK (4096)

#define RB(l) rb_buf[(rb_head+(l))%RB_SIZE]
#define RBE(l) rb_buf_energy[(rb_head+(l))%RB_SIZE]
#define Q(l) Quantize(l)
#define RB_SIZE 20000

bool running = true;
char *terminate_on_lap = NULL;
struct timeval t_start;
int codewords[32768];


void handler(int arg) { running = false; }


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


inline uint8_t SwapBits2(uint8_t a)
{
  return (uint8_t)
    (((a * 0x0802LU & 0x22110LU) | (a * 0x8020LU & 0x88440LU)) * 0x10101LU >> 16);
}



void fec_init(void)
{
  int kk;
  for (kk = 0; kk < 32768; kk ++)
    codewords[kk] = -1;

  int data;
  for (data = 0; data < 1024; data ++) {
    uint32_t fec = 0;
    int data_bit_jj;
    for (data_bit_jj = 0; data_bit_jj < 10; data_bit_jj ++) {
      int data_in = 0;
      if ((data >> data_bit_jj) & 0x1)
        data_in = 1;
      uint32_t fec_out = (fec >> 4) & 0x1;
      uint32_t data_out = data_in ^ fec_out;
      uint32_t fec_adder =
        (data_out << 4) |
        (data_out << 2) |
        (data_out << 0);
      fec = (fec << 1) & 0x1f;
      fec = fec ^ fec_adder;
    }
    fec = (SwapBits2(fec) >> 3);
    uint32_t codeword = (fec << 10) | data;
    if (codeword & ~0x7fff) {
      fprintf(stderr, "Error during hamming map generation\n");
      exit(1);
    }
    if (codewords[codeword] != -1) {
      fprintf(stderr, "map already filled\n");
      exit(1);
    }
    codewords[codeword] = data;
    int codeword_bit_jj;
    for (codeword_bit_jj = 0; codeword_bit_jj < 15; codeword_bit_jj ++) {
      uint32_t codeword_err = codeword ^ (1 << codeword_bit_jj);
      if (codeword_err & ~0x7fff) {
        fprintf(stderr, "Error during hamming map generation(2)\n");
        exit(1);
      }
      if (codewords[codeword_err] != -1) {
        fprintf(stderr, "map already filled(2)\n");
        exit(1);
      }
      codewords[codeword_err] = data;
    }
  }
}


bool crc16(uint32_t uap, int length, uint8_t *buffer)
{
  uint32_t crc = uap;
  int byte_kk, bit_jj;
  printf("Computing CRC16 step-by-step\n");
  for (byte_kk = 0; byte_kk < length; byte_kk ++) {
    for (bit_jj = 0; bit_jj < 8; bit_jj ++) {
      uint32_t data_in = (buffer[byte_kk] >> bit_jj) & 0x1;
      uint32_t crc_out = (crc >> 15) & 0x1;
      uint32_t crc_in = crc_out ^ data_in;
      crc = (crc << 1) & 0xffff;
      uint32_t crc_adder =
        (crc_in << 12) |
        (crc_in << 5) |
        (crc_in << 0);
      crc = crc ^ crc_adder;
    }
    printf("with byte %d: %04hX\n", byte_kk, crc);
  }
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


inline int length (uint64_t word) { return (_length (word, 64, 0)); }

uint64_t compute_remainder (uint64_t input, uint64_t divisor)
{
  int divisor_length = length (divisor);
  int input_length = length (input);

  if (divisor_length + input_length > 63) {
    return input;
  }

  // line up input
  input = input << divisor_length;

  while (length (input) >= divisor_length) {
    uint64_t tmp = divisor << (length (input) - divisor_length);
    input = input ^ tmp;
  }

  return input;
}

std::unordered_map<uint32_t, lap_node> lap_map;

bool find_list = false;
std::unordered_map<uint32_t, int> lap_to_find;
int number_of_lap_to_find;

lap_node::lap_node()
{
  lap = 0;
  memset(uaps, 0, sizeof(uaps));
  state = LAP_STATE_NEW;
  ts = 0;
  processed_packets = 0;
  brute_forced = false;
  brute_force_failures = 0;
  brute_force_noinit = 0;
  brute_force_noinfo = 0;
  packet_too_close_confirmed = 0;
  packet_too_close_notconfirmed = 0;
  broken_uap = 0;
  packet_inthepast = 0;
}

lap_node::lap_node(uint32_t _lap): lap_node()
{
  lap = _lap;
}

void lap_node::set_status(lap_status_t _state)
{
  state = _state;
}

lap_status_t lap_node::get_status(void)
{
  return state;
}

void lap_node::set_uap_data(int index, bool _valid, uint32_t _uap, uint32_t _clk1, uint32_t _clk2, int _clk_index)
{
  if(index < 0 || index > 32) {
    fprintf(stderr, "Invalid uap index, aborting...\n");
    exit(1);
  }
  uaps[index].uap_valid = _valid;
  uaps[index].uap = _uap;
  uaps[index].clk1 = _clk1;
  uaps[index].clk2 = _clk2;
  uaps[index].clk_index = _clk_index;
}

void lap_node::get_uap_data(int index, bool *_valid, uint32_t *_uap, uint32_t *_clk1, uint32_t *_clk2, int *_clk_index)
{
  if(index < 0 || index > 32) {
    fprintf(stderr, "Invalid uap index, aborting...\n");
    exit(1);
  }
  *_valid = uaps[index].uap_valid;
  *_uap = uaps[index].uap;
  *_clk1 = uaps[index].clk1;
  *_clk2 = uaps[index].clk2;
  *_clk_index = uaps[index].clk_index;
}

void lap_node::set_ts(long long _ts)
{
  prevts = ts;
  ts = _ts;
}

long long lap_node::get_ts(void)
{
  return ts;
}

void lap_node::increase_processed_packets (void)
{
  processed_packets ++;
}

int lap_node::get_processed_packets (void)
{
  return processed_packets;
}

void lap_node::bf_failed (void)
{
  brute_force_failures ++;
}

int lap_node::get_bf_failures (void)
{
  return brute_force_failures;
}

void lap_node::bf_cannot_init (void)
{
  brute_force_noinit ++;
}

int lap_node::get_bf_noinit (void)
{
  return brute_force_noinit;
}

void lap_node::bf_noinfo (void)
{
  brute_force_noinfo ++;
}

int lap_node::get_bf_noinfo (void)
{
  return brute_force_noinfo;
}

void lap_node::new_packet_too_close (bool confirmed)
{
  if (confirmed)
    packet_too_close_confirmed ++;
  else
    packet_too_close_notconfirmed ++;
}

void lap_node::get_packet_too_close (int *confirmed, int *notconfirmed)
{
  *confirmed = packet_too_close_confirmed;
  *notconfirmed = packet_too_close_notconfirmed;
}

void lap_node::count_broken_uap (void)
{
  broken_uap ++;
}

int lap_node::get_broken_uap (void)
{
  return broken_uap;
}

void lap_node::count_packet_inthepast (void)
{
  packet_inthepast ++;
}

int lap_node::get_packet_inthepast (void)
{
  return packet_inthepast;
}

void lap_node::set_bruteforce_done()
{
  brute_forced = true;
}

void lap_node::dump (void)
{
  fprintf(stdout, "lap = %u\n", lap);
  fprintf(stdout, "brute_forced = %d\n", (int) brute_forced);
  fprintf(stdout, "ts = %llX\n", (long long) ts);
  fprintf(stdout, "lap = %06X\n", lap);
  fprintf(stdout, "%d %d %d %d %d %d %d %d\n",
      processed_packets,
      brute_force_failures,
      brute_force_noinit,
      brute_force_noinfo,
      packet_too_close_confirmed,
      packet_too_close_notconfirmed,
      broken_uap,
      packet_inthepast);
}


typedef struct __attribute__((__packed__)) bdr_pseudoheader {
  uint8_t rf_channel;
  uint8_t PAD1[3];
  uint8_t transport_rate;
  uint8_t PAD2[3];
  uint32_t lap;
  float debug_energy;
  uint32_t whitened_packet_header;
  uint16_t flags;
} bdr_pseudoheader_t;

pcap_t *dev = NULL;
pcap_dumper_t *capfile = NULL;

#define SNAPLEN 3800

typedef struct packet_buffer {
  long long ts;
  uint8_t chan;
  float energy;
  uint8_t rate;
  uint32_t lap;
  uint32_t wph;
  int payload_length;
  uint8_t payload[SNAPLEN];
  struct packet_buffer *next;
  struct packet_buffer *prev;
} packet_buffer_t;

packet_buffer_t *head = NULL;
packet_buffer_t *last = NULL;

void packet_add(long long ts, uint8_t chan, float energy,
    uint8_t rate, uint32_t lap, uint32_t wph,
    int payload_length, uint8_t *payload)
{
#define TOLERANCE_US 10
  packet_buffer_t **cur = &head;
  while (*cur) {
    // check if ts are similar
    long long curts = (*cur)->ts;
    if (ts > curts - TOLERANCE_US &&
      ts < curts + TOLERANCE_US) {
      fprintf(stdout, "found existing packet in tolerance, checking energy\n");
      if (energy > (*cur)->energy) {
        fprintf(stdout, "replacing packet with higher energy\n");
        (*cur)->ts = ts;
        (*cur)->chan = chan;
        (*cur)->energy = energy;
        (*cur)->rate = rate;
        (*cur)->lap = lap;
        (*cur)->wph = wph;
        (*cur)->payload_length = payload_length;
        memcpy((*cur)->payload, payload, payload_length);
      } else {
        fprintf(stdout, "don't add packet with lower energy\n");
      }
      return;
    }

    if (ts < curts) {
      fprintf(stdout, "adding new packet before existing one\n");
      break;
    }

    cur = & ((*cur)->next);
  }

  if (! *cur) fprintf(stdout, "adding new packet at tail\n");
  packet_buffer_t *newpacket = new packet_buffer_t;
  newpacket->next = *cur;
  *cur = newpacket;
  (*cur)->ts = ts;
  (*cur)->chan = chan;
  (*cur)->energy = energy;
  (*cur)->rate = rate;
  (*cur)->lap = lap;
  (*cur)->wph = wph;
  (*cur)->payload_length = payload_length;
  memcpy((*cur)->payload, payload, payload_length);
}


void flush_bdr_packets(void)
{
  if (dev == NULL) {
    dev = pcap_open_dead(DLT_BLUETOOTH_BREDR_BB, SNAPLEN);
  }
  if (!dev) {
    fprintf(stderr, "Cannot open dead device\n");
    exit(1);
  }
  if (capfile == NULL) {
    capfile = pcap_dump_open(dev, "fica.pcap");
  }
  if (!capfile) {
    fprintf(stderr, "Cannot open file for writing\n");
    exit(1);
  }

  packet_buffer_t *cur = head;

  if (!cur) return;

  fprintf(stdout, "flushing packets to pcap\n");
  while(cur) {
    packet_buffer_t *next = cur->next;
    long long timenow_sec_us = cur->ts;
    int ts_sec = timenow_sec_us / 1000000;
    int ts_usec = timenow_sec_us % 1000000;

    uint8_t buffer[SNAPLEN + sizeof(bdr_pseudoheader_t)];
    bdr_pseudoheader_t *ph = (bdr_pseudoheader_t *) buffer;
    struct pcap_pkthdr pcaphdr;

    pcaphdr.ts.tv_sec = ts_sec;
    pcaphdr.ts.tv_usec = ts_usec;
    pcaphdr.len = sizeof(bdr_pseudoheader_t) + cur->payload_length;
    pcaphdr.caplen = pcaphdr.len;
    memset(ph, 0, sizeof(*ph));
    ph->rf_channel = cur->chan;
    ph->transport_rate = cur->rate;
    ph->lap = cur->lap;
    ph->whitened_packet_header = cur->wph;
    ph->flags = 0;
    memcpy(buffer + sizeof(*ph), cur->payload, cur->payload_length);  
    pcap_dump((u_char *) capfile, &pcaphdr, buffer);

    delete cur;
    cur = next;
  }

  head = NULL;
}


BTSDR::BTSDR()
{
  chan = 37;
  srate = 4;
  rb_head = -1;
  samples_suspended = 0;
  samples_processed = 0;
  RB_init();
  carrier = CARRIER_LO;
  issync = 0;
  droptraffic = 0;
  power = 0;

  uint32_t kk;
  for(kk = 0; kk < 256; kk ++) {
    uint8_t a = kk & 0xFF;
    uint8_t swapped = (uint8_t)(((a * 0x0802LU & 0x22110LU) | (a * 0x8020LU & 0x88440LU)) * 0x10101LU >> 16);
    swapmap[kk] = ((uint32_t) swapped & 0xFF);
  }
}


BTSDR::BTSDR(int channel)
{
  chan = channel;
  srate = 2;
  rb_head = -1;
  samples_suspended = 0;
  samples_processed = 0;
  RB_init();
  issync = 0;
  droptraffic = 0;

  uint32_t kk;
  for(kk = 0; kk < 256; kk ++) {
    uint8_t a = kk & 0xFF;
    uint8_t swapped = (uint8_t)(((a * 0x0802LU & 0x22110LU) | (a * 0x8020LU & 0x88440LU)) * 0x10101LU >> 16);
    swapmap[kk] = ((uint32_t) swapped & 0xFF);
  }
}


int BTSDR::delay = 0;
int BTSDR::syncsn[2] = {-1, -1};
int BTSDR::syncts[2] = {-1, -1};
bool BTSDR::calibrated = false;


void BTSDR::setChannel(int channel)
{
  chan = channel;
}

void BTSDR::setCarrier(float _carrier)
{
  carrier = _carrier;
}

BTSDR::~BTSDR() {
}

void BTSDR::RB_init(void) {
  rb_buf = (int16_t *)malloc(RB_SIZE * 2);
  rb_buf_energy = (float *) malloc(RB_SIZE * sizeof(float));
  memset(rb_buf, RB_SIZE * 2, 0);
}

void BTSDR::RB_inc(void) {
  rb_head++;
  rb_head = (rb_head) % RB_SIZE;
}

inline bool BTSDR::Quantize(int32_t l) {
  return RB(l * srate) > g_threshold;
}

inline float BTSDR::QE(int32_t l) {
  return RBE(l * srate);
}

uint8_t BTSDR::SwapBits(uint8_t a) {
  return (uint8_t) swapmap[a];
}


void BTSDR::ExtractBytes(int l, uint8_t* buffer, int count) {
  int t;
  for (t = 0; t < count; t++) {
    buffer[t] = ExtractByte(l + t * 8);
  }
}

uint8_t BTSDR::ExtractByte(int l) {
  uint8_t byte = 0;
  int c;
  for (c = 0; c < 8; c++) byte |= Q(l + c) << (7 - c);
  return byte;
}

// this function extracts a sequence of bits encoded using 2/3 FEC over DBPSK
int BTSDR::ExtractPayloadFec23(int l, int maxlen_byte, uint8_t *payload, uint32_t *whitener)
{
  int maxlen_bit = maxlen_byte * 8;
  int block_of_10 = maxlen_bit / 10;

  memset(payload, 0, maxlen_byte);

  if ((maxlen_bit % 10) != 0)
    block_of_10 ++;


  int errors = 0;

  int block_kk;
  for (block_kk = 0; block_kk < block_of_10; block_kk ++) {
    int block_bit_jj;
    uint32_t fec = 0;
    fprintf(stdout, "BLOCK: ");
    uint32_t codeword_rx = 0;
    for (block_bit_jj = 0; block_bit_jj < 10; block_bit_jj ++) {
      int bit_number = block_kk * 10 + block_bit_jj;
      int coded_bit_number = block_kk * 15 + block_bit_jj;
      int byte_number = bit_number / 8;
      int byte_bit_number = bit_number % 8;
      int data_in = 0;
      uint32_t bit_read = 0;
      if (Q(l + coded_bit_number)) {
        fprintf(stdout, "1");
        data_in = 1;
        bit_read = 1;
      } else {
        fprintf(stdout, "0");
      }
      codeword_rx |= (bit_read << block_bit_jj);

      uint32_t fec_out = (fec >> 4) & 0x1;
      uint32_t data_out = data_in ^ fec_out;
      uint32_t fec_adder =
        (data_out << 4) |
        (data_out << 2) |
        (data_out << 0);
      fec = (fec << 1) & 0x1f;
      fec = fec ^ fec_adder;
    }
    fprintf(stdout, "|");
    uint32_t fec_rx = 0;
    for (block_bit_jj = 10; block_bit_jj < 15; block_bit_jj ++) {
      int coded_bit_number = block_kk * 15 + block_bit_jj;
      uint32_t bit_read = 0;
      if (Q(l + coded_bit_number)) {
        bit_read = 1;
        fprintf(stdout, "1");
      } else {
        fprintf(stdout, "0");
      }
      fec_rx <<= 1;
      fec_rx |= bit_read;
      codeword_rx |= (bit_read << block_bit_jj);
    }
    fprintf(stdout, "=");
    int jj;
    uint32_t fec_copy = fec;
    for (jj = 0; jj < 5; jj ++) {
      if (fec_copy & 0x10) {
        fprintf(stdout, "1");
      } else {
        fprintf(stdout, "0");
      }
      fec_copy <<= 1;
    }
    fprintf(stdout, "|%02X==%02X", fec_rx, fec);

    if(fec != fec_rx) {
      fprintf(stdout, " E");
      // errors ++;
    } else {
      fprintf(stdout, " C");
    }

    fprintf(stdout, "|Codeword=%08X ", codeword_rx);
    if (codewords[codeword_rx] != -1) {
      fprintf(stdout, " C");

      // recompute payload, we should add a mechnism for
      // understanding we corrected a single error
      uint32_t correct_data = codewords[codeword_rx];
      for (block_bit_jj = 0; block_bit_jj < 10; block_bit_jj ++) {
        int bit_number = block_kk * 10 + block_bit_jj;
        // int coded_bit_number = block_kk * 15 + block_bit_jj;
        int byte_number = bit_number / 8;
        int byte_bit_number = bit_number % 8;
        int bit_whitened = 0;
        if (correct_data & (1 << block_bit_jj)) {
          bit_whitened = 1;
        }

        // dewhiten this bit
        uint32_t whitener_out = ((*whitener) >> 6) & 0x1;
        uint32_t whitener_shifted = ((*whitener) << 1) & 0x7f;
        *whitener = whitener_shifted ^ (whitener_out | (whitener_out << 4));
        payload[byte_number] |= (whitener_out << byte_bit_number);
      }
    } else {
      errors ++;
      fprintf(stdout, " E");
    }

    fprintf(stdout, "\n");
  }

  return errors;
}


// this function extracts header if FEC is perfect, then it bruteforce all possible clk
// to dewhiten the header and finally bruteforce uap until a valid hec code is found
int BTSDR::ExtractHeaderBruteForce(int l, uint32_t *head, uint32_t *clks, int *clk_found, uint32_t uap)
{
  int kk, jj;
  int errors = 0;
  uint32_t header = 0;
  for (kk = 0; kk < 54; kk += 3) {
    int s0 = 0;
    int s1 = 0;
    for (jj = 0; jj < 3; jj ++) {
      if (Q(l + kk + jj))
        s1 ++;
      else
        s0 ++;
    }
    header >>= 1;
    if (s1 == 0 || s0 == 0) errors ++;
    if (s1 > s0) {
      header |= 0x20000;
    }
  }

  *head = header;

  if (errors != 18) return -1;

  // first received bit is at the left (LSB) so header is composed of these fields
  // MSB 8bit |   10bit LSB
  //   +HEC+|...LT_ADDR

  // brute force the clock
  uint32_t clk;
  *clk_found = 0;
  for (clk = 0; clk < 64; clk ++) {
    uint32_t header_dewhiten = header;
    uint32_t whitener;
    whitener = (clk & 0x3f) | 0x40;

    // dewhiten the header with this clock
    for (kk = 0; kk < 18; kk ++) {
      // according to page 432 Figure 7.9 we use the MSB of the whitener
      // then we update it
      uint32_t whitener_out = (whitener >> 6) & 0x1;
      uint32_t whitener_shifted = (whitener << 1) & 0x7f;
      whitener = whitener_shifted ^ (whitener_out | (whitener_out << 4));

      // we dewhiten from LSB to MSB as LSB is transmitted first
      header_dewhiten = header_dewhiten ^ (whitener_out << kk);
    }

    // recompute the HEC over the dewhitened sequence and the known UAP
    uint32_t lsfr = uap;
    for (kk = 0; kk < 10; kk ++) {
      uint32_t lsfr_out = (lsfr >> 7) & 0x1;
      uint32_t data_in = (header_dewhiten >> kk) & 0x1;
      uint32_t lsfr_in = (lsfr_out ^ data_in);
      uint32_t lsfr_adder =
        (lsfr_in << 7) |
        (lsfr_in << 5) |
        (lsfr_in << 2) |
        (lsfr_in << 1) |
        (lsfr_in << 0);
      lsfr = (lsfr << 1) & 0xff;
      lsfr = lsfr ^ lsfr_adder;
    }

    // compare computed hec with received one
    // first received bit is in header_dewhiten bit 10 (from 0), last received in bit 17 (from 0)
    // computed first bit to be transmitted is in bit 7, last to be transmitted is in bit 0
    for (kk = 0; kk < 8; kk ++) {
      uint32_t bit_rx = (header_dewhiten >> (10 + kk)) & 0x1;
      uint32_t bit_tx = (lsfr >> (7 - kk)) & 0x1;
      if (bit_rx != bit_tx)
        break;
    }

    if (kk == 8) {
      clks[*clk_found] = clk;
      (*clk_found) ++;
    }
  }
  return *clk_found;
}

typedef enum {
  BDR_PACKET_NULL = 0,
  BDR_PACKET_POLL = 1,
  BDR_PACKET_DM1 = 3,
} bdr_packet_type;

bool BTSDR::DetectPreamble(void) {
  int transitions = 0;
  int c;

  if ((Q(4) > Q(3) && Q(2) > Q(1) && Q(0) > Q(1)) |
    (Q(3) > Q(4) && Q(1) > Q(2) && Q(1) > Q(0))) {

    // preambtsdr|   sync_word   |trailer
    // preambtsdr|code|lap|barker|trailer
    // 0     |4   |38 |62  |68
    uint64_t barker = SwapBits(ExtractByte(62));
    barker = barker & 0x3f;

    if(barker == 0x13 || barker == 0x2C) {
      // process access code to determine whether it is
      // CAC => trailer is present
      // DAC,GIAC,DIAC => trailer is present only with FHS packet
      // Standard says:
      // "When used as self-contained messages without a header,
      //  the DAC and IAC do not include the trailer bits and
      //  are of length 68 bits"
      // It also says:
      // "There is one general IAC (GIAC) for general inquiry
      //  operations and there are 63 dedicated IACs (DIACs) for
      //  dedicated inquiry operations."

      // device address is 48bit and is divided into
      // MSB    ---->    LSB
      // NAP(2B) UAP(1B) LAP(3B)

      // 1) extract lap
      uint64_t lap =
        SwapBits(ExtractByte(54)) << 16 |
        SwapBits(ExtractByte(46)) << 8 |
        SwapBits(ExtractByte(38));

      uint64_t code =
        (((uint64_t) SwapBits(ExtractByte( 4))) <<  0) |
        (((uint64_t) SwapBits(ExtractByte(12))) <<  8) |
        (((uint64_t) SwapBits(ExtractByte(20))) << 16) |
        (((uint64_t) SwapBits(ExtractByte(28))) << 24) |
        (((uint64_t) SwapBits(ExtractByte(36))) << 32);
      code = code & 0x3FFFFFFFFLLU;

      // access word
      uint64_t aw = (barker << 58) | (lap << 34) | code;

      // use lap to recreate aw from scratch, do not use barker, set it
      // accordingly to extracted lap
      uint64_t barker_true = 0;
      if ((lap & 0x800000) != 0) barker_true = 0x13;
      else barker_true = 0x2C;

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
          fprintf(stdout, "newlap = %06X on channel %d\n", _lap, chan);
          lap_map[_lap] = lap_node(_lap);
        }
      }
      else return 0;

#define INVALID_CLK_INDEX -1
#define DELTA_TS_SAME_THRESHOLD 40  // this should depend on frame length!
#define DELTA_TS_SLOT_THRESHOLD 620 // should be 625, give margin
#define SLOT_DURATION 625.0
#define ERROR_THRESHOLD 0.03
      uint32_t header = 0;
      uint32_t clk_tabtsdr[64];
      long long timenow_sec_us = (samples_processed - RB_SIZE) / (NYQFREQ_INT / DECIMA);

      lap_map[_lap].increase_processed_packets ();

      switch(lap_map[_lap].get_status()) {
        case LAP_STATE_NEW:
        {
          fprintf(stdout, "-------- %lld\n", timenow_sec_us);
          fprintf(stdout, "[%d] lap %X, trying to init\n", chan, _lap);
          int a, jj;
          uint32_t uap;
          lap_map[_lap].set_ts(timenow_sec_us);
          int valid_uaps = 0;
          for (uap = 0; uap < 256; uap ++) {
            int clk_found = ExtractHeaderBruteForce(72, &header, clk_tabtsdr, &a, uap);
            if (clk_found <= 0) {
              // don't log but keep trace for this lap
              lap_map[_lap].bf_failed ();
              continue;
            }
            else if (clk_found != 2) {
              // log this, it should never happen
              fprintf(stdout, "[%d] number of clk values invalid, aborting...\n", chan);
              exit(1);
            }
            lap_map[_lap].set_uap_data(valid_uaps, true, uap,
              clk_tabtsdr[0], clk_tabtsdr[1], INVALID_CLK_INDEX);
#ifdef DEBUGBDR
            fprintf(stdout, "[%d] UAP = %u: ", chan, uap);
            for (jj = 0; jj < clk_found; jj ++) fprintf(stdout, "%d, ", clk_tabtsdr[jj]);
            fprintf(stdout, "\n");
#endif
            valid_uaps ++;
          }
          if (valid_uaps != 32) {
            // if this happens we keep LAP_STATE_NEW, do not log but remember it
            fprintf(stdout, "[%d] cannot init\n", chan);
            lap_map[_lap].bf_cannot_init();
          } else {
            fprintf(stdout, "[%d] initialised\n", chan);
            lap_map[_lap].set_status(LAP_STATE_BRUTE_FORCING);
          }
          return 0;
        }
        case LAP_STATE_BRUTE_FORCING:
        {
          fprintf(stdout, "-------- %lld\n", timenow_sec_us);
      long long tmpdeltats = timenow_sec_us - lap_map[_lap].get_ts();
          if (tmpdeltats < 0) {
            fprintf(stdout, "[%d] found packet in the past for lap %X (delta = %lld), skipping...\n", chan, _lap, tmpdeltats);
            lap_map[_lap].count_packet_inthepast ();
            return 0;
          }
          fprintf(stdout, "[%d] again lap %X after %lld us, trying previously found UAP values\n", chan, _lap, tmpdeltats);
          if (llabs(lap_map[_lap].get_ts() - timenow_sec_us) < DELTA_TS_SAME_THRESHOLD) {
      // this might happen for packet received on side channels
            fprintf(stdout, "[%d] frame too close to another, just verify...\n", chan);
            int a, jj, clk_index;
            int confirmed_uap = 0;
            int valid_uap = 0;
            for (jj = 0; jj < 32; jj ++) {
              uint32_t clks[2];
              uint32_t uap;
              bool uap_valid;
              lap_map[_lap].get_uap_data(jj, &uap_valid, &uap, &clks[0], &clks[1], &clk_index);
              if(uap_valid == false)
                continue;
              valid_uap ++;
              int clk_found = ExtractHeaderBruteForce(72, &header, clk_tabtsdr, &a, uap);
              if (clk_found != 2)
                continue;
              if (clk_tabtsdr[0] != clks[0] || clk_tabtsdr[1] != clks[1])
                continue;
              confirmed_uap ++;
            }
            if (confirmed_uap == valid_uap)
              lap_map[_lap].new_packet_too_close (true);
            else
              lap_map[_lap].new_packet_too_close (false);
            fprintf(stdout, "[%d] confirmed %d over %d uaps, not doing anything on it\n", chan, confirmed_uap, valid_uap);
            return 0;
          }
          else if (llabs(lap_map[_lap].get_ts() - timenow_sec_us) < DELTA_TS_SLOT_THRESHOLD) {
            // we cannot handle this situation, so simply tell it happened and return
            fprintf(stdout, "[%d] frame could be generated by connection with same lap, cannot handle, aborting session...\n", chan);
            lap_map.erase(_lap);
            return 0;
          }
          else {
            // measure how far away we are from being a multiple of 625
            long long prevts = lap_map[_lap].get_ts();
            long long deltats = timenow_sec_us - prevts;
            float deltats_float = (float) deltats;
            float periods = deltats_float / 625.0;
            float periods_round = roundf(periods);
            float error = fabsf(periods - periods_round);
            if (error > ERROR_THRESHOLD) {
              fprintf(stdout, "[%d] error too big (%f), removing the session and restarting...\n", chan, error);
              lap_map[_lap].dump ();
              lap_map.erase(_lap);
              return 0;
            }
            int a, jj, clk_index;
            int count_valid_uap = 0;
            int count_broken_uap = 0;
            for (jj = 0; jj < 32; jj ++) {
              uint32_t clks[2];
              uint32_t uap;
              bool uap_valid;
              lap_map[_lap].get_uap_data(jj, &uap_valid, &uap, &clks[0], &clks[1], &clk_index);
              if(uap_valid == false)
                continue;
              int clk_found = ExtractHeaderBruteForce(72, &header, clk_tabtsdr, &a, uap);
              if (clk_found != 2) {
                // if we get here means either:
                // 1. connection changed, we should remove the old one and restart
                // 2. one frame was received corrupted
#ifdef DEBUGBDR
                fprintf(stdout, "[%d] broken uap %d\n", chan, uap);
#endif
                count_broken_uap ++;
                continue;
              }
              if (clk_tabtsdr[0] == clks[0] && clk_tabtsdr[1] == clks[1]) {
#ifdef DEBUGBDR
                fprintf(stdout, "[%d] cannot evaluate this uap, same clks as before, leaving untouched\n", chan);
#endif
                count_valid_uap ++;
                continue;
              }
              int slot = ((int) periods_round) % 64;
              // if we already found which one of the two is valid, then test only that one
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
                  int slot_guessed = (clk_tabtsdr[qnew] - clks[qold]) % 64;
                  if (slot == slot_guessed) {
                    qnew_valid = qnew;
                    matched ++;
                  }
                }
              }
              if (matched > 1) {
#ifdef DEBUGBDR
                fprintf(stdout, "[%d] Multiple valid combinations of clks for UAP %d, not deciding yet...\n", chan, uap);
#endif
                lap_map[_lap].set_uap_data(jj, true, uap, clk_tabtsdr[0], clk_tabtsdr[1], INVALID_CLK_INDEX);
                count_valid_uap ++;
              }
              else if (matched > 0) {
#ifdef DEBUGBDR
                fprintf(stdout, "[%d] keeping uap = %d (%d %d)\n", chan, uap, clk_tabtsdr[0], clk_tabtsdr[1]);
#endif
                lap_map[_lap].set_uap_data(jj, true, uap, clk_tabtsdr[0], clk_tabtsdr[1], qnew_valid);
                count_valid_uap ++;
              }
              else {
                lap_map[_lap].set_uap_data(jj, false, 0, 0, 0, INVALID_CLK_INDEX);
              }
            }
            if (count_valid_uap == 0 && count_broken_uap > 0) {
              lap_map[_lap].count_broken_uap ();
              fprintf(stdout, "[%d] This frame was likely broken, no uap changed, do not update ts\n", chan);
              return 0;
            }
            else if (count_valid_uap == 0) {
              fprintf(stdout, "[%d] No more valid uap for this lap, restarting\n", chan);
              lap_map[_lap].set_status(LAP_STATE_RESOLVED); // does this make sense??
              struct timeval tsdiff;
              struct timeval ts;
              gettimeofday(&ts, NULL);
              compute_tsdiff(&ts, &t_start, &tsdiff);
              fprintf(stdout, "[%d] UNRES %d.%d06 (%lld) lap = %06X\n", chan, (int) tsdiff.tv_sec, (int) tsdiff.tv_usec, timenow_sec_us, _lap);
              lap_map[_lap].dump ();
              lap_map.erase(_lap);
              return 0;
            }
            else if (count_valid_uap <= 2) {
              fprintf(stdout, "[%d] Only %d uap left, switching to state resolved\n", chan, count_valid_uap);
              lap_map[_lap].set_status(LAP_STATE_RESOLVED);
        lap_map[_lap].set_bruteforce_done ();
              int jj;
              fprintf(stdout, "\n\n***** RESOLVED %06X on channel %d\n", _lap, chan);
              for (jj = 0; jj < 32; jj ++) {
                uint32_t clks[2];
                uint32_t uap;
                bool isvalid;
                int clock_index;
                lap_map[_lap].get_uap_data(jj, &isvalid, &uap, &clks[0], &clks[1], &clock_index);
                struct timeval ts;
                gettimeofday(&ts, NULL);
                struct timeval tsdiff;
                compute_tsdiff(&ts, &t_start, &tsdiff);
                if (isvalid) {
                  fprintf(stdout, "RES %d.%d06 (%lld) lap = %06X, uap = %u, clks %u-%u, clk_index = %d, channel = %d\n", (int) tsdiff.tv_sec, (int) tsdiff.tv_usec, timenow_sec_us, _lap, uap, clks[0], clks[1], clock_index, chan);
                }
              }
              lap_map[_lap].dump ();
              fprintf(stdout, "\n\n*****************\n");
              if (terminate_on_lap) {
                uint32_t _terminate_on_lap = 0;
                sscanf(terminate_on_lap, "%X", &_terminate_on_lap);
                if(_terminate_on_lap == _lap) {
                  running = 0;
                }
              }
              else if (find_list) {
                if (lap_to_find.find(_lap) != lap_to_find.end ()) {
                  if (lap_to_find[_lap] != 1) {
                    fprintf (stderr, "it seems this lap was already found, aborting...\n");
                    exit (1);
                  }
                  fprintf (stdout, "lap %X hold %d in db to find, switching to 0\n",
                    _lap, lap_to_find[_lap]);
                  lap_to_find[_lap] = 0;
                  number_of_lap_to_find --;
                  if (number_of_lap_to_find == 0) {
                    fprintf (stdout, "all laps found, exiting...\n");
                    running = 0;
                  }
                  else {
                    fprintf (stdout, "still %d lap(s) to go\n", number_of_lap_to_find);
                  }
                }
              }
            }
            else {
              fprintf(stdout, "[%d] Found %d valid uap for this lap\n", chan, count_valid_uap);
            }
            lap_map[_lap].set_ts(timenow_sec_us);
            return 0;
          }
          return 0;
        }
        case LAP_STATE_RESOLVED:
        {
          // Don't do anything (for now)
          return 0;
        }
      }
    }
  }
  return transitions == 4 && abs(g_threshold) < 15500;
}

int32_t BTSDR::ExtractThreshold(void) {
  int32_t threshold = 0;
  int c;
  for (c = 0; c < 8 * srate; c++) {
    threshold += (int32_t)RB(c);
  }
  return (int32_t)threshold / (8 * srate);
}


void BTSDR::Receiver(int samples, int *isample, int *qsample, int processed_44MSsamples, pcapmerger &ps)
{
  int kk;
  for(kk = 0; kk < samples; kk ++) {
    double phase, dphase;
    phase = atan2(float(qsample[kk]), float(isample[kk]));
    dphase = phase - last_phase;

    if (dphase < -M_PI) dphase += 2 * M_PI;
    if (dphase > M_PI) dphase -= 2 * M_PI;

    uint16_t sample = uint16_t(dphase / M_PI*UINT16_MAX);
    RB_inc();
    RB(0) = (int)sample;
    RBE(0) = ((float) (qsample[kk])) * ((float) (qsample[kk])) + ((float) (isample[kk])) * ((float) (isample[kk]));

    if(samples_suspended > 2)
      samples_suspended --;
    else {
      srate = 4;
      g_threshold = ExtractThreshold();
    }

    samples_processed ++;

    last_phase = phase;
  }
}


class Channelizer {
public:
  Channelizer() : decima(0) {};
  void init(float, int, int16_t*);
  int FeedSamples(int samples, int8_t *complex_source, int *sampled_real, int *sampled_imag);

private:
  float F0;
  int Fchan;
  float DeltaF;
  int16_t *decimator;
  int16_t playbuffer_real[PLAYBUFSIZE];
  int16_t playbuffer_imag[PLAYBUFSIZE];
  int playbuffer_offset;
  int16_t SIN[NYQFREQ_INT];
  int16_t COS[NYQFREQ_INT];
  int decima;
  int sample;
  int sampleperiod;
  float TC;
};

void Channelizer::init(float Fcenter, int channel, int16_t *_decimator)
{
  F0 = Fcenter;

  Fchan = 2402 + channel;
  if (channel < 0 || channel > 79) {
    fprintf(stderr, "Invalid BDR channel\n");
    exit(-1);
  }

  memset(playbuffer_real, sizeof(playbuffer_real), 0);
  memset(playbuffer_imag, sizeof(playbuffer_imag), 0);

  DeltaF = Fchan - F0;
  TC = 1.0 / NYQFREQ;

  for(int kk = 0; kk < NYQFREQ_INT; kk ++) {
    float theta = -2.0 * PI * float(DeltaF) * float(kk) * TC;
    SIN[kk] = (int16_t) 32000.0 * sin(theta);
    COS[kk] = (int16_t) 32000.0 * cos(theta);
  }

  playbuffer_offset = 0;

  decimator = _decimator;
  decima = 0;
  sample = 0;
  sampleperiod = 0;
  TC = 1.0 / NYQFREQ;
}

int Channelizer::FeedSamples(int samples, int8_t *complex_source, int *sampled_real, int *sampled_imag)
{
  // (real + j * imag) * exp(j * theta)
  // (real + j * imag) * (cos(theta) + j * sin(theta))
  // real * cos(theta) - imag * sin(theta) +
  // j * (imag * cos(theta) + real * sin(theta)
  int sampled_samples = 0;
  int kk;

  for(kk = 0; kk < samples; kk ++) {
    int real2, imag2;
    if(DeltaF == 0) {
      real2 = int32_t(complex_source[2 * kk + 0]);
      imag2 = int32_t(complex_source[2 * kk + 1]);
    }
    else {
      real2 = (int32_t(complex_source[2 * kk + 0]) * int32_t(COS[sampleperiod]) -
           int32_t(complex_source[2 * kk + 1]) * int32_t(SIN[sampleperiod]));
      imag2 = (int32_t(complex_source[2 * kk + 1]) * int32_t(COS[sampleperiod]) +
           int32_t(complex_source[2 * kk + 0]) * int32_t(SIN[sampleperiod]));

      // worst case is
      // *) -128 * -32000 - 128 * -32000 = 8160000
      // after shifting values are in range [-250, 250]
      real2 = real2 >> 15;
      imag2 = imag2 >> 15;
    }

    sampleperiod ++;
    if(sampleperiod == NYQFREQ_INT)
      sampleperiod = 0;

    sample ++;
    playbuffer_real[playbuffer_offset] = int16_t(real2);
    playbuffer_imag[playbuffer_offset] = int16_t(imag2);

    int jj;
    int newreal = 0;
    int newimag = 0;

    // perform fir filtering
    for(jj = 0; jj < (PLAYBUFSIZE - 1) / 2; jj ++) {
      int index1 = (playbuffer_offset - jj);
      int index2 = (playbuffer_offset + 1 + jj);
      if(index1 < 0) index1 += PLAYBUFSIZE;
      if(index2 >= PLAYBUFSIZE) index2 -= PLAYBUFSIZE;
      newreal += ((int)decimator[jj]) * ((int)(playbuffer_real[index1] + playbuffer_real[index2]));
      newimag += ((int)decimator[jj]) * ((int)(playbuffer_imag[index1] + playbuffer_imag[index2]));
    }
    jj = (PLAYBUFSIZE - 1) / 2;
    int index1 = playbuffer_offset - jj;
    if(index1 < 0) index1 += PLAYBUFSIZE;
    newreal += ((int)decimator[jj]) * ((int)playbuffer_real[index1]);
    newimag += ((int)decimator[jj]) * ((int)playbuffer_imag[index1]);

    if((decima % DECIMA) == 0) {
      sampled_real[sampled_samples] = newreal;
      sampled_imag[sampled_samples] = newimag;
      sampled_samples ++;
    }

    decima ++;
    playbuffer_offset = (playbuffer_offset + 1) % PLAYBUFSIZE;
  }

  return sampled_samples;
}

void loadfiltercoefficients(int16_t *fir)
{
  // created with matlab command "fir1(30, 1/44)"
  float firfloat[PLAYBUFSIZE] = {
    4.1549144e-03, 4.8030862e-03, 6.5465541e-03, 9.3832073e-03,
    1.3248968e-02, 1.8019316e-02, 2.3514459e-02, 2.9507896e-02,
    3.5737940e-02, 4.1921526e-02, 4.7769538e-02, 5.3002719e-02,
    5.7367242e-02, 6.0648983e-02, 6.2685649e-02, 6.3376004e-02,
    6.2685649e-02, 6.0648983e-02, 5.7367242e-02, 5.3002719e-02,
    4.7769538e-02, 4.1921526e-02, 3.5737940e-02, 2.9507896e-02,
    2.3514459e-02, 1.8019316e-02, 1.3248968e-02, 9.3832073e-03,
    6.5465541e-03, 4.8030862e-03, 4.1549144e-03,
  };

  // with the chosen filter, maximum after multiplying by 450000 is 31223
  for(int kk = 0; kk < 31; kk ++) {
    fir[kk] = int16_t(firfloat[kk] * 450000.0);
  }
}

// TODO: can we turn thin into a macro?
// should we check this is a valid channel?
bool ischanlo(int chan)
{
  if(chan <= 39) {
    return true;
  }
  return false;
}


// Entry point
int main(int argc, char *argv[])
{
  setvbuf(stdout, NULL, _IONBF, 0);

  gettimeofday(&t_start, NULL);
  fprintf(stdout, "START %d.%d06\n", (int) t_start.tv_sec, (int) t_start.tv_usec);

  fec_init();

  int delay = 0;
  if (argc >= 4)
    delay = atoi (argv[3]);
  fprintf (stdout, "delay = %d\n", delay);

  int skip = 0;
  if (argc >= 5)
    skip = atoi (argv[4]);
  fprintf (stdout, "skip = %d\n", skip);

  if (argc >= 6) {
    if (argv[5][0] == 'L') {
      if (argc < 7) {
        fprintf (stderr, "Missing number of connections to find\n");
        exit (-1);
      }
      number_of_lap_to_find = atoi (argv[6]);
      if (number_of_lap_to_find < 1) {
        fprintf (stderr, "Invalid number of lap to find\n");
        exit (-1);
      }
      find_list = true;
      FILE *f2find = fopen(argv[5] + 1, "rt");
      if (!f2find) {
        fprintf (stderr, "Cannot open file with lap to find\n");
        return -1;
      }
      int lineread = 0;
      while (1) {
        char linebuf[256];
        if (fgets(linebuf, sizeof(linebuf), f2find) == NULL) break;
        lineread ++;
        struct ether_addr *addr = ether_aton (linebuf);
        if (addr == NULL) {
          fprintf (stderr, "Error in line %d\n", lineread);
          return -1;
        }
        uint32_t lap = addr->ether_addr_octet[3] << 16 |
                 addr->ether_addr_octet[4] << 8 |
                 addr->ether_addr_octet[5];
        lap_to_find[lap] = 1;
      }
      fclose (f2find);
      if (((int) lap_to_find.size ()) < number_of_lap_to_find) {
        fprintf (stderr, "Number of lap to find bigger than lap in list\n");
        exit (-1);
      }
      fprintf (stdout, "Read %d lap(s), searching for %d\n", (int) lap_to_find.size (), number_of_lap_to_find);
    } else {
      terminate_on_lap = argv[5];
    }
  }

  pcapmerger pc;
  if(pc.create("trace.pcap")) {
    fprintf(stderr, "Cannot create pcap saver\n");
    return -1;
  }

  signal(SIGINT, &handler);

  if(argc < 3) {
    fprintf(stderr, "Missing arguments\n");
    return -1;
  }

  FILE *fptr_l = NULL;
  if((fptr_l = fopen(argv[1], "rb")) == NULL) {
    fprintf(stderr, "ERRORE\n");
    return -1;
  }
  FILE *fptr_h = NULL;
  if((fptr_h = fopen(argv[2], "rb")) == NULL) {
    fprintf(stderr, "ERRORE\n");
    return -1;
  }

  int16_t decimator[PLAYBUFSIZE];
  loadfiltercoefficients(decimator);

  // Bluetooth channels are 79
  BTSDR btsdr[79];
  Channelizer chan[79];
  // Low part of the spectrum
  for(int bc = 0; bc <= 39; bc ++) {
    btsdr[bc].setChannel(bc);
    btsdr[bc].setCarrier(CARRIER_LO);
    chan[bc].init(CARRIER_LO, bc, decimator);
  }
  // High part of the spectrum
  for(int bc = 40; bc <= 78; bc ++) {
    btsdr[bc].setChannel(bc);
    btsdr[bc].setCarrier(CARRIER_HI);
    chan[bc].init(CARRIER_HI, bc, decimator);
  }

  int8_t input_buffer_l[SAMPLES_PER_BLOCK * 2];
  int8_t input_buffer_h[SAMPLES_PER_BLOCK * 2];
  int sampled_real[SAMPLES_PER_BLOCK];
  int sampled_imag[SAMPLES_PER_BLOCK];

  int processed_44MSsamples = 0;

  if (delay != 0) {
    if (delay > 0) {
      while (delay > 0) {
        int a = fread(input_buffer_l, 2, 44, fptr_l);
        delay--;
      }
    } else {
      while (delay < 0) {
        int a = fread(input_buffer_h, 2, 44, fptr_h);
        delay++;
      }
    }
  }

  if (skip != 0) {
    while (skip > 0) {
      int a = fread(input_buffer_l, 2, 44, fptr_l);
      int b = fread(input_buffer_h, 2, 44, fptr_h);
      skip--;
    }
  }

  long long cum_sampleslo = 0;
  long long cum_sampleshi = 0;

  int my_lap = 0x9F8C90; //0x9353D6;
  lap_map[my_lap] = lap_node(my_lap);
  lap_map[my_lap].set_status(LAP_STATE_RESOLVED);
  lap_map[my_lap].set_bruteforce_done ();

  while(running) {

    int samples_l = fread(input_buffer_l, 2, SAMPLES_PER_BLOCK, fptr_l);
    int samples_h = fread(input_buffer_h, 2, SAMPLES_PER_BLOCK, fptr_h);

    // if different number of samples per band, we quit
    if(samples_l <= 0 || samples_h <= 0 || samples_l != samples_h) break;

    int8_t *complex_samples;
    int new_samples;

    // Process Bluetooth channels
    for(int bc = 0; bc < 79; bc ++) {
      if(ischanlo(bc)) {
        complex_samples = input_buffer_l;
        new_samples = chan[bc].FeedSamples(samples_l, complex_samples, sampled_real, sampled_imag);
        if(new_samples > 0) {
          btsdr[bc].Receiver(new_samples, sampled_real, sampled_imag, processed_44MSsamples, pc);
        }
        if (bc == 0) cum_sampleslo = cum_sampleslo + new_samples;
      }
      else {
        complex_samples = input_buffer_h;
        new_samples = chan[bc].FeedSamples(samples_h, complex_samples, sampled_real, sampled_imag);
        if(new_samples > 0) {
          btsdr[bc].Receiver(new_samples, sampled_real, sampled_imag, processed_44MSsamples, pc);
        }
        if (bc == 78) cum_sampleshi = cum_sampleshi + new_samples;
      }
    }

    flush_bdr_packets();

    processed_44MSsamples += samples_l;
  }

  long long timenow_sec_us_lo = (cum_sampleslo - RB_SIZE) / (NYQFREQ_INT / DECIMA);
  long long timenow_sec_us_hi = (cum_sampleshi - RB_SIZE) / (NYQFREQ_INT / DECIMA);
  fprintf(stdout, "QQQQQQQQQQ FINAL TIME %lld %lld\n", timenow_sec_us_lo, timenow_sec_us_hi);
  fprintf(stdout, "********** FINAL DUMP\n");
  std::unordered_map<uint32_t, lap_node>::iterator it;
  for (it = lap_map.begin (); it != lap_map.end (); it ++) {
    fprintf(stdout, "---- Dumping %X\n", it->first);
    it->second.dump ();
  }

  pcap_dump_close(capfile);

  pc.flush();
  pc.finalise();

  return 0;
}
