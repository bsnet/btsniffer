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

#ifndef BTDECODER_H
#define BTDECODER_H

#include <stdint.h>
#include <vector>
#include <functional>
#include <unordered_map>
#include <sys/time.h>

typedef enum {
    LAP_STATE_NEW = 0,
    LAP_STATE_BRUTE_FORCING = 1,
    LAP_STATE_RESOLVED = 2,
} lap_status_t;


class lap_node
{
public:
	lap_node();
	lap_node(uint32_t lap);
	// bool is_bruteforce_done(void);
	// void set_bruteforce_done(void);
	void set_tstart(long long t) {tstart = t;}
	long long get_tstart() {return tstart;}
	void dump(void);
    lap_status_t get_status(void);
	void set_status(lap_status_t);
	void set_uap_data(int index, bool valid, uint32_t uap, uint32_t clk1, uint32_t clk2, int clk_index);
	void get_uap_data(int index, bool *valid, uint32_t *uap, uint32_t *clk1, uint32_t *clk2, int *clk_index);
	void set_ts(long long);
	long long get_ts();
	void increase_processed_packets (void);
	int get_processed_packets (void);
	void bf_failed (void);
	int get_bf_failures (void);
	void bf_cannot_init (void);
	int get_bf_noinit (void);
	void bf_noinfo (void);
	int get_bf_noinfo (void);
	void new_packet_too_close (bool);
	void get_packet_too_close (int *confirmed, int *notcofirmed);
	void count_broken_uap (void);
	int get_broken_uap (void);
	void count_packet_inthepast (void);
	int get_packet_inthepast (void);
private:
	lap_status_t state;
	uint32_t lap;
	struct uap_data {
		bool uap_valid;
		uint32_t uap;
		uint32_t clk1;
		uint32_t clk2;
		int clk_index;
	};
	struct uap_data uaps[32];
	long long ts;
	long long prevts;
	bool brute_forced;
	int processed_packets;
	int brute_force_failures;
	int brute_force_noinit;
	int brute_force_noinfo;
	int packet_too_close_confirmed;
	int packet_too_close_notconfirmed;
	int broken_uap;
	int packet_inthepast;
	long long tstart;
};

#endif
