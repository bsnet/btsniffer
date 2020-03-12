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

#include <string>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <pthread.h>
#include <csignal>
#include <stdio.h>
#include <unordered_map>
#include "lapnode.hpp"

lap_node::lap_node()
{
	lap = 0;
	memset(uaps, 0, sizeof(uaps));
	state = LAP_STATE_NEW;
	ts = 0;
	processed_packets = 0;
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

void lap_node::set_uap_data(int index, bool _valid, uint32_t _uap,
                            uint32_t _clk1, uint32_t _clk2, int _clk_index)
{
	if(index < 0 || index > 32) {
		std::cerr << "ERROR: Invalid UAP index " << index
		          << ", aborting..." << std::endl;
		exit(1);
	}
	uaps[index].uap_valid = _valid;
	uaps[index].uap = _uap;
	uaps[index].clk1 = _clk1;
	uaps[index].clk2 = _clk2;
	uaps[index].clk_index = _clk_index;
}

void lap_node::get_uap_data(int index, bool *_valid, uint32_t *_uap,
                            uint32_t *_clk1, uint32_t *_clk2, int *_clk_index)
{
	if(index < 0 || index > 32) {
		std::cerr << "ERROR: Invalid UAP index " << index
		          << ", aborting..." << std::endl;
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
	if (confirmed) packet_too_close_confirmed ++;
	else packet_too_close_notconfirmed ++;
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
