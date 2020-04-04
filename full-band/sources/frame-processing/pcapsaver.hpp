/*
 * Copyright 2019-2020 Marco Cominelli
 * Copyright 2017-2020 Francesco Gringoli
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

#ifndef _PCAP_MERGER_
#define _PCAP_MERGER_

#include <pcap.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include <arpa/inet.h>
#include <stdint.h>
#include <math.h>

struct ble_hdr {
        uint8_t         channel;
        uint8_t         sig_power;  
        uint8_t         noise_power;
        uint8_t         offenses;
        uint32_t        ref_addr;
        uint16_t        flags;
} __attribute__((__packed__));

#define	BLE_DEWHITENED	(1 << 0)
#define	CRC_CHECKED	(1 << 10)
#define	CRC_VALID	(1 << 11)
#define BLE_MODULATION_MASK	(7 << 7)
#define BLE_2MB		(1 << 7)
#define BLE_CODED_125	(2 << 7)
#define	BLE_CODED_500	(4 << 7)

class packet {
private:
#define MAX_PAYLOAD 255
#define AA_LENGTH 4
	struct pcap_pkthdr pcaphdr;
	uint8_t blob[sizeof(struct ble_hdr) + AA_LENGTH + MAX_PAYLOAD];
	bool hasheader;
	bool haspayload;
	packet *next;
public:
	packet() : next(NULL), hasheader(false), haspayload(false) {};
	void setnext(packet *p);
	packet *getnext(void);
	void setheader(int channel, uint32_t aa, struct timeval *ts, uint32_t ptype);
	void setpayload(int length, uint8_t *payload);
	struct pcap_pkthdr *getpcaphdr(void);
	uint8_t *getdata(void);
	void setdelay(int us);
	friend bool operator<(const packet &p1, const packet &p2);
	struct timeval *getts(void);
};

class pcapmerger {
private:
	pcap_dumper_t *capfile;
	pcap_t *dev;
	packet *head1;
	packet *tail1;
	packet *head2;
	packet *tail2;
public:
	pcapmerger() : capfile(NULL), dev(NULL), head1(NULL),
		tail1(NULL), head2(NULL), tail2(NULL) {};
	int create(const char *filename);
	int addpacket1(struct timeval *ts, int channel,
		       uint32_t aa, int length, uint8_t *payload, int ptype);
	int addpacket2(struct timeval *ts, int channel,
		       uint32_t aa, int length, uint8_t *payload, int ptype);
	void finalise(void);
	void changedelay(int delay);
	void flush(void);
};

#endif
