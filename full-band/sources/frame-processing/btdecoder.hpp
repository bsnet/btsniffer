/*
 * Copyright 2019-2020 Marco Cominelli
 * Copyright 2017-2020 Francesco Gringoli
 * Copyright 2014 Omri Iluz
 * Copyright 2012 Jiang Wei
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

#include <stdint.h>
#include <vector>
#include <functional>
#include <unordered_map>

#include "pcapsaver.hpp"

typedef enum {
    LAP_STATE_NEW = 0,
    LAP_STATE_BRUTE_FORCING = 1,
    LAP_STATE_RESOLVED = 2,
} lap_status_t;

class lap_node {
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
  public:
    lap_node();
    lap_node(uint32_t lap);
    void set_bruteforce_done(void);
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
};

class BTSDR {
public:
	BTSDR();
	BTSDR(int channel);
	~BTSDR();

	void setChannel(int channel);
	void setCarrier(float carrier);
	void Receiver(int samples, int *isample, int *qsample, int processed_44MSsamples, pcapmerger &ps);
	void setDelay(int usec_delay);

	int issync; // to be removed
	int droptraffic;

private:
        // std::unordered_map<uint32_t, lap_node> lap_map;
	uint32_t	swapmap[256];

        static int syncsn[2];
        static int syncts[2];
	static int delay;
	static bool calibrated;

	float carrier;
	uint8_t chan;
	int32_t g_threshold; // Quantization threshold
	int srate;           // sample rate downconvert ratio
	double last_phase;

	int samples_suspended;
	long long samples_processed;

	int rb_head;
	int16_t *rb_buf;
        float *rb_buf_energy;
	/* Init Ring Buffer */
	void RB_init(void);
	/* increment Ring Buffer Head */
	void RB_inc(void);
	/* Access Ring Buffer location l */

	int Q9;

	uint8_t SwapBits(uint8_t a);

	bool Quantize(int32_t l);
	float QE(int32_t l);

	int32_t ExtractThreshold(void);

	bool DetectPreamble(void);
	bool DetectPreamble2M(void);
	bool DetectPreambleCODED(void);

	uint8_t inline ExtractByte(int l);

	void ExtractBytes(int l, uint8_t* buffer, int count);

	bool DecodeBTLEPacket(pcapmerger &ps, int rate);
	bool DecodeBTLEPacketCODED(pcapmerger &ps);

	void btle_reverse_whiten(uint8_t chan,uint8_t* data, uint8_t len);

	uint32_t btle_reverse_crc(const uint8_t* data, uint8_t len, uint8_t* dst);

	int ExtractHeaderBruteForce(int, uint32_t*, uint32_t*, int*, uint32_t);
	int ExtractHeader(int, uint32_t*, uint32_t, uint32_t*);
	int ExtractHeaderAndPayload(int, uint32_t*);
	bool VerifyHEC(uint32_t header, uint8_t uap);
	int ExtractPayloadFec23(int l, int maxlen_byte, uint8_t *payload, uint32_t *whitener);

	float power;
};
