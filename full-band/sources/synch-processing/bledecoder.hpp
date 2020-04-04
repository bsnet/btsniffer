/*
 * Copyright 2019-2020 Marco Cominelli
 * Copyright 2017-2020 Francesco Gringoli
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

#ifndef BLEDECODER_HPP
#define BLEDECODER_HPP

#include <stdint.h>
#include <vector>
#include <functional>

#define MAX_LE_SYMBOLS 64
#define PLAYBUFSIZE    31

#define LE_ADV_AA 0x8E89BED6
#define LE_ADV_AA_inv 0x6B7D9171

#define ADV_IND         0
#define ADV_DIRECT_IND  1
#define ADV_NONCONN_IND 2
#define SCAN_REQ        3
#define SCAN_RSP        4
#define CONNECT_REQ     5
#define ADV_SCAN_IND    6

#define CARRIER_LO 2421.5
#define CARRIER_HI 2460.5

#define PI 3.141592653589793
#define NYQFREQ 44.0
#define NYQFREQ_INT 44
#define DECIMA 11

typedef enum {
    MOD_1MB,
} mod_type_t;


class BLESDR {
public:
    BLESDR();
    BLESDR(int channel);
    ~BLESDR();

    void setChannel(int channel);
    void setCarrier(float carrier);
    void Receiver(int samples, int *isample, int *qsample, int processed_44MSsamples);
    void setDelay(int usec_delay);

    int issync; // to be removed
    int droptraffic;

private:
    uint32_t swapmap[256];

    static int syncsn[2];
    static int syncts[2];
    static int delay;
    static bool calibrated;

    float carrier;
    uint8_t chan;
    int32_t g_threshold; // Quantization threshold
    int srate; // sample rate downconvert ratio
    double last_phase;

    int samples_suspended;
    int samples_processed;

    int rb_head;
    int16_t *rb_buf;
    void RB_init(void); // init ring buffer
    void RB_inc(void);  // increment buffer head

    int Q9;

    uint8_t SwapBits(uint8_t a);
    bool Quantize(int16_t l);
    int32_t ExtractThreshold(void);
    bool DetectPreamble(void);
    uint8_t inline ExtractByte(int l);
    void ExtractBytes(int l, uint8_t* buffer, int count);
    bool DecodeBLEPacket(int rate);
    void btle_reverse_whiten(uint8_t chan,uint8_t* data, uint8_t len);
    uint32_t btle_reverse_crc(const uint8_t* data, uint8_t len, uint8_t* dst);
    bool HandlePacket(mod_type_t mod, uint64_t packet_addr_l, uint8_t *packet_data, int packet_length);
};

#endif //BLEDECODER_HPP
