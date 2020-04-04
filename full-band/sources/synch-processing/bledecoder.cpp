/*
 *  Copyright 2019-2020 Marco Cominelli
 *  Copyright 2017-2020 Francesco Gringoli
 *  Copyright 2014      Omri Iluz
 *  Copyright 2012      Jiang Wei
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

#include "bledecoder.hpp"

#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include <math.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <assert.h>

#define RB(l) rb_buf[(rb_head+(l))%RB_SIZE]
#define Q(l) Quantize(l)
#define RB_SIZE 20000

bool running = true;

void handler(int arg)
{
    running = false;
}

BLESDR::BLESDR()
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

    uint32_t kk;
    for(kk = 0; kk < 256; kk ++) {
        uint8_t a = kk & 0xFF;
        uint8_t swapped = (uint8_t)(((a * 0x0802LU & 0x22110LU) | (a * 0x8020LU & 0x88440LU)) * 0x10101LU >> 16);
        swapmap[kk] = ((uint32_t) swapped & 0xFF);
    }
}

BLESDR::BLESDR(int channel)
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

int BLESDR::delay = 0;
int BLESDR::syncsn[2] = {-1, -1};
int BLESDR::syncts[2] = {-1, -1};
bool BLESDR::calibrated = false;

int32_t BLESDR::ExtractThreshold(void) {
    int32_t threshold = 0;
    int c;
    for (c = 0; c < 8 * srate; c++) {
        threshold += (int32_t)RB(c);
    }
    return (int32_t)threshold / (8 * srate);
}

void BLESDR::setChannel(int channel)
{
    chan = channel;
}

void BLESDR::setCarrier(float _carrier)
{
    carrier = _carrier;
}

BLESDR::~BLESDR() {
}

void BLESDR::RB_init(void) {
    rb_buf = (int16_t *)malloc(RB_SIZE * 2);
    memset(rb_buf, RB_SIZE * 2, 0);
} 
void BLESDR::RB_inc(void) {
    rb_head++;
    rb_head = (rb_head) % RB_SIZE;
}

inline bool BLESDR::Quantize(int16_t l) {
    return RB(l * srate) > g_threshold;
}

uint8_t BLESDR::SwapBits(uint8_t a) {
    //return (uint8_t)(((a * 0x0802LU & 0x22110LU) | (a * 0x8020LU & 0x88440LU)) * 0x10101LU >> 16);
    return (uint8_t) swapmap[a];
}


void BLESDR::ExtractBytes(int l, uint8_t* buffer, int count) {
    int t;
    for (t = 0; t < count; t++) {
        buffer[t] = ExtractByte(l + t * 8);
    }
}

uint8_t BLESDR::ExtractByte(int l) {
    uint8_t byte = 0;
    int c;
    for (c = 0; c < 8; c++) byte |= Q(l + c) << (7 - c);
    return byte;
}

bool BLESDR::DetectPreamble(void) {
    int transitions = 0;
    int c;

    /* preamble sequence is based on the 9th symbol (either 0x55555555 or 0xAAAAAAAA) */
    if (Q(9)) {
        for (c = 0; c < 8; c++) {
            transitions += Q(c) > Q(c + 1);
        }
        Q9 = 1;
    }
    else {
        for (c = 0; c < 8; c++) {
            transitions += Q(c) < Q(c + 1);
        }
        Q9 = 2;
    }
    return transitions == 4 && abs(g_threshold) < 15500;
}


void BLESDR::Receiver(int samples, int *isample, int *qsample, int processed_44MSsamples)
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
        if(samples_suspended > 2)
            samples_suspended --;
        else {
        srate = 4;
            g_threshold = ExtractThreshold();
            if (DetectPreamble()) {
                DecodeBLEPacket(1);
            }
        }

        samples_processed ++;

        last_phase = phase;
    }
}

#define MAX_CONN 10
uint32_t aa_found[MAX_CONN];
uint32_t crc_found[MAX_CONN];
int conn_found = 0;

int search_aa(uint32_t aa)
{
    int kk;
    for(kk = 0; kk < conn_found; kk ++) {
        if(aa_found[kk] == aa) {
            return kk;
        }
    }
    return -1;
}

int add_aa(uint32_t aa, uint32_t crc)
{
    if(conn_found == MAX_CONN) {
        fprintf(stderr, "Maximum number of connections found, skipping this one...\n");
        return -1;
    }
    aa_found[conn_found] = aa;
    crc_found[conn_found] = crc;
    conn_found ++;
    return conn_found - 1;
}

void dump_aa()
{
    int kk;
    for(kk = 0; kk < conn_found; kk ++) {
        printf("%d %08X - %06X\n", kk, (uint32_t) aa_found[kk], (uint32_t) crc_found[kk]);
    }
}


bool BLESDR::HandlePacket(mod_type_t mod, uint64_t packet_addr_l, uint8_t *packet_data,
                          int packet_length)
{
    uint8_t crc[3];
    uint32_t packet_addr;
    uint32_t packet_crc;
    uint32_t calced_crc;
    int c;

    btle_reverse_whiten(chan, packet_data, 2 + packet_length + 3);
    
    if (packet_addr_l == LE_ADV_AA) {  // Advertisement packet
        packet_addr = LE_ADV_AA;
        crc[0] = crc[1] = crc[2] = 0x55;
    }
    else {
        int index = search_aa(packet_addr_l);
        if(index == -1)
            return false;

        uint32_t pcrc = crc_found[index];
        crc[0] = (pcrc >> 16) & 0xFF;
        crc[1] = (pcrc >> 8) & 0xFF;
        crc[2] = (pcrc >> 0) & 0xFF;
    }

    /* calculate packet crc */
    calced_crc = btle_reverse_crc(packet_data, packet_length + 2, crc);

    packet_crc = 0;
    for (c = 0; c < 3; c++) packet_crc = (packet_crc << 8) | packet_data[packet_length + 2 + c];

    /* BLE valid packet found, dump information */
    if (packet_crc != calced_crc)
        return false;
    
    int i = 0;
    char text[17];
    uint8_t packet_data_sw[255];
    int totlen = packet_length + 2 + 3;
    int timenow = (samples_processed - RB_SIZE)/ (NYQFREQ_INT / DECIMA);

    // create swapped packet in memory
    for (i = 0; i < totlen; i ++) {
        packet_data_sw[i] = SwapBits(packet_data[i]);
    }

    // suspend processing all the samples composing this packet
    samples_suspended = (5 + 2 + packet_length + 3) * 8 * srate;

    // perform this only on the chosen channels IF this is an uncoded synch packet
    if(mod == MOD_1MB && issync == 1) {
        uint8_t synch_packet[] = {0x40, 0x14,                // header
                                  0x00, 0x11, 0x22, 0x33, 0x44, 0x55,    // adv mac address
        };

        if(!memcmp(packet_data_sw, synch_packet, sizeof(synch_packet))) {

        printf("synch packet\n");

            uint32_t sn = 0;
            for (i = 0; i < 4; i++) {
                sn <<= 8;
                sn |= packet_data_sw[i + 8]; 
            }
            if(carrier == CARRIER_LO) {
                syncsn[0] = sn;
                syncts[0] = timenow;
                if(sn == syncsn[1]) {
                    delay = syncts[0] - syncts[1];
                    printf("delay = %d\n", delay);
                    if(calibrated == false) {
                        printf("CALIBRATED!\n");
                        calibrated = true;
                    }
                }
            }
            else {
                syncsn[1] = sn;
                syncts[1] = timenow;
                if(sn == syncsn[0]) {
                    delay = syncts[0] - syncts[1];
                    printf("delay = %d\n", delay);
                    if(calibrated == false) {
                        printf("CALIBRATED!\n");
                        calibrated = true;
                    }
                }
            }
        }
    }

    printf("-----------------------------\n");
    int ptype = 0;
    switch(mod) {
    case MOD_1MB:
        printf("Modulation = 1Mb/s, ");
    break;
    default:
        assert(0);
    break;
    }
    
    printf("Chan %d, Carrier %.02f, Time = %dus, Q9 = %d, AA = %08X\n", chan, carrier, timenow, Q9, (uint32_t) packet_addr_l);
    for (i = 0; i < totlen; i++) {
        if((i % 16) == 0)
            printf("%04X: ", i);
        uint8_t ch = packet_data_sw[i];
        printf("0x%02X ", ch);
        if(isgraph(ch) || isalpha(ch) || isdigit(ch))
            text[i % 16] = ch;
        else
            text[i % 16] = '.';
        if(((i + 1) % 16) == 0) {
            text[16] = 0;
            printf("%s\n", text);
        }
    }
    text[(i % 16)] = 0;
    int jj = i;
    for( ; (jj % 16) != 0; jj ++)
        printf("     ");
    if((i % 16) != 0)
        printf("%s\n", text);

    printf("\n");

    if(chan >= 37 && (packet_data_sw[0] & 0xf) == 5) {
        uint64_t newaa = (packet_data_sw[14] << 0) |
            (packet_data_sw[15] << 8) |
            (packet_data_sw[16] << 16) |
            (packet_data_sw[17] << 24);
        uint32_t newcrc = (packet_data_sw[18] << 0) |
            (packet_data_sw[19] << 8) |
            (packet_data_sw[20] << 16);
        printf("FOUND -> newaa = %08X, newcrc = %06X\n", (uint32_t) newaa, (uint32_t) newcrc);
        add_aa(newaa, newcrc);
    }

    // maybe we can use the packet instead of dropping it
    // to enhance the reception probability at least in this channel
    if(droptraffic == 1) {
        return false;
    }

    int timeusec = timenow % 1000000;
    int timesec = (timenow - timeusec) / 1000000;
    struct timeval ts = {timesec, timeusec};
    if(carrier == CARRIER_LO) {
        printf("adding to CARRIER_LO, channel = %d, ts = %d.%06d\n",
            chan, (int) ts.tv_sec, (int) ts.tv_usec);
    } else if(carrier == CARRIER_HI) {
        int timenowdelayed = timenow + delay;
        int timeusec = timenowdelayed % 1000000;
        int timesec = (timenowdelayed - timeusec) / 1000000;
        struct timeval ts = {timesec, timeusec};
        printf("adding to CARRIER_HI, channel = %d, ts = %d.%06d (timenow = %d, delay = %d)\n",
            chan, (int) ts.tv_sec, (int) ts.tv_usec, timenow, delay);
    }
}


bool BLESDR::DecodeBLEPacket(int plen)
{
    int c;
    uint8_t packet_data[500];
    int packet_length;
    uint64_t packet_addr_l;
    uint32_t packet_addr;
    uint8_t packet_header_arr[2];

    // extract address
    packet_addr_l = 0;
    for (c = 0; c < 4; c++)
        packet_addr_l |= ((uint64_t)SwapBits(ExtractByte((c + plen) * 8))) << (8 * c);

    // prefiltering
    if (packet_addr_l != LE_ADV_AA && search_aa(packet_addr_l) == -1)
        return false;

    // extract pdu header: skip first five bytes (PREAMBLE + AA)
    ExtractBytes((4 + plen) * 8, packet_header_arr, 2);

    // whiten header only so we can extract pdu length
    btle_reverse_whiten(chan, packet_header_arr, 2);
    packet_length = SwapBits(packet_header_arr[1]) & 0xFF;

    if (packet_length < 2 || packet_length > 255)
        return false;

    // extract data
    ExtractBytes((4 + plen) * 8, packet_data, 2 + packet_length + 3);
    return HandlePacket(MOD_1MB, packet_addr_l, packet_data, packet_length);
}

void BLESDR::btle_reverse_whiten(uint8_t chan,uint8_t* data, uint8_t len) {

    uint8_t  i;
    uint8_t lfsr = SwapBits(chan) | 2;
    while (len--) {
        for (i = 0x80; i; i >>= 1) {

            if (lfsr & 0x80) {

                lfsr ^= 0x11;
                (*data) ^= i;
            }
            lfsr <<= 1;
        }
        data++;
    }
}

uint32_t BLESDR::btle_reverse_crc(const uint8_t* data, uint8_t len, uint8_t* dst) {

    uint8_t v, t, d;
    uint32_t crc = 0;
    while (len--) {

        d = SwapBits(*data++);
        for (v = 0; v < 8; v++, d >>= 1) {

            t = dst[0] >> 7;

            dst[0] <<= 1;
            if (dst[1] & 0x80) dst[0] |= 1;
            dst[1] <<= 1;
            if (dst[2] & 0x80) dst[1] |= 1;
            dst[2] <<= 1;


            if (t != (d & 1)) {

                dst[2] ^= 0x5B;
                dst[1] ^= 0x06;
            }
        }
    }
    for (v = 0; v < 3; v++) crc = (crc << 8) | dst[v];
    return crc;
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

    if(channel < 11)
        Fchan = 2404 + channel * 2;
    else if(channel < 37)
        Fchan = 2428 + (channel - 11) * 2;
    else if(channel == 37)
        Fchan = 2402;
    else if(channel == 38)
        Fchan = 2426;
    else if(channel == 39)
        Fchan = 2480;
    else {
        fprintf(stderr, "Invalid channel, defaulting to 37\n");
        Fchan = 2402;
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
            newreal +=
                ((int)decimator[jj]) * ((int)(playbuffer_real[index1] + playbuffer_real[index2]));
            newimag +=
                ((int)decimator[jj]) * ((int)(playbuffer_imag[index1] + playbuffer_imag[index2]));
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
    // created by matlab command "fir1(30, 1/22)"
    float firfloat[PLAYBUFSIZE] = {
        2.1800045e-03, 2.8429114e-03, 4.2951167e-03, 6.7272005e-03,
        1.0256546e-02, 1.4909062e-02, 2.0608805e-02, 2.7176859e-02,
        3.4339924e-02, 4.1748120e-02, 4.9000605e-02, 5.5676774e-02,
        6.1370255e-02, 6.5722563e-02, 6.8453284e-02, 6.9383941e-02,
        6.8453284e-02, 6.5722563e-02, 6.1370255e-02, 5.5676774e-02,
        4.9000605e-02, 4.1748120e-02, 3.4339924e-02, 2.7176859e-02,
        2.0608805e-02, 1.4909062e-02, 1.0256546e-02, 6.7272005e-03,
        4.2951167e-03, 2.8429114e-03, 2.1800045e-03,
    };

    // with the chosen filter, maximum after multiplying by 450000 is 31223
    for(int kk = 0; kk < 31; kk ++) {
    fir[kk] = int16_t(firfloat[kk] * 450000.0);
    }
}

int chan2index(int chan)
{
    int index = -1;
    if(chan == 37) {
        index = 0;
    } else if(chan >= 0 && chan <= 10) {
        index = chan + 1;
    } else if(chan == 38) {
        index = 12;
    } else if(chan == 39) {
        index = 39;
    } else {
        index = chan + 2;
    }

    assert(index != -1);

    return index;
}


int index2chan(int index)
{
    int chan = -1;
    if(index == 0) {
        chan = 37;
    } else if(index <= 11) {
        chan = index - 1;
    } else if(index == 12) {
        chan = 38;
    } else if(index == 39) {
        chan = 39;
    } else {
        chan = index - 2;
    }

    assert(chan != -1);

    return chan;
}


int main(int argc, char *argv[])
{
    int delay = 0;
    if (argc >= 4)
        delay = atoi (argv[3]);
    fprintf (stdout, "delay = %d\n", delay);

    int skip = 0;
    if (argc >= 5)
        skip = atoi (argv[4]);
    fprintf (stdout, "skip = %d\n", skip);

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

    BLESDR ble[41];
    Channelizer chan[41];
#define LO_CHAN(ch) (ch)
#define HI_CHAN(ch) (ch + 1)
    int bci; // ble channel index
    int bc;  // ble channel
    // Low band
    bc = index2chan(19);
    ble[LO_CHAN(19)].setChannel(bc);
    ble[LO_CHAN(19)].setCarrier(CARRIER_LO);
    chan[LO_CHAN(19)].init(CARRIER_LO, bc, decimator);
    // High band
    bc = index2chan(19);
    ble[HI_CHAN(19)].setChannel(bc);
    ble[HI_CHAN(19)].setCarrier(CARRIER_HI);
    chan[HI_CHAN(19)].init(CARRIER_HI, bc, decimator);

    // channel 19 on LO spectrum used for synch, keep traffic
    ble[LO_CHAN(19)].issync = 1;
    ble[LO_CHAN(19)].droptraffic = 0;

    // channel 19 on HI spectrum used for synch, drop traffic (duplicate otherwise)
    ble[HI_CHAN(19)].issync = 1;
    ble[HI_CHAN(19)].droptraffic = 1;

    unsigned int newaa = 0xB8324717;
    unsigned int newcrc = 0xD149D1;
    add_aa(newaa, newcrc);

#define SAMPLES_PER_BLOCK (4096)
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
        }
        else {
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

    while(running) {

        int samples_l = fread(input_buffer_l, 2, SAMPLES_PER_BLOCK, fptr_l);
        int samples_h = fread(input_buffer_h, 2, SAMPLES_PER_BLOCK, fptr_h);

        // if we do not have the same number of samples per band we break
        if(samples_l <= 0 || samples_h <=0 || samples_l != samples_h)
            break;

        int8_t *complex_samples;
        int new_samples;

        // Set channel index to synch channel (default = 19)
        bci = 19;
        // process synch channel in LO band
        complex_samples = input_buffer_l;
        new_samples = chan[LO_CHAN(bci)].FeedSamples(samples_l,
            complex_samples, sampled_real, sampled_imag);
        if(new_samples > 0) {
            ble[LO_CHAN(bci)].Receiver(new_samples,
                sampled_real, sampled_imag, processed_44MSsamples);
        }

        // process synch channel in HI band
        complex_samples = input_buffer_h;
        new_samples = chan[HI_CHAN(bci)].FeedSamples(samples_l,
            complex_samples, sampled_real, sampled_imag);
        if(new_samples > 0) {
                ble[HI_CHAN(bci)].Receiver(new_samples,
                    sampled_real, sampled_imag, processed_44MSsamples);
        }

        processed_44MSsamples += samples_l;
    }

    return 0;
}
