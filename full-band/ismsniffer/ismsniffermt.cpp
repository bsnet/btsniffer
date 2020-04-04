//
// Copyright 2019-2020 Marco Cominelli
// Copyright 2017-2020 Francesco Gringoli
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <uhd/types/tune_request.hpp>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include <csignal>
#include <complex>
#include <pthread.h>
#include <stdlib.h>
#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>

#include "blegenerator.hpp"

#define noSTANDARD

#define HW_SAMPLE_SIZE 4 /* we want to receive sc16 from hardware, but now is sc8, so 2 bytes */
#define ALIGNMENT 64     /* just in case gcc wants to use AVX2 */
#define SAMPLE_LENGTH 500000
#define SAMP_MARGIN 3000

namespace po = boost::program_options;
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

struct rxthreadpars {
    uhd::usrp::multi_usrp::sptr *usrp;
    std::string file;
    size_t samps_per_buff;
    double total_time;
    bool continue_on_bad_packet;
    volatile int semA[2];
    volatile int semB[2];
    uint8_t *bankA = NULL;
    uint8_t *bankB = NULL;
    int samplesA, samplesB;
};

void pack_high8_baseline(uint8_t *__restrict__ dst, const uint16_t *__restrict__ src, size_t bytes)
{
    uint8_t *end_dst = dst + bytes;
    do {
        *dst++ = *src++ >> 8;
    } while(dst < end_dst);
}

void pack_high8_indexed_src(uint8_t *__restrict__ dst, const uint16_t *__restrict__ src, size_t bytes)
{
    uintptr_t end_dst = (uintptr_t)(dst + bytes);
    uintptr_t srcu = (uintptr_t)src, dstu = (uintptr_t)dst;

    ptrdiff_t src_dst_offset = srcu - 2*dstu;
    do {
        __m128i v0 = _mm_loadu_si128((__m128i*)(dstu*2+src_dst_offset));
        __m128i v1 = _mm_loadu_si128((__m128i*)(dstu*2+src_dst_offset)+1);
        // indexed addressing modes can't stay micro-fused with VPAND
        // so this doesn't help as much with AVX if using offset loads + AND
        v0 = _mm_srli_epi16(v0, 8);
        v1 = _mm_srli_epi16(v1, 8);
        __m128i pack = _mm_packus_epi16(v0, v1);
        _mm_storeu_si128((__m128i*)dstu, pack);
        dstu += 16;
    } while(dstu < end_dst);
}

void *store_thread_routine(void *pars)
{
    struct rxthreadpars *rxpar = (struct rxthreadpars *) pars;
    const std::string file = rxpar->file;
    std::ofstream outfile;
    outfile.open(file.c_str(), std::ofstream::binary);

    int kk = 0;
    int totsampA = 0;
    int totsampB = 0;

    uint8_t __attribute__ ((aligned (ALIGNMENT))) buffer[HW_SAMPLE_SIZE / 2 * SAMPLE_LENGTH];

    while(not stop_signal_called) {
        int bank = (kk % 2);
	while(__sync_lock_test_and_set(&rxpar->semB[bank], 1));
        uint8_t *curbank = (bank == 0 ? rxpar->bankA : rxpar->bankB);
        int samples = (bank == 0 ? rxpar->samplesA : rxpar->samplesB);
	if(bank == 0) totsampA += samples;
	else totsampB += samples;

	// reduce data by compressing sc16 into sc8
	pack_high8_indexed_src(buffer, (uint16_t *) curbank, samples * HW_SAMPLE_SIZE / 2);

if(samples != SAMPLE_LENGTH) printf("dd %d\n", samples);

        if (outfile.is_open())
            outfile.write((const char*) buffer, samples * HW_SAMPLE_SIZE / 2);
	__sync_lock_release(&rxpar->semA[bank]);
	kk ++;
    }

    if (outfile.is_open())
        outfile.close();

    printf("stored %d, %d\n", totsampA, totsampB);
    std::cout << "store thread finished" << std::endl;

    return NULL;
}

void *rx_thread_routine(void *pars)
{
    struct rxthreadpars *rxpar = (struct rxthreadpars *) pars;
    uhd::usrp::multi_usrp::sptr usrp = *(rxpar->usrp);
    const std::string cpu_format = "sc16";
    size_t samps_per_buff = rxpar->samps_per_buff;
    double time_requested = rxpar->total_time;
    bool continue_on_bad_packet = rxpar->continue_on_bad_packet;

    // sleep a little bit
    usleep(100000);

    //create a receive streamer
    uhd::stream_args_t stream_args(cpu_format, cpu_format);
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);
    uhd::rx_metadata_t md;
    bool overflow_message = true;

    //setup streaming
    uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
    size_t num_requested_samples = 0;
    stream_cmd.num_samps = num_requested_samples;
    stream_cmd.stream_now = true;
    stream_cmd.time_spec = uhd::time_spec_t();
    rx_stream->issue_stream_cmd(stream_cmd);

    boost::system_time start = boost::get_system_time();
    unsigned long long ticks_requested = (long)(time_requested * (double)boost::posix_time::time_duration::ticks_per_second());
    boost::posix_time::time_duration ticks_diff;
    boost::system_time last_update = start;
    unsigned long long last_update_samps = 0;

    typedef std::map<size_t,size_t> SizeMap;

    int kk = 0;
    int totsampA = 0;
    int totsampB = 0;

    // equivalent to sleep ~ 1 second
    int byte_to_skip = 40000000;

    while(not stop_signal_called) {
	int bank = (kk % 2);
	while(__sync_lock_test_and_set(&rxpar->semA[bank], 1));
	uint8_t *curbank = (bank == 0 ? rxpar->bankA : rxpar->bankB);
        boost::system_time now = boost::get_system_time();

	bool enable_size_map = false;

        int samp_in_buff = 0;

        while(samp_in_buff < samps_per_buff - SAMP_MARGIN) {
            size_t num_rx_samps = rx_stream->recv(curbank + samp_in_buff,
						  samps_per_buff - samp_in_buff,
						  md, 3.0, enable_size_map);
            if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
                std::cout << boost::format("Timeout while streaming") << std::endl;
                std::string error = str(boost::format("Receiver error: timeout"));
		throw std::runtime_error(error);
            }
            if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
                if (overflow_message) {
                    overflow_message = false;
                    std::cerr << boost::format(
                        "Got an overflow indication. Please consider the following:\n"
                        "  Your write medium must sustain a rate of %fMB/s.\n"
                        "  Dropped samples will not be written to the file.\n"
                        "  Please modify this example for your purposes.\n"
                        "  This message will not appear again.\n"
                    ) % (usrp->get_rx_rate()*sizeof(std::complex<char>)/1e6);
                }
            }
            else if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
                std::string error = str(boost::format("Receiver error: %s") % md.strerror());
                if (continue_on_bad_packet){
                    std::cerr << error << std::endl;
                }
                else
                    throw std::runtime_error(error);
            }

	    if(byte_to_skip > 0) {
                byte_to_skip -= num_rx_samps;
		continue;
	    }
	    samp_in_buff += num_rx_samps;
	}

        ticks_diff = now - start;
        if (ticks_requested > 0 &&
	    (unsigned long long)ticks_diff.ticks() > ticks_requested)
	    break;

	if(bank == 0) {
	    rxpar->samplesA = samp_in_buff;
	    totsampA += (int) samp_in_buff;
	}
	else {
	    rxpar->samplesB = samp_in_buff;
	    totsampB += (int) samp_in_buff;
	}

	__sync_lock_release(&rxpar->semB[bank]);
	kk ++;
    }

    stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
    rx_stream->issue_stream_cmd(stream_cmd);

    //finished
    printf("rxA = %d, rxB = %d\n", totsampA, totsampB);
    std::cout << std::endl << "rx stream done!" << std::endl << std::endl;

    return NULL;
}

struct txthreadpars {
    uhd::usrp::multi_usrp::sptr *usrp;
    size_t samps_per_buff;
    double msdelay;
};

void *tx_thread_routine(void *pars)
{
    struct txthreadpars *txpar = (struct txthreadpars *) pars;
    uhd::usrp::multi_usrp::sptr usrp = *(txpar->usrp);
    const std::string cpu_format = "sc8";
    size_t samps_per_buff = txpar->samps_per_buff;

    // create a transmit streamer
    uhd::stream_args_t stream_args(cpu_format, cpu_format);
    uhd::tx_streamer::sptr tx_stream = usrp->get_tx_stream(stream_args);

    // create the stream
    uhd::tx_metadata_t md;
    md.start_of_burst = false;
    md.end_of_burst = false;

    int counter = 0;
#define ZERO_PAD_BYTES 100

    do {
#define SYNC_CHANNEL 17
        std::vector<int8_t> buff = generatesamples(SYNC_CHANNEL, counter, ZERO_PAD_BYTES);
        size_t num_tx_samps = buff.size() / 2;
        md.end_of_burst = true;
        int ssent = (int) tx_stream->send(&buff.front(), num_tx_samps, md);
        std::cout << "Sent synch" << std::endl;
        boost::this_thread::sleep(boost::posix_time::milliseconds(1000)); // txpar->msdelay));
	counter ++;
    } while(not stop_signal_called);

    //finished
    std::cout << std::endl << "tx stream done!" << std::endl << std::endl;
    return NULL;

}

typedef boost::function<uhd::sensor_value_t (const std::string&)> get_sensor_fn_t;

bool check_locked_sensor(std::vector<std::string> sensor_names, const char* sensor_name, get_sensor_fn_t get_sensor_fn, double setup_time){
    if (std::find(sensor_names.begin(), sensor_names.end(), sensor_name) == sensor_names.end())
        return false;

    boost::system_time start = boost::get_system_time();
    boost::system_time first_lock_time;

    std::cout << boost::format("Waiting for \"%s\": ") % sensor_name;
    std::cout.flush();

    while (true) {
        if ((not first_lock_time.is_not_a_date_time()) and
                (boost::get_system_time() > (first_lock_time + boost::posix_time::seconds(setup_time))))
        {
            std::cout << " locked." << std::endl;
            break;
        }
        if (get_sensor_fn(sensor_name).to_bool()){
            if (first_lock_time.is_not_a_date_time())
                first_lock_time = boost::get_system_time();
            std::cout << "+";
            std::cout.flush();
        }
        else {
            first_lock_time = boost::system_time();	//reset to 'not a date time'

            if (boost::get_system_time() > (start + boost::posix_time::seconds(setup_time))){
                std::cout << std::endl;
                throw std::runtime_error(str(boost::format("timed out waiting for consecutive locks on sensor \"%s\"") % sensor_name));
            }
            std::cout << "_";
            std::cout.flush();
        }
        boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    std::cout << std::endl;
    return true;
}

bool check_locked_rx(uhd::usrp::multi_usrp::sptr usrp)
{
    double setup_time = 1;
    check_locked_sensor(usrp->get_rx_sensor_names(0),
			"lo_locked",
			boost::bind(&uhd::usrp::multi_usrp::get_rx_sensor, usrp, _1, 0),
			setup_time);

    return true;
}

bool check_locked_tx(uhd::usrp::multi_usrp::sptr usrp)
{
    //Check Ref and LO Lock detect
    std::vector<std::string> sensor_names;
    sensor_names = usrp->get_tx_sensor_names(0);
    if (std::find(sensor_names.begin(), sensor_names.end(), "lo_locked") != sensor_names.end()) {
        uhd::sensor_value_t lo_locked = usrp->get_tx_sensor("lo_locked",0);
        std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string() << std::endl;
        UHD_ASSERT_THROW(lo_locked.to_bool());
    }

    return true;
}

// set up a single device specified by an optional args parameter and rx/tx gains
// setup rx on A:A and use RX2 for receiving and TX/RX for transmitting (if requested)
uhd::usrp::multi_usrp::sptr usrp_setup(std::string chainname,
				       std::string args,
				       int rxspectrum_part,
				       bool setup_tx,
				       double rxgain,
				       double txgain)
{
    // create a usrp device
    std::cout << boost::format("Setting up chain %s ") % chainname << std::endl;
    std::cout << boost::format("Creating the usrp device with: %s...") % args << std::endl;
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);
    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

    // lock mboard clocks
    std::string ref = "internal";
    usrp->set_clock_source(ref);

    // set the sample rate
    double rxrate = 44000000;
    std::cout <<
	boost::format("Setting RX Rate for %s: %f Msps...") % chainname % (rxrate / 1e6) <<
	std::endl;
    usrp->set_rx_rate(rxrate);
    std::cout <<
	boost::format("Actual RX Rate for %s: %f Msps...") % chainname % (usrp->get_rx_rate() / 1e6) <<
	std::endl;

    // set the center frequency
#ifdef STANDARD
    double rxfreq = 2420000000.0;
    if(rxspectrum_part > 0)
	rxfreq = 2460000000.0;
#else
    double rxfreq = 2421500000.0;
    if(rxspectrum_part > 0)
	rxfreq = 2460500000.0;
#endif // STANDARD
    std::cout <<
	boost::format("Setting RX Freq for %s: %f MHz...") % chainname % (rxfreq / 1e6) <<
	std::endl;
    uhd::tune_request_t tune_request(rxfreq);
    usrp->set_rx_freq(tune_request);
    std::cout <<
	boost::format("Actual RX Freq for %s: %f MHz...") % chainname % (usrp->get_rx_freq() / 1e6) <<
	std::endl;

    //set the rf gain
    std::cout <<
	boost::format("Setting RX Gain for %s: %f dB...") % chainname % rxgain <<
	std::endl;
    usrp->set_rx_gain(rxgain);
    std::cout << boost::format("Actual RX Gain for %s: %f dB...") % chainname % usrp->get_rx_gain() <<
	std::endl;

    //set the antenna
    std::string rxant = "RX2";
    usrp->set_rx_antenna(rxant);

    if(setup_tx == false) return usrp;

    //set the sample rate for txing
    double txrate = 2000000;
    std::cout <<
	boost::format("Setting TX Rate for %s: %f Msps...") % chainname % (txrate / 1e6) <<
	std::endl;
    usrp->set_tx_rate(txrate);
    std::cout <<
	boost::format("Actual TX Rate for %s: %f Msps...") % chainname % (usrp->get_tx_rate() / 1e6) <<
	std::endl;

    //set the center frequency
    double txfreq = 2440000000.0;
    std::cout <<
	boost::format("Setting TX Freq for %s: %f MHz...") % chainname % (txfreq / 1e6) <<
	std::endl;
    tune_request = uhd::tune_request_t(txfreq);
    usrp->set_tx_freq(tune_request);
    std::cout <<
	boost::format("Actual TX Freq for %s: %f MHz...") % chainname % (usrp->get_tx_freq() / 1e6) <<
	std::endl;

    //set the rf gain
    std::cout <<
	boost::format("Setting TX Gain for %s: %f dB...") % chainname % txgain <<
	std::endl;
    usrp->set_tx_gain(txgain);
    std::cout << boost::format("Actual TX Gain for %s: %f dB...") % chainname % usrp->get_tx_gain() <<
	std::endl;

    //set the antenna
    std::string txant = "TX/RX";
    usrp->set_tx_antenna(txant);

    return usrp;
}


int UHD_SAFE_MAIN(int argc, char *argv[])
{
    uhd::set_thread_priority_safe();

    //variables to be set by po
    std::string loargs, hiargs, lorxfile, hirxfile;
    size_t samples_per_buffer = SAMPLE_LENGTH; 
    double lorxgain, hirxgain;
    double txgain;
    bool continue_on_bad_packet = true;
    double rx_total_time;
    bool have_lo = false;
    bool have_hi = false;
    bool runsynch = true;

    //setup the program options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
	    ("runlo", "operate lo band")
        ("loargs", po::value<std::string>(&loargs)->default_value(""), "multi uhd device address args for lo band")
        ("lorxgain", po::value<double>(&lorxgain)->default_value(50), "rx gain for the lo band RF chain")
        ("lorxfile", po::value<std::string>(&lorxfile)->default_value("usrp_samples_lo.dat"), "name of the file to write lo band binary samples to")
	    ("runhi", "operate hi band")
        ("hiargs", po::value<std::string>(&hiargs)->default_value(""), "multi uhd device address args for hi band")
        ("hirxgain", po::value<double>(&hirxgain)->default_value(50), "rx gain for the hi band RF chain")
        ("hirxfile", po::value<std::string>(&hirxfile)->default_value("usrp_samples_hi.dat"), "name of the file to write hi band binary samples to")
        ("duration", po::value<double>(&rx_total_time)->default_value(0), "total number of seconds to receive")
	("txgain", po::value<double>(&txgain)->default_value(50), "tx gain for the RF chain")
	("nosynch", "do not transmit synch frame");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("nosynch")) {
        runsynch = false;
    }

    //print the help message
    if (vm.count("help")) {
        std::cout << boost::format("Dual band capable receiver to fileUHD RX samples to file. %s") % desc << std::endl;
        std::cout
            << std::endl
            << "This application streams data from one or two USRP devices to a file.\n"
            << std::endl;
        return ~0;
    }

    if (vm.count("runlo")) have_lo = true;

    if (vm.count("runhi")) have_hi = true;

    if(have_lo == false && have_hi == false) {
	std::cout << "No band set, exiting..." << std::endl;
        return 0;
    }

    uhd::usrp::multi_usrp::sptr usrplo;
    if(have_lo == true)
        usrplo = usrp_setup("chainlo", loargs, 0, runsynch, lorxgain, txgain);

    uhd::usrp::multi_usrp::sptr usrphi;
    if(have_hi == true)
        usrphi = usrp_setup("chainhi", hiargs, 1, false, hirxgain, 0);

    double setup_time = 1;
    boost::this_thread::sleep(boost::posix_time::seconds(setup_time)); //allow for some setup time

    //check Ref and LO Lock detect
    if(have_lo == true) {
        check_locked_rx(usrplo);
        if (runsynch) {
            check_locked_tx(usrplo);
        }
    }
    if(have_hi == true)
	check_locked_rx(usrphi);

    boost::this_thread::sleep(boost::posix_time::seconds(1)); //allow for some setup time

    //set sigint if user wants to interrupt
    std::signal(SIGINT, &sig_int_handler);
    std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

    pthread_t rxlothread, txlothread, storelothread;
    struct rxthreadpars rxlopar;
    struct txthreadpars txlopar;
    if(have_lo == true) {
        rxlopar.usrp = &usrplo;
        rxlopar.file = lorxfile;
        rxlopar.samps_per_buff = samples_per_buffer;
        rxlopar.total_time = rx_total_time;
        rxlopar.continue_on_bad_packet = continue_on_bad_packet;

        if(posix_memalign((void **) &rxlopar.bankA, ALIGNMENT, HW_SAMPLE_SIZE * samples_per_buffer) != 0 ||
           posix_memalign((void **) &rxlopar.bankB, ALIGNMENT, HW_SAMPLE_SIZE * samples_per_buffer) != 0) {
	    std::cerr << "Cannot allocate memory." << std::endl;
            return 0;
	}
	rxlopar.semA[0] = 0; rxlopar.semA[1] = 0;
	rxlopar.semB[0] = 1; rxlopar.semB[1] = 1;

        txlopar.usrp = &usrplo;
        txlopar.samps_per_buff = samples_per_buffer;

        pthread_create(&rxlothread, NULL, rx_thread_routine, (void *) &rxlopar);
        pthread_create(&storelothread, NULL, store_thread_routine, (void *) &rxlopar);

        if (runsynch) {
            pthread_create(&txlothread, NULL, tx_thread_routine, (void *) &txlopar);
        }
    }

    pthread_t rxhithread, storehithread;
    struct rxthreadpars rxhipar;
    if(have_hi == true) {
	rxhipar.usrp = &usrphi;
	rxhipar.file = hirxfile;
	rxhipar.samps_per_buff = samples_per_buffer;
	rxhipar.total_time = rx_total_time;
	rxhipar.continue_on_bad_packet = continue_on_bad_packet;

        if(posix_memalign((void **) &rxhipar.bankA, ALIGNMENT, HW_SAMPLE_SIZE * samples_per_buffer) != 0 ||
           posix_memalign((void **) &rxhipar.bankB, ALIGNMENT, HW_SAMPLE_SIZE * samples_per_buffer) != 0) {
	    std::cerr << "Cannot allocate memory." << std::endl;
            return 0;
	}
        rxhipar.semA[0] = 0; rxhipar.semA[1] = 0;
        rxhipar.semB[0] = 1; rxhipar.semB[1] = 1;

	pthread_create(&rxhithread, NULL, rx_thread_routine, (void *) &rxhipar);
	pthread_create(&storehithread, NULL, store_thread_routine, (void *) &rxhipar);
    }

    if(have_lo == true) {
        pthread_join(rxlothread, NULL);
        if (runsynch) {
            pthread_join(txlothread, NULL);
        }
	pthread_join(storelothread, NULL);
    }
    if(have_hi == true) {
        pthread_join(rxhithread, NULL);
	pthread_join(storehithread, NULL);
    }

    return EXIT_SUCCESS;
}
