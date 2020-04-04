//
// Copyright 2020 Marco Cominelli
// Copyright 2017 Francesco Gringoli
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

#include "blegenerator.hpp"

namespace po = boost::program_options;
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

struct txthreadpars {
    uhd::usrp::multi_usrp::sptr *usrp;
    size_t samps_per_buff;
    double msdelay;
    int channel;
};

void *tx_thread_routine(void *pars)
{
    struct txthreadpars *txpar = (struct txthreadpars *) pars;
    uhd::usrp::multi_usrp::sptr usrp = *(txpar->usrp);
    const std::string cpu_format = "sc8";
    size_t samps_per_buff = txpar->samps_per_buff;
    int channel = txpar->channel;

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
        std::vector<int8_t> buff = generatesamples(channel, counter, ZERO_PAD_BYTES);
        size_t num_tx_samps = buff.size() / 2;
        md.end_of_burst = true;
        int ssent = (int) tx_stream->send(&buff.front(), num_tx_samps, md);
        std::cout << "Sent synch" << std::endl;
        boost::this_thread::sleep(boost::posix_time::milliseconds(1000));
    counter ++;
    } while(not stop_signal_called);

    //finished
    std::cout << std::endl << "tx stream done!" << std::endl << std::endl;
    return NULL;

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

double ble_chan2freq(int channel)
{
    double freq = 2402;
        if(channel < 11)
                freq = 2404 + channel * 2;
        else if(channel < 37)
                freq = 2428 + (channel - 11) * 2;
        else if(channel == 37)
                freq = 2402;
        else if(channel == 38)
                freq = 2426;
        else if(channel == 39)
                freq = 2480;
        else {
                fprintf(stderr, "Invalid channel, defaulting to 37\n");
                freq = 2402;
        }

    return freq * 1e6;
}

// set up a single device for transmitting at the requested channel
uhd::usrp::multi_usrp::sptr usrp_setup(std::string chainname,
                       std::string args,
                       double txgain,
                       int blechannel)
{
    double txfreq = ble_chan2freq(blechannel);

    // create a usrp device
    std::cout << boost::format("Setting up chain %s ") % chainname << std::endl;
    std::cout << boost::format("Creating the usrp device with: %s...") % args << std::endl;
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);
    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

    // lock mboard clocks
    std::string ref = "internal";
    usrp->set_clock_source(ref);

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
    std::cout <<
    boost::format("Setting TX Freq for %s: %f MHz...") % chainname % (txfreq / 1e6) <<
    std::endl;
    uhd::tune_request_t tune_request = uhd::tune_request_t(txfreq);

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
    std::string args;
    double txgain;
    int channel;
    size_t samples_per_buffer = 10000;

    //setup the program options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&args)->default_value(""), "multi uhd device address args")
        ("channel", po::value<int>(&channel)->default_value(37), "ble channel")
        ("txgain", po::value<double>(&txgain)->default_value(50), "tx gain for the RF chain");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    //print the help message
    if (vm.count("help")) {
        std::cout << boost::format("ble advertiser. %s") % desc << std::endl;
        std::cout << std::endl;
        return ~0;
    }

    uhd::usrp::multi_usrp::sptr usrp;
    usrp = usrp_setup("chain", args, txgain, channel);

    double setup_time = 1;
    boost::this_thread::sleep(boost::posix_time::seconds(setup_time)); //allow for some setup time

    //check Ref and LO Lock detect
    check_locked_tx(usrp);

    boost::this_thread::sleep(boost::posix_time::seconds(1)); //allow for some setup time

    //set sigint if user wants to interrupt
    std::signal(SIGINT, &sig_int_handler);
    std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

    pthread_t txthread;
    struct txthreadpars txpar;
    txpar.usrp = &usrp;
    txpar.samps_per_buff = samples_per_buffer;
    txpar.channel = channel;
    pthread_create(&txthread, NULL, tx_thread_routine, (void *) &txpar);
    pthread_join(txthread, NULL);

    return EXIT_SUCCESS;
}

