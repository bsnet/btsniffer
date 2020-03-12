# btsniffer
## Tool for discovering and de-anonymizing Classic Bluetooth connections

This software can be used to detect Bluetooth Classic connections and discover the UAP of the master device in the piconet. This software has been tested on Ubuntu Linux v16.04 using Software-Defined Radios (SDR) manufactured by Ettus Research, namely the USRP B210.

## What's in this repository

The programs here run on a host computer connected to SDR platforms.

* **single-board**: this system can run in real-time on a host computer, capturing 8 or 16 Bluetooth channels concurrently in order to quickly detect ongoing connections.

* **full-band**: code will be soon available 

## Setup

To compile and run a program contained in this repository, please have a look at the README file in the corresponding subdirectory.
