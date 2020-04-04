# btsniffer
## Tool for discovering and de-anonymizing Classic Bluetooth connections

The programs in this repository implement two experimental wide-band Bluetooth sniffers based on Software-Defined Radios (SDR). The goal is to evaluate the performance of these systems when used as tools to discover and de-anonymize Classic Bluetooth connections (i.e. discover the UAP of the master device). The software has been tested on Ubuntu Linux v16.04, using as SDR platforms the Ettus USRP B210.

## What's in this repository

This repository contains the code for running two different systems. Both of them run on a host computer connected to an SDR frontend.

* **full-band**: this experimental system includes a set of programs used to capture an 80-MHz trace of the 2.4-GHz ISM band and to process this trace offline. Two USRP B210 are required to run this system.

* **single-board**: this system can run in real-time on a host computer, capturing 8 or 16 Bluetooth channels concurrently in order to detect all Bluetooth sessions. Only one USRP B210 is needed to run this system.

## Setup

To compile and run the programs contained in this repository, please have a look at the README file in the corresponding subdirectories.
