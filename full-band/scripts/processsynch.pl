#!/usr/bin/perl

my $ts = 0;

# depending on which of LO/HI comes first we have two situations:
# 1)
# Modulation = 1Mb/s, Chan 17, Carrier 2460.50, Time = 364236us, Q9 = 1, AA = 8E89BED6
# 0000: 0x40 0x14 0x00 0x11 0x22 0x33 0x44 0x55 0x00 0x00 0x0B 0x35 0x00 0x00 0x00 0x00 @..."3DU...5....
# 0010: 0x00 0x00 0x00 0x00 0x00 0x00 0x93 0xA7 0xC7                                    .........
#
# synch packet
# delay = 128
# CALIBRATED!
#
#### or
# 2)
# Modulation = 1Mb/s, Chan 17, Carrier 2421.50, Time = 566322us, Q9 = 1, AA = 8E89BED6
# 0000: 0x40 0x14 0x00 0x11 0x22 0x33 0x44 0x55 0x00 0x00 0x03 0x48 0x00 0x00 0x00 0x00 @..."3DU...H....
# 0010: 0x00 0x00 0x00 0x00 0x00 0x00 0xE7 0x40 0x0A                                    .......@.
#
# adding to CARRIER_LO, channel = 17, ts = 0.566322
# synch packet
# delay = 0
# CALIBRATED!


while(<STDIN>) {
  # if($_ =~ /, ts = ([\d\.]+)/) {
  if($_ =~ /Modulation \= .*Time \= (\d+)us/) {
    $ts = $1 / 1000000.0;
  } elsif($_ =~ /\Adelay = ([-\d]+)/ and $ts > 0) {
    print "$ts $1\n";
  }
}
