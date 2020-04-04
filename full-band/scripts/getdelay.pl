#!/usr/bin/perl

if(scalar @ARGV < 2) {
  print "Missing parameters\n";
  print "Usage: getdelay.pl /path/to/synchtime.log starttime_us\n";
  exit(1);
}

my $filename = $ARGV[0];
my $starttime = $ARGV[1];
open(my $fh, $filename)
  or die "Could not open file '$filename' $!";
my $initialskip = 0;
my $initialdelay = 0;
while(<$fh>) {
  if($_ =~ /([\d\.]+)\s+([\-\d]+)/) {
    my $skipus = $1 * 1000000;
    my $delayus = $2;
    if($starttime <= $skipus) {
      if($initialdelay == 0) {
        print "$initialskip $delayus\n";
      } else {
        print "$initialskip $initialdelay\n";
      }
      exit(0);
    }
    $initialskip = $skipus;
    $initialdelay = $delayus;
  }
}

print "$initialskip $initialdelay\n";
