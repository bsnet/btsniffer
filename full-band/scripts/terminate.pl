#!/usr/bin/perl
#
#***** RESOLVED 219F87
#RES 5.33091206 (30065) lap = 219F87, uap = 74, clks 12-51, clk_index = 0
#RES 5.33091506 (30065) lap = 219F87, uap = 156, clks 19-44, clk_index = 1
#
#
#*****************
#
if(scalar @ARGV < 1) {
  print "Missing parameters\n";
  print "Usage: decoder .. | ./terminate.pl lap\n";
  exit(1);
}

my $addr = "";
my $found = 0;
while(<STDIN>) {
  if($_ =~ /\A\*\*\*\*\* RESOLVED (.{6})/) {
    $addr = $1;
    if (lc($addr) eq lc($ARGV[0])) {
      $found = 1;
    }
  }
  elsif($found == 1 and $_ =~ /\ARES ([\d\.]+)\s+\(\d+\)\s+lap\s+\=\s+(.{6}),\s+uap\s+\=\s+(\d+)/) {
    print "$_";
  }
  # 6 224 0 0 1 0 0 0
  elsif($found == 1 and $_ =~ /(\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+) (\d+)/) {
    print "$_";
  }
  elsif($found == 1 and $_ =~ /\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/) {
    $found = 0;
  }
}
