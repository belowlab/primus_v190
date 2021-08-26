#! /usr/bin/perl

use strict;

my $val_1 = 'R241364515';  # FID1
my $val_2 = 'R253585921';  # FID2

my $new_PI_HAT = `python /data100t1/home/grahame/projects/compadre/primus/primus_old/lib/perl_modules/PRIMUS/ersa_2_pihat.py --fid1 $val_1 --fid2 $val_2`;
# my $new_PI_HAT = `python /data100t1/home/grahame/projects/compadre/primus/primus_old/lib/perl_modules/PRIMUS/ersa_2_pihat.py`;

print $new_PI_HAT;
print "\n";