#perl program to update the station layout numbers in an array file
use strict;

###############################
# Main program
################################
my $debug=0;

# check command line and required environment variables
if ($ENV{'SIM'} eq "") { die("SIM environement variable must be set up. source the sim_setup.csh script"); }

if (scalar(@ARGV) < 1){
    die "Usage: $0 filename\n";
}
if ($ARGV[$#ARGV] eq '-d') {
    $debug=1;
    printf("Debugging on\n");
    pop @ARGV;
}

my $ant=0;

open(FP,$ARGV[0]) || die "Cannot open file $ARGV[0]\n";

while(<FP>) {
    my $line = $_;
    chomp $line;
    printf  "%s%03d\n",$line, $ant;
    $ant++;
}
close FP;
