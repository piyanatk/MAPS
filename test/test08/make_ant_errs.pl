# perl program to make a bunch of antenna layout files with gain and phase errors.

use strict;

# Make Gaussian random deviates (2 of them) based on Box-Muller transform of the 
# uniform deviates provided by rand.
sub gaussRand() {
    my $r = 0.0;
    my $phi = rand()*2.0*pi();
    while(($r = rand()) == 0.0) {}; # get a uniform random num > 0
    my $arg = sqrt(-2.0*log($r));
    return $arg*cos($phi),$arg*sin($phi);
}

sub writeLayoutFile_dipole_gp{
    my $name = shift();
    my $gainerr = shift();
    my $phaserr = shift();
    my $r1;
    my $r2;
    # create the file
    open(FP,">".$name) || die "Cannot open file $name";
    
    #write the dipole info
    printf FP "# lfd basic short dipole over groundplane.\n# parameters are: N,E,Up,type,gain,phase,PA,height\n";
    
    ($r1,$r2) = gaussRand();
    printf FP "-1.6065 -1.6065 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "-1.6065 -0.5355 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "-1.6065  0.5355 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "-1.6065  1.6065 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "-0.5355 -1.6065 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "-0.5355 -0.5355 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "-0.5355  0.5355 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "-0.5355  1.6065 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "0.5355  -1.6065 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "0.5355  -0.5355 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "0.5355   0.5355 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "0.5355   1.6065 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "1.6065  -1.6065 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "1.6065  -0.5355 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "1.6065   0.5355 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;
    ($r1,$r2) = gaussRand();
    printf FP "1.6065   1.6065 0.0 2 %.5g\t% .5g\t0.0 0.35\n",1.0+$r1*$gainerr,$r2*$phaserr;

    #close file
    close FP;
}

###############################
# Main program
################################
my $debug=0;

# check command line and required environment variables
if ($ENV{'SIM'} eq "") { die("SIM environement variable must be set up. source the sim_setup.csh script"); }

if (scalar(@ARGV) < 2){
    die "Usage: $0 num_ant rel_gain_err abs_phase_err_degs -d (debug)\n";
}
if ($ARGV[$#ARGV] eq '-d') {
    $debug=1;
    printf("Debugging on\n");
    pop @ARGV;
}
if ($ARGV[0] < 1 || $ARGV[0] > 1024) {
    die("bad number of antennas: $ARGV[0] \n");
}

my $gainerr = $ARGV[1];
my $phaserr = $ARGV[2];

for(my $ant=0;$ant<$ARGV[0]; $ant++) {
    # make a filename
    my $filename = sprintf("%s%03d.layout","lfd_dipole_gp",$ant);
    # make a random instance 
    if($debug) { printf("Making file: %s\n",$filename); }
    writeLayoutFile_dipole_gp($filename,$gainerr,$phaserr);
}
