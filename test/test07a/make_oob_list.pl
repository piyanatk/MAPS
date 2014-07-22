#!/usr/bin/perl

use strict;
my $ra;
my $dec;

for ($ra=-1.5; $ra<1.5; $ra+=0.1) {
  for($dec=-45.0; $dec<-7.5; $dec+=2.0) {
    printf("%.2f %.2f 1 0 0 0\n",$ra,$dec);
  }
}
