#!/usr/bin/perl -w
use strict;

#generate random 8 digits string

my $len=shift @ARGV;
my $stringlen=8;
if(defined $len && length($len)>0) {
	$stringlen=$len;
}

my $string;

for(my $num=0;$num<$stringlen;$num++) {
	$string.=get_rand();
}

print STDERR $string,"\n";


sub get_rand {
	my @chars=("A".."Z","a".."z",0..9,"\$","\%","\&","-","*","\@");
	my $index   = rand @chars;
    my $element = $chars[$index];
	return $element;
}	

