#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Email::MIME;
use Email::Sender::Simple qw(sendmail); #get some bug warning here

########
#Prerequisites
########


########
#Interface
########


my $version="0.1";


my $usage="

add-user_omicsoft
version: $version
Usage: add-user_omicsoft [parameters]

Description: Add user, generate random password, then send email


Parameters:

    --user              Full name
    --email             Email address


";

	
unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts



########
#Parameters
########

my $user;
my $email;

GetOptions(
	"user=s" => \$user,
	"email=s" => \$email,
);


my $tempfolder="temp";


########
#Running
########

my $userid;

if($email=~/(\w+)\@/) {
	$userid=$1;
}

my $password;

for(my $num=0;$num<8;$num++) {
	$password.=get_rand();
}


########
#Run oscript
########


open(OUT,">$tempfolder/$userid.oscript") || die $!;

print OUT "

Begin ExecuteCommand /Namespace=Server;
Server \"tcp://falco.burnham.org:8065\" /UserID=jyin /Password=20120729;
Command ManageUsers_Add;
Options
\"
UserID=$userid
UserGroups={standard users}
Password=$password
FullName =$user
Email=$email
\";
End;

";

close OUT;

system("sudo /opt/mono-4.0.4/bin/mono /home/omicsoft/oshell/oshell.exe --runscript /opt/arrayserver/ $tempfolder/$userid.oscript $tempfolder/ /opt/mono-4.0.4/bin/mono > $tempfolder/$userid\_run.log");


########
#Send email
########

# first, create your message

my $message = Email::MIME->create(
  header_str => [
    From    => 'jyin@sbpdiscovery.org',
    To      => $email,
    Subject => 'Omicsoft Array Studio Server Account',
  ],
  attributes => {
    encoding => 'quoted-printable',
    charset  => 'ISO-8859-1',
  },
  body_str => "
Dear $user,

An Omicsoft ArrayStudio Server account has been created for you.

Your user ID is $userid. Your password is $password

Best regards,

Jun Yin
SBP Bionformatics Core

",
);

# send the message

sendmail($message);




########
#Functions
########

sub get_rand {
	my @chars=("A".."Z","a".."z",0..9,"\$","\%","\&","-","*","\@");
	my $index   = rand @chars;
    my $element = $chars[$index];
	return $element;
}	

