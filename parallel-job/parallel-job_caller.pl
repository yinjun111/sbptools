#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#CutAdapt+FASTQC+RSEM+STAR


########
#Prerequisites
########

#none


########
#Interface
########


my $version="0.1";


my $usage="

parallel-job
version: $version
Usage: sbptools parallel-job

Description: In Falco, BSR Linux server, use screen+parallel to parallelly running jobs in background.
             In Firefly, BSR HPC cluster, use multiple controled qsub sessions for paralleling.

Parameters:

    --in|-i           task file(s), shell script with one command per line

    Output files
    --wo|-o           working output directory. Default as folder of input script 
    --eo              SGE error message output directory. Default as folder of your input script
    --oo              SGE output message. Default as input script

    Control the tasks
    --ppn|-p          No. of precesses for each task
    --mem|-m          Memory usage, e.g. 10G
    --queue|-q        Queue in cluster
    --task|-t         No. of tasks or paralleling process to cluster
    --hold|-d         JOb dependency, e.g. -hold_jid in qsub
    --tandem|--td     Only used when multiple taskes files in input
                           Each task needs to be finished before the next one is ran
	--name|-n         Prefix name of task
    --runmode|-r      
    --env|-e          Use your own runing envir, e.g. to source ~/.bashrc
	
	
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

my $infile;
my $ppn=4;
my $mem;
my $task=20;
my $runmode=0;
my $verbose=1;
my $queue="";
my $env=0;
my $hold;
my $tandem=0;
my $name;

#output
my $wo;
my $eo;
my $oo;

GetOptions(
	"in|i=s" => \$infiles,
	"name|n=s" => \$name,
	"wo|o=s" => \$wo,
	"eo=s" => \$eo,
	"oo=s" => \$oo,
	"ppn|p=s" => \$ppn,
	"mem|m=s" => \$mem,
	"task|t=s" => \$task,
	"queue|q=s" => \$queue,
	"runmode"=> \$runmode,
	"hold|d=s" => \$hold,
	"tandem|d"=> \$tandem,
	"env|e" => \$env,
	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);

#my $logfile="parallel-job_run.log";
my $logfile=$outfile; #need to check log file ...
$logfile=~s/\.txt/_parallel-job_run.log/;

#write log file
#open(LOG, ">$logfile") || die "Error write $logfile. $!";
#my $now=current_time();
#print LOG "perl $0 $params\n\n";
#print LOG "Start time: $now\n\n";
#print LOG "Current version: $version\n\n";
#print LOG "\n";


########
#Welcome message
########

#welcome message

my @userreals=getpwuid($<);
my $user=$userreals[0];
my @userattrs=getpwnam($user);
my @groupattrs=getgrgid($userattrs[3]);


print STDERR "\nWelcome $userattrs[6]($user) from $groupattrs[0] to Firefly!\n";
#print LOG "\nWelcome $userattrs[6]($user) from $groupattrs[0] to Firefly!\n";


#queue selection?
#what queues do we have?


######
#Process input file
######

my @infiles=split(",",$infile);
my @infile_abspaths=map abs_path_dir($_), @infiles;
my @infolders=map abs_path_dir($_), @infile_abspaths;
my @qjnames=map file_short_name($_),@infiles;

if(defined $name) {
	@names=split(",",$name); #name of jobs
}


if(defined $name && @names!=@infiles) {
	print STDERR "\nERROR: ",scalar(@infiles), " task files found in -i ",join(",",@infiles)," different from ",scalar(@names)," names in -n ",join(",",@names), "\n\n";
	exit;
}



####
#Generate scripts
####

#non-tandam may use more tasks,e.g. task x files


my @previous_jobs;
my @submit_scripts;

for(my $filenum=0;$filenum<@infiles;$filenum++) {
	my $infile=$infiles[$filenum];
	my $infolder=$infolders[$filenum];
	my $qjname;
	
	if(defined $name) {
		#use predefined job name
		$qjname=$names[$filenum];
	}
	else {
		$qjname=$qjnames[$filenum];
	}
	
	my @current_jobs; #record submitted job names for current file
	my ($current_wo,$current_eo,$current_oo,$current_so);
	
	#output folder
	if(defined $wo) {
		unless(-e $wo) {
			mkdir($wo);
		}
		$current_wo=abs_path($wo);
	}
	else {
		$current_wo=$infolder;
	}
	
	#
	if(defined $eo) {
		unless(-e $eo) {
			mkdir($eo);
		}
		$current_eo=abs_path($eo);
	}
	else {
		$current_eo="$current_wo/$qjname\_submit";
	}
	
	#
	if(defined $oo) {
		unless(-e $oo) {
			mkdir($oo);
		}
		$current_oo=abs_path($oo);
	}
	else {
		$current_oo="$current_wo/$qjname\_submit";
	}	
	
	my $submit_script="$current_wo/$qjname\_submit.sh";
	$current_so="$current_wo/$qjname\_submit"; #script output
	push @submit_scripts,$submit_script;
	
	unless(-e $current_so) {
		system("mkdir $current_so");
	}
	
	#params for this job
	my $job_param={
		#so,wo,eo,oo,mem,que,ppn,env
		wo=>$current_wo,
		eo=>$current_eo,
		oo=>$current_oo,
		so=>$current_so,
		mem=>$mem,
		que=>$queue,
		ppn=>$ppn,
		env=>$env
	};
	
	#open log file #add up by multipe files
	open(LOG,">>$current_wo/$logfile") || die $!;
	
	print LOG "perl $0 $params\n\n";
	print LOG "parallel-job version $version\n\n";
	my $starttime=localtime();
	print LOG "Program started at $starttime.\n\n";
	
	print LOG "\nWelcome $userattrs[6]($user) from $groupattrs[0] to Firefly!\n";
	#print LOG "Selecting $queue \n";
	
	print STDERR "\n\n" if $verbose;
	
	#######
	#Read input file
	#######
	
	#output separate scripts
	open(IN,$infile) || die $!;
	my @command_lines;
	
	while(<IN>) {
		tr/\r\n//d;
		next if $_=~/^#/; #skip comments
		next if length($_)==0; #skip empty lines
		push @command_lines,$_;
	}
	close IN;
	
	print STDERR scalar(@command_lines)," command lines found in $infile.\n\n" if $verbose;
	print LOG scalar(@command_lines)," command lines found in $infile.\n\n";
	
	#split command lines by # of tasks
	my @split_command_lines=split_jobs(\@command_lines,min($task,scalar(@command_lines)));
	
	print STDERR "Split command lines into ",min($task,scalar(@command_lines))," tasks.\n\n" if $verbose;
	
	open(OUT,">$submit_script") || die $!;
	#generate qj files
	for(my $num=0;$num<@split_command_lines;$num++) {
		#qj file name
		my $commandline_script="$qjname\_".form_num($num+1,scalar(@split_command_lines)).".sh";
		#qsub cmd
		print OUT "qsub $current_so/$commandline_script\n";
		
		#write qj task
		#tandem submission controlled by hold_ids
		
		if(($tandem && @previous_jobs) || (defined $hold && length($hold)>1)) {
			#with either multiple submission, or defined hold ids
			my @holdids;
			if($tandem && @previous_jobs) {
				push @hold_ids,@previous_jobs;
			}
			if(defined $hold && length($hold)>1) {
				push @holdids, split(",",$hold);
			}
			
			$job_params->{"hold_ids"}=join(",",@holdids);
		}
		
		#function to write the script
		write_qj_task($split_command_lines[$num],$commandline_script,$job_params);
		
		push @current_jobs,$commandline_script;
	}
	close OUT;
	
	push @previous_jobs,@current_jobs;
}



####
#Run job submission script
####

#current run mode is, -r turn on or off

my $runcommand=join(";",map {"sh ".$_} @submit_scripts);

if($runmode) {
	print STDERR "Submitting server jobs, running $runcommand.\n\n" if $verbose;
	print LOG "Submitting server jobs, running $runcommand.\n\n";
	system($runcommand);
	print STDERR "Done.\n\n" if $verbose;
	print LOG "Done.\n\n";
}
else {
	print STDERR "In server, run $runcommand.\n\n" if $verbose;
	print LOG "In server, run $runcommand.\n\n";
}

close LOG;



########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}

sub split_jobs {
	#split the task file into several files
	my ($command,$task)=@_;
	my @splitcommand;
	
	for(my $num=0;$num<@{$command};$num++) {
		if(defined $splitcommand[$num % $task]) {
			$splitcommand[$num%$task]=$splitcommand[$num%$task]."\n".$command->[$num];
		}
		else {
			$splitcommand[$num%$task]=$command->[$num];
		}
	}
	
	return @splitcommand;
}

sub write_qj_task {
	my ($command_line,$commandline_script,$params)=@_;
	
	#global var
	#so,wo,eo,oo,mem,que,ppn,env
	
	#split every line
	open(OUT2,">".$params->{"so"}."/$commandline_script") || die $!;
	print OUT2 "#!/bin/sh

##specify which shell is used
#\$ -S /bin/bash
##job name
#\$ -N $commandline_script\n";

	print OUT2 "## Queue name\n#\$ -q ",$params->{"que"},"\n";
	print OUT2 "## Parallel environment to defind # of cores\n#\$ -pe smp ",$params->{"ppn"},"\n";
	
	if(defined $params->{"mem"}) {
		print OUT2 "#\$ -l h_vmem=",$params->{"mem"},"\n";
	}

print OUT2
"#\$ -R y
## Set time limit
#\$ -l h_rt=1600:00:00\n";
print OUT2 "## Set working directory\n#\$ -wd ",$params->{"wo"},"\n";

if(defined $params->{"hold_ids"}) {
	print OUT2 "#\$ -hold_jid ",$params->{"hold_ids"},"\n";
}

print OUT2 "## Stdout and stderr log files\n";
print OUT2 "#\$ -o ",$params->{"oo"},"/$commandline_script.out.txt\n";
print OUT2 "#\$ -e ",$params->{"eo"},"/$commandline_script.err.txt\n";

print OUT2 "##########################################\n";

	if($params->{"env"}) {
		print OUT2 "source ~/.bashrc\n";
	}
	
	print OUT2 $command_line,"\n";
	
	print OUT2 "##########################################\n";
	close OUT2;
}


sub form_num {
	my ($num,$maxnum)=@_;
	my $maxnum_len=length($maxnum)>1?length($manxnum):2;
	
	return(0 x ($maxnum_len-length($num))).$num;
}

sub abs_path_dir {
	my $path=shift @_;
	my $path_dir;
	if($path=~/[^\/]+$/) {
		$path_dir=$`;
	}
	return $path_dir;
}

sub file_short_name {
	my $filename=shift @_;
	my $filename_short;
	
	if($filename=~/([^\/]+)\.\w+$/) {
		$filename_short=$1;
	}
	elsif($filename=~/([^\/]+)$/) {
		$filename_shortj=$1;
	}
	else {
		$filename_short=$filename;
	}
	
	return $filename_short;
}


