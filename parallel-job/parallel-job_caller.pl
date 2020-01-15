#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use List::Util qw(min);


#CutAdapt+FASTQC+RSEM+STAR
########
#Updates
########

#0.11 change procs to ppn, procs is still usable but hidden


########
#Prerequisites
########

#none


########
#Interface
########


my $version="0.11";


my $usage="

parallel-job
version: $version
Usage: sbptools parallel-job -i yourscripts.sh -n yourscripts -o jobsubmissionfolder -t 5 --nodes 1 --procs 4 -m 10gb -r
Or simpley type: sbptools parallel-job -i yourscripts.sh -r

Description: In Firefly, BSR HPC cluster, use multiple controled qsub sessions for paralleling.

Parameters:

    --in|-i           shell script file(s) with one command per line

    --name|-n         Prefix name of the task. Default as your script name
	
    Output files
    --wo|-o           working output directory. Default as folder of your first input script
    --eo              SGE error message output directory. Default as folder of your input script
    --oo              SGE output message output directory. Default as folder of your input script

    Control the tasks submited to the cluster
    --task|-t         No. of tasks submitted to cluster for each script file[10]

    For each task, there are two ways of specifying the computing resource,
      but you can't mix --nodes and --ncpus together.
	A) by specifying number of nodes and process
    --nodes           No. of nodes for each task
    #--procs           No. of processors for each task
    --ppn             No. of processes for each task	
	B) by specifying the total number of cpus
    --ncpus           No. of cpus for each task for tasks can't use multiple nodes

    --mem|-m          Memory usage for each process, e.g. 10gb


    --tandem|--td     Only used when multiple task files in input
                           Each task file needs to be successfully ran before the next one can be started


    --runmode|-r      
    --env|-e          Use your own runing envir, e.g. to source ~/.bashrc
	
	
";

#In Falco, BSR Linux server, use screen+parallel to parallelly running jobs in background.


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts


########
#Parameters
########

my $infiles;
my $nodes;
my $ppn;
my $ncpus;
my $mem;
my $task=10;
my $procs;
my $runmode=0;
my $verbose=1;
#my $queue="";
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
	"ncpus=s" => \$ncpus,
	"nodes=s" => \$nodes,	
	"procs=s" => \$procs,
	"ppn=s" => \$ppn,	
	"mem|m=s" => \$mem,
	"task|t=s" => \$task,
	#"queue|q=s" => \$queue,
	"runmode"=> \$runmode,
	#"hold|d=s" => \$hold, #can be implemented if needed to #PBS -W 
	"tandem"=> \$tandem,
	"env|e" => \$env,
	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);


my $logfile="parallel-job_run.log";
my $submitscriptfile="parallel-job_submit.sh";



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


########
#Check parameters
########

#either --nodes and --procs, or --ncpus, but not both
#this step is handled by PBS


######
#Process input file
######

#read infiles and folders
my @infiles=split(",",$infiles);
my @infile_abspaths=map {abs_path($_)} @infiles;
my @infolders=map abs_path_dir($_), @infile_abspaths;
my @qjnames=map file_short_name($_),@infiles;


my @names;

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

my $firstfolder; #for overall submit scripit

my @previous_jobs;
my @submit_scripts;

for(my $filenum=0;$filenum<@infiles;$filenum++) {
	
	#for each file in the file list

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
	if(defined $wo ) {
		unless(-e $wo) {
			mkdir($wo);
		}
		$current_wo=abs_path($wo);
		
		$firstfolder=$wo;
	}
	else {
		$current_wo=$infolder;
		
		unless(defined $firstfolder && length($firstfolder)>0) {
			$firstfolder=$current_wo;
		}
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
	
	$current_so="$current_wo/$qjname\_submit"; #script output
	
	
	my $submit_script="$current_wo/$qjname\_submit.sh";
	open(OUT,">$submit_script") || die "Error writing $submit_script. $!";
	
	push @submit_scripts,$submit_script;
	
	unless(-e $current_so) {
		system("mkdir $current_so");
	}
	
	#params for this job
	my $job_params={
		#so,wo,eo,oo,mem,que,procs,env
		wo=>$current_wo,
		eo=>$current_eo,
		oo=>$current_oo,
		so=>$current_so,
		mem=>$mem,
		#que=>$queue,
		ncpus=>$ncpus,
		nodes=>$nodes,
		ppn=>$ppn,		
		procs=>$procs,		
		env=>$env
	};
	
	#open log file #add up by multipe files
	open(LOG,">>$current_wo/$logfile") || die "Error writing $current_wo/$logfile. $!";
	
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
	open(IN,$infile) || die "Error reading $infile. $!";
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
	
	
	#generate qj files
	for(my $num=0;$num<@split_command_lines;$num++) {
		#qj file name
		my $commandline_script="$qjname\_".form_num($num+1,scalar(@split_command_lines)).".sh";
		
		write_qj_task($split_command_lines[$num],$commandline_script,$job_params);
		
		#qsub cmd
		unless($tandem) {
			print OUT "qsub $current_so/$commandline_script\n";
		}
		else {
			#write qj task
			#tandem submission controlled by -W depend=afterok:job1:job2:job3
			
			my $qjnamevar="$qjname\_".form_num($num+1,scalar(@split_command_lines));
			
			#remove special char from qjname variables
			$qjnamevar=~s/\W/_/g;
			
			if(@previous_jobs) {
				print OUT $qjnamevar,"=\$(qsub -W depend=afterok:",join(":",map {"\$".$_} @previous_jobs)," $current_so/$commandline_script)\n";
			}
			else {
				print OUT $qjnamevar,"=\$(qsub $current_so/$commandline_script)\n";				
				
			}
			
			push @current_jobs,$qjnamevar;
			
		}
	}
	
	push @previous_jobs,@current_jobs;
}

close OUT;


####
#Run job submission script
####

#merge all files for --tandem option
system("cat ".join(" ",@submit_scripts)." > $firstfolder/$submitscriptfile");

#my $runcommand=join(";",map {"sh ".$_} @submit_scripts);
my $runcommand="sh $firstfolder/$submitscriptfile";

#current run mode is, -r turn on or off
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
	#so,wo,eo,oo,mem,que,procs,env
	
	#return pbs file name
	
	#split every line
	open(OUT2,">".$params->{"so"}."/$commandline_script") || die "Error writing ".$params->{"so"}."/$commandline_script";
	print OUT2 "#!/bin/sh

##specify which shell is used
#PBS -S /bin/bash
##job name
#PBS -N $commandline_script\n";
	
#queue function skip for now. Only default queue implemented in Firefly
#print OUT2 "## Queue name\n#PBS -q ",$params->{"que"},"\n";


#computing resource setting
if(defined $params->{"nodes"}) {
	if(defined $params->{"ppn"}) {
		print OUT2 "## Parallel environment to defind # of cores\n#PBS -l nodes=",$params->{"nodes"},":ppn=",$params->{"ppn"},"\n";
	}
	else {
		print OUT2 "## Parallel environment to defind # of cores\n#PBS -l nodes=",$params->{"nodes"},"\n";
	}
}

if(defined $params->{"procs"}) {
	print OUT2 "## Parallel environment to defind # of procs\n#PBS -l procs=",$params->{"procs"},"\n";
}

if(defined $params->{"ncpus"}) {
	print OUT2 "## Parallel environment to defind # of cpus\n#PBS -l ncpus=",$params->{"ncpus"},"\n";
}

if(defined $params->{"mem"}) {
	print OUT2 "#PBS -l mem=",$params->{"mem"},"\n"; 
}

#time limit
print OUT2
"## Set time limit
#PBS -l walltime=1600:00:00\n";

#work directory 
#print OUT2 "## Set working directory\n#PBS -wd ",$params->{"wo"},"\n";

print OUT2 "## Stdout and stderr log files\n";
print OUT2 "#PBS -o ",$params->{"oo"},"/$commandline_script.out.txt\n";
print OUT2 "#PBS -e ",$params->{"eo"},"/$commandline_script.err.txt\n";

print OUT2 "##########################################\n";

	if($params->{"env"}) {
		print OUT2 "source ~/.bashrc\n";
	}
	
	#set work directory
	print OUT2 "cd ",$params->{"wo"},"\n\n";

	print OUT2 $command_line,"\n";
	
	print OUT2 "##########################################\n";
	close OUT2;
	
	return $params->{"so"}."/$commandline_script";
	
}


sub form_num {
	my ($num,$maxnum)=@_;
	my $maxnum_len=length($maxnum)>1?length($maxnum):2;
	
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
		$filename_short=$1;
	}
	else {
		$filename_short=$filename;
	}
	
	return $filename_short;
}


