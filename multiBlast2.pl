#!/usr/bin/perl
use strict;
use File::Path;
our $THREAD_NUM=4;          # default thread number
our $BLAST_SINGLE_THREAD=4; # default single blast thread number
our $THREADS4SYSTEM=1;      # default threads for a system

our $MAINmEMORY=64;      # main memory size in gigabytes
our $MEMORY4ONEtHREAD=2; # required memory for a thread in Giga bytes

our $BLASTALL=`which blastall`;
our $MEGABLAST=`which megablast`;
our $MEGABLASTp=`which blastn`;
our $GZIP=`which gzip`;
our $GUNZIP=`which gunzip`;
our $CAT=`which cat`;
our $TMPDIR="/var/tmp";

our $BASENAME="tmpResult";
our $MAXrAND=10000;

&main;

sub main{
    my $pidList;
    my $catList;
    my %options;
    my $dirname;
    my @threads;
    %options=&interface();
    $dirname=&makedir();
#    print stderr "Will out $dirname\n";
    @threads=&countProcessors(\%options);
    ($pidList, $catList)=&divideQuery(\%options, $dirname, \@threads);
    &wait4result($pidList, $catList, \%options);
    &concatenate($catList, \%options);
    File::Path::rmtree([$dirname]) or die("$!:$dirname");
}

sub interface{
    my $i;
    my $key;
    my $value;
    my %options;
    my %exceptions=("-h" => 1, "-help" => 1, "-version" => 1
		    , "-show_gis" => 1, "-html" => 1, "-lcase_masking" => 1
		    , "-no_greedy" => 1, "-ungapped" => 1
		    , "-parse_deflines" => 1, "-remote" => 1, "-use_sw_tback" => 1);
    foreach ($i=0; $i<@ARGV; $i+=2){
	$key=$ARGV[$i];
	$value=$ARGV[$i+1];
	if(exists($exceptions{$key})){
	    $options{$key}="";
	    $i--;
	}else{
	    $options{$key}=$value;
	}
	if($key eq "--help"){
	    &help;
	}
    }

    if(defined($options{"--force_threads"})){
	if($options{"--force_threads"} < 1){
	    print stderr "--force_threads should be natural number\n";
	    &help;
	}else{
	    $options{"--force_threads"}*=1;
	}
	$BLAST_SINGLE_THREAD=$options{"--force_threads"};
    }
    if(defined($options{"--force_paralell"})){
	if($options{"--force_paralell"} < 1){
	    print stderr "--force_paralell should be natural number\n";
	    &help;
	}else{
	    $options{"--force_paralell"}*=1;
	}
	$THREAD_NUM=$options{"--force_paralell"};
    }
    if(defined($options{"--silent"}) && $options{"--silent"} ne "T"){
	delete($options{"--silent"});
    }
    return %options;
}

sub help{
    my $file=__FILE__;
    print stderr << "EOF";
$file
This program will divide the query and run blast paralell.
To run megablast (legacy) use "-p megablast".
If -p is not specified or set as "-p megablast+" 
megablast+ (default of blastn on blast+) will be executed.

=== Special options ===
 --help show this help
 --force_threads  [integer] fix threads to this number
 --force_paralell [integer] fix number of parallelization
EOF
    ;
    exit;
}

sub countProcessors{ #($options)
    my $options=$_[0];
    my $command;
    my $wholeThreads;
    my $threads;
    my $memory;
    my $singleThread;
    my $currentSingleThread;
    my @result;
    $command="grep 'processor' /proc/cpuinfo | wc -l";
    $wholeThreads=`$command`;
    $wholeThreads*=1;

    $command='free -g | grep "^Mem:" | sed "s/Mem:\s*//" | sed "s/\s.*//"';
    $memory=`$command`;
    $memory=int($memory / $MEMORY4ONEtHREAD + $THREADS4SYSTEM);

    if($wholeThreads > $memory){
	$wholeThreads=$memory;
    }

    for($singleThread=$BLAST_SINGLE_THREAD; $singleThread > 0; $singleThread--){
	if(defined($options->{"--force_threads"})){
	    $currentSingleThread=$options->{"--force_threads"};
	}else{
	    $currentSingleThread=$singleThread;
	}
	if(defined($options->{"--force_paralell"})){
	    $threads=$options->{"--force_paralell"};
	}else{
	    $threads=int(($wholeThreads - $THREADS4SYSTEM) / $currentSingleThread);
	}
	if($threads > 0){
	    @result=($threads, $currentSingleThread);
	    last;
	}
    }

    return @result;
}

sub makedir{ #(\%options);
    my $i;
    my $dirname;
    do{
	$i=rand($MAXrAND);
	$dirname=sprintf("%s/%s_%d", $TMPDIR, $BASENAME, $i);
    }while( -d "$dirname");
    mkdir $dirname;
    return $dirname;
}

sub divideQuery{ #(\%options, $dirname, \@threads);
    my $options=$_[0];
    my $dirname=$_[1];
    my $threads=$_[2]; # @(num of parallel, num of theads)
    my $line;
    my $outfile="";
    my $salt;
    my $i;
    my $j;
    my %runInfo;
    my $pid;
    my $maxnum=5000;
    my @catList;
    my $fileSize=0;
    my $readSize=0;
    my $targetSize=0;
    my $singleThread;
    my $file="";
    
    if(!defined($options->{"--silent"})){
	printf(stderr "Number of threads: %d\n", $threads->[0]);
    }
    if(defined($options->{"-i"})){
	$file=$options->{"-i"};
    }elsif(defined($options->{"-query"})){
	$file=$options->{"-query"};
    }
    if($file ne ""){
	if(-r $file ){
	    $fileSize= -s $file;
	    $targetSize=int($fileSize/$threads->[0]);
	    if(!defined($options->{"--silent"})){
		print stderr "Input is a file. Estimated size: $fileSize\n";
		print stderr "Target size will be $targetSize\n";
	    }
	    open FIN, "$file" or die("Could not open query file: $file\n");
	}
    }else{
	*FIN=*STDIN;
	$targetSize=700*1000;
    }

    $i=1;
    $outfile="$dirname/query_$i";
    open FOUT, ">$outfile" or die("Could not save divided query: $outfile\n");
    for($readSize=0; $line=<FIN>;){
	$readSize+=length($line);
	if($line=~/^>/){
	    if($readSize > $targetSize){
		$readSize = $readSize % $targetSize;
		close FOUT;
		if($outfile ne "$dirname/dummy"){
		    $pid=&fork2blast($outfile, $options, $threads->[1]);
		    sleep(1);
		    $runInfo{$pid}[0]=$i;
		    $runInfo{$pid}[1]="${outfile}.out";

		    if(scalar(keys %runInfo) >= $threads->[0]){
			&wait4oneResult(\%runInfo, \@catList, $options);
		    }
		}
		$i++;
		$outfile="$dirname/query_$i";
		open FOUT, ">$outfile" 
		    or die("Could not save divided query: $outfile\n");
	    }
	}
	print FOUT $line;
    }
    close FOUT;
    $pid=&fork2blast($outfile, $options, $threads->[1]);
    $runInfo{$pid}[0]=$i;
    $runInfo{$pid}[1]="${outfile}.out";
    close FIN;
    return (\%runInfo, \@catList);
}

sub fork2blast{
    my $query=$_[0];
    my $options=$_[1];
    my $singleThread=$_[2];
    my $i;
    my $pid;
    my %pidList;
    $pid=fork();
    if($pid==-1){
	&fork_error();
	exit;
    }elsif($pid==0){
	&runBlast($query, $options, $singleThread);
	exit;
    }
    return $pid;
}

sub fork_error{
    print stderr "Detected error on forking\n";
    exit;
}

sub runBlast{
    my $query=$_[0];
    my $options=$_[1];
    my $singleThread=$_[2];
    my $command;
    my $key;
    my $value;
    my $out="$query.out.gz";
    my $mode="legacy";
    if(!defined($options->{"-p"}) || $options->{"-p"} eq "megablast+"){
	$command="$MEGABLASTp -query $query -num_threads $singleThread";
	$mode="plus";
    }elsif(!defined($options->{"-p"}) || $options->{"-p"} eq "megablast"){
	$command="$MEGABLAST -i $query -a $singleThread";
    }else{
	$command=sprintf("$BLASTALL -i $query -a $singleThread -p %s"
			 , $options->{"-p"});
    }
    foreach $key (keys %{$options}){
	$value=$options->{$key};
	if($key ne "-i" && $key ne "-o" && $key ne "-a" && $key ne "-p"
	    && $key ne "-query" && $key ne "-out" && $key ne "-num_threads"
	    && $key ne "--force_threads" && $key ne "--force_paralell"
	    && $key ne "--silent"){
	    $command.=" $key $value";
	}
    }
    $command.=" | $GZIP > $out";
    if(!defined($options->{"--silent"})){
	print stderr "$command\n";
    }
    `$command`;
    exit;
}

sub wait4result{#(\%pidList, \@catList, \%options);
    my $pidList=$_[0];
    my $catList=$_[1];
    my $options=$_[2];
    my $pid;
    while(scalar(keys(%{$pidList})) > 0){
	&wait4oneResult($pidList, $catList, $options);
	sleep(1);
    }
}

sub wait4oneResult{#(\%pidList, \@catList, $options);
    my $pidList=$_[0];
    my $catList=$_[1];
    my $options=$_[2];
    my $pid;
    my $pid=wait();
    if(exists($pidList->{$pid})){
	if(!defined($options->{"--silent"})){
#	    print stderr "pid $pid $pidList->{$pid}[0] $pidList->{$pid}[1] finished\n";
	    printf(stderr "Pid %d %s %s finished. %d jobs remaining,\n"
		   , $pid, $pidList->{$pid}[0], $pidList->{$pid}[1]
		   , scalar(keys %{$pidList}));
	}
	$catList->[$pidList->{$pid}[0]]=$pidList->{$pid}[1];
	delete($pidList->{$pid});
    }else{
	print stderr "unknown pid $pid finished\n";
    }
}

sub concatenate{ #(@catList, $options);
    my $list=$_[0];
    my $options=$_[1];
    my $i;
    my $number=scalar(@{$list});
    my $file;
    my $line;
    my $outfile="";
    if(defined($options->{"-o"})){
	$outfile=$options->{"-o"};
    }elsif(defined($options->{"-out"})){
	$outfile=$options->{"-out"};
    }
    if($outfile ne ""){
	open FOUT, ">$outfile" or die("Could not save to the file: $outfile\n");
    }else{
	*FOUT=*STDOUT;
    }
    for($i=1; $i<$number; $i++){
	$file=sprintf("%s.gz", $list->[$i]);
	if(!defined($options->{"--silent"})){
	    print stderr "$i: opening; $file\n";
	}
	if(open FIN2, "$CAT $file | $GUNZIP | " ){
	    while($line=<FIN2>){
		if($outfile ne ""){
		    print FOUT $line;
		}else{
		    print $line;
		}
	    }
#	    print stderr "$i: done!\n";
	    close FIN2;
#	    sleep 1;
#	    print stderr "$i: closed!\n";
	}else{
	    die("Could not open result file: $file\n");
	}
    }
    close FOUT;
}
