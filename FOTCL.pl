#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Wed Aug 28 14:26:20 2013
###################################
use strict;
use warnings;
use FindBin;
use Getopt::Long;

# Feature Occupancy Transition between Cell Lines (FOTCL)

my $start_time = time();

=head1 Program Description:
Given a list of genome interval of interest, try to figure out which features (TF/DNase/Histone) changes 
significant between to cell lines (loss of gain)

=cut


my $USAGE_DESCRIPTION = "Usage: perl $0 -i input.bed -o output.txt -a cell1 -b cell2 -flank 0
    -i  input file in bed format (only the first 3 colums will be used)
    -o  ouput file
    -a  cell line 1
    -b  cell line 2
    --flank overlap with feature allowing gaps, default 0, no gaps\n\n";

############# PARAMETER VALUES ##############
my $INPUT_FILENAME = "";      # genome intervals of interest
my $OUTPUT_FILENAME = "";     
my $CELL1 = "";               # cell line1 name from ENCODE
my $CELL2 = "";               # cell line2 name from ENCODE
my $HELP=0;
my $FLANK_SIZE = 0;

############# GLOBAL VARIABLES ##############
my %tfdnase = ();             #cell=>tf=>cell|tf|acc|lab|treatment=>filename 
my %features = ();
my %featureOverlapsForPrinting = ();
my %result = ();

my $CELL1_TOTAL = 0;   #how many intervals there are total in your cell1 file
my $CELL2_TOTAL = 0;   #how many intervals there are total in your cell2 file
my $INPUTTOTAL=0;

############ MAIN EXECUTION ##############
parseArguments();

prepareFeatures();

# Take CELL1 as cell1
my $cell1;
open($cell1,">cell1_tmp.bed3") or die "Can't create temp output file cell1_tmp.bed3\n";
prepareInInputFile($cell1,10000,$tfdnase{$CELL1},$CELL1);
close $cell1;

# Take CELL2 AS cell2
my $cell2;
open($cell2,">cell2_tmp.bed3") or die "Can't create temp output file cell2_tmp.bed3\n";
prepareInInputFile($cell2,10000,$tfdnase{$CELL2},$CELL2);
close $cell2;

sortSubsetAndBackgroundFiles();

intersectAndWriteToFile();

calculateFotclScores();

printResults();

`rm cell2.bed3`;
`rm cell1.bed3`;
`rm input.bed3`;

&printRunTime();

############## SUBROUTINES ###############
# Given cell1 and cell2, checking which TF/DNase are available
# Read from data/files.txt file
sub prepareFeatures{
   my $infile="$FindBin::Bin/data/files.txt";
   #my $infile="$FindBin::Bin/data/BroadHistone_files.txt";
   if(! -e $infile){
       die "TF/DNase/Histone peaks haven't been initialized, please use initialized.pl to download TF/DNase/Histone from ENCODE\n";
   }
   open(IN,$infile) or die $!;
   my %hash;
   my %antibody;
   while(<IN>){
       my($filename,@other) = split /\t|;/;
       next if ($filename!~/Peak\.gz$/);
       my $cell="";
       my $lab="";
       my $antibody="DNase"; #default is DNase
       my $treatment="";
       my $dccAccession="";
       foreach my $a (@other){
            if($a=~/cell=(.*)/){
                $cell=uc($1);
            }elsif($a=~/lab=(.*)/){
                $lab=$1;
            }elsif($a=~/antibody=(.*)/){
                $antibody=$1;
            }elsif($a=~/treatment=(.*)/){
                $treatment=$1;
            }elsif($a=~/dccAccession=(.*)/){
                $dccAccession=$1;
            }
       }
       next if ($cell eq "");
       next if (! -e "$FindBin::Bin/data/$filename");
       my $key=join "|",($cell,$antibody,$dccAccession,$lab,$treatment);
       $antibody{$antibody}=1;
       $hash{$cell}->{$antibody}->{$key}="$FindBin::Bin/data/$filename";
   }
   close IN;
   
   # Checking which TF and DNase are both avaiable for $CELL1 and $CELL2
   my %cell1=();
   my %cell2=();
   %cell1=%{$hash{$CELL1}} if (exists($hash{$CELL1})); 
   %cell2=%{$hash{$CELL2}} if (exists($hash{$CELL2})); 
   my @tfs = grep {exists($cell1{$_})} keys %cell2;  # common tfs

   # write out data avalible matrix;
   open(OUT,">DataResourceMatrix.csv") or die $!;
   print OUT join ",",("",sort keys %antibody);
   print OUT "\n";
   foreach my $cell (sort keys %hash){
         print OUT $cell;
        foreach my $antibody (sort keys %antibody){
            my $flag=0;
            $flag=1 if exists($hash{$cell}->{$antibody});
            print OUT ",$flag";
        }
        print OUT "\n";
   }
   close OUT;

   if(scalar(@tfs)){
        foreach my $tf (@tfs){
            $tfdnase{$CELL1}->{$tf}=$cell1{$tf};
            $tfdnase{$CELL2}->{$tf}=$cell2{$tf};
        }
   }else{
        # No common TF/DNase, exit
        print STDERR "There are no TF/DNase/Histone availble for both '$CELL1' and '$CELL2', Please check the DataResourceMatrix.csv \n";
        exit 1;
   }
}

=head1 Expect the following format for the intersectBed file
chr20	307680	307739	chr20	307680	307739	NHEK_H3K27me3_encodeHistone
=cut
sub processOverlapFile{
	my $fileToReadIn = shift(@_);
	my $whichOverlap = shift(@_);
	open(IN,"<",$fileToReadIn) or die "Can't open $fileToReadIn\n";
	while(<IN>)
	{
		chomp;
		my @fields = split/\t/;
		my $numFields = scalar(@fields);
		die "Expecting 7 intersectBed output fields, instead got this: $_\n" unless ($numFields == 7);
		my $feature = $fields[3];
        my ($cell,$tf,$acc,$lab,$treat) = split /\|/,$feature;
		die "Expecting the last field of the overlap field to be a name of a feature.  Instead got this: $feature\n" unless $feature =~ /\w/i;
#####################################
# PUT IN IF YOU SEE A FEATURE THAT YOU ARE TRYING TO KEEP
#####################################
		my $testInterval = $fields[0] . "_" . $fields[1] . "_" . $fields[2];
        my $inputInterval = join "_",(@fields[4..6]);

		next if (exists($features{$tf}->{$cell}->{$feature}->{$inputInterval}));
		if (!exists($features{$tf})) {die "For some reason, this feature wasn't in featurelist file: $feature\n";}
		else {
			$features{$tf}->{$cell}->{$feature}->{$inputInterval} = 1;
			$features{$tf}->{$cell}->{$feature}->{$whichOverlap}++;
		} 
	}
	close IN;
}


sub calculateFotclScores {
	print STDERR "Calculating enrichment scores...";
    
    foreach my $tf (%features){
        foreach my $key1 (keys %{$features{$tf}->{$CELL1}}){
            my $cell1total = $features{$tf}->{$CELL1}->{$key1}->{TOTAL};
            my $cell1overlap = $features{$tf}->{$CELL1}->{$key1}->{'CELL1_OVERLAP_NUMBER'};
            my $filename1 = $features{$tf}->{$CELL1}->{$key1}->{'FILENAME'};
            $cell1overlap=0 unless($cell1overlap);
            foreach my $key2(keys %{$features{$tf}->{$CELL2}}){
                my $cell2total = $features{$tf}->{$CELL2}->{$key2}->{'TOTAL'};
                my $filename2 = $features{$tf}->{$CELL2}->{$key2}->{'FILENAME'};
                my $cell2overlap = $features{$tf}->{$CELL2}->{$key2}->{'CELL2_OVERLAP_NUMBER'};
                $cell2overlap=0 unless($cell2overlap);
                my $str=join "\t",($tf,$INPUTTOTAL,$cell1overlap,$cell2overlap,$cell1total,$cell2total,$key1,$key2,$filename1,$filename2);
                $result{$str}=1;
            }
        }
    }
	print STDERR "done!\n";
}

#quick subroutine to print a header onto the output files
sub printResults
{
	open(OUT,">",$OUTPUT_FILENAME) or die "Can't open $OUTPUT_FILENAME\n";
	print OUT "Feature\tInputTotal\tOverlaps with $CELL1\tOverlaps with $CELL2\t$CELL1 Total\t$CELL2 Total\t$CELL1 dataset\t$CELL2 dataset\t$CELL1 filename\t$CELL2 filename\n";
	foreach my $feature (sort keys %result){
	    print OUT $feature;
        print OUT "\n";
    }
	close OUT;
}

sub intersectAndWriteToFile {
	print STDERR "Overlapping Input with cell $CELL1 ...";
	`intersectBed -wa -wb -sorted -a cell1.bed3 -b input.bed3 > cell1_overlap.temp`;
	print STDERR "done!\n";	
	print STDERR "Processing cell $CELL1  overlaps...";
	processOverlapFile("cell1_overlap.temp","CELL1_OVERLAP_NUMBER");
	print STDERR "done!\n";
	`rm cell1_overlap.temp`;
	
	print STDERR "Overlapping Input with cell $CELL2...";
	`intersectBed -wa -wb -sorted -a cell2.bed3 -b input.bed3 > cell2_overlap.temp`;
	print STDERR "done!\n";
	print STDERR "Processing cell $CELL2 overlaps...";
	processOverlapFile("cell2_overlap.temp","CELL2_OVERLAP_NUMBER");
	print STDERR "done!\n";
	`rm cell2_overlap.temp`;
}

sub sortSubsetAndBackgroundFiles
{
	print STDERR "Sorting ...";
	`sort -k 1,1 -k2,2n cell1_tmp.bed3 > cell1.bed3`;
	`rm cell1_tmp.bed3`;
	
	`sort -k 1,1 -k2,2n cell2_tmp.bed3 > cell2.bed3`;
	`rm cell2_tmp.bed3`;

    `cut -f1-3 $INPUT_FILENAME |sort -k 1,1 -k2,2n|uniq >input.bed3`;
     $INPUTTOTAL=`wc -l input.bed3 |awk '{print \$1}'`;
     chomp($INPUTTOTAL);
	 print STDERR "done!\n";
}


sub parseAndStoreFeatureList
{
	my $input = shift(@_);
	$features{$input}->{"CELL1_OVERLAP_NUMBER"} = 0;
	$features{$input}->{"CELL2_OVERLAP_NUMBER"} = 0;
}

# Read in input file 
# First the file handle, second is the denominator, 3rd is the hash reference
sub prepareInInputFile
{
	my $fh = shift(@_);
	my $COUNTER_DENOMINATOR = shift(@_); #every this many lines, print output to user (a dot)
	my $ref= shift(@_); #signal for which parsing subroutine to use
    my $cell=shift(@_);
	
	print STDERR "Prepare TFs/DNase/Histone of $cell ..";
	my $counter = 0;
	
	foreach my $tf (sort keys %{$ref}){
        foreach my $key (sort keys %{$ref->{$tf}}){
            open(IN,"gunzip -c $ref->{$tf}->{$key} |") or die $!;
            my $num=0;
            while(<IN>){
		        chomp;
		        if ($counter % $COUNTER_DENOMINATOR == 0) {print STDERR ".";}
	            my @fields = split/\t/;
            	my $chr = $fields[0];
            	my $start = $fields[1];
                my $stop = $fields[2];
	
                #GENERAL INTERVAL FORMATTING CHECK
	            next unless ( ($chr =~ /^chr/) && ($start =~ /^\d+$/) && ($stop =~ /^\d+$/) && ($stop >= $start));
	            $start -= $FLANK_SIZE;
	            $stop += $FLANK_SIZE;
	            print $fh "$chr\t$start\t$stop\t$key\n";
		        $counter++;	
                $num++;
	        }
            my ($cell,@ooo) = split /\|/,$key;
            $features{$tf}->{$cell}->{$key}->{'TOTAL'}=$num;
            my $filename=$ref->{$tf}->{$key};
            $filename=~s/.*\///g;
            $filename=~s/\.narrowPeak\.gz//g;
            $features{$tf}->{$cell}->{$key}->{'FILENAME'}=$filename;
            close IN;
        }
    }
	print STDERR "done!\n";
}

#parse arguments for the program
sub parseArguments{
	
    unless(GetOptions(
            "i=s"=>\$INPUT_FILENAME,
            "o=s"=>\$OUTPUT_FILENAME,
            "a=s"=>\$CELL1,
            "b=s"=>\$CELL2,
            "flank=i"=>\$FLANK_SIZE,
            "h|help"=>\$HELP)){
        die  $!,"\n$USAGE_DESCRIPTION\n";
    }
    if($HELP){
        print $USAGE_DESCRIPTION,"\n";
        exit(0);
    }

    $CELL1=uc($CELL1);
    $CELL2=uc($CELL2);

    #unless($INPUT_FILENAME) {die "Program quit because it didn't recognize a cell1 file from the input.\nYou can call FOTCL.pl with no arguments for usage details.\n";}
    #unless ($OUTPUT_FILENAME) {die "Program quit because it didn't recognize an output filename from the input.\nYou can call FOTCL.pl with no arguments for usage details.\n";}
    #unless ($CELL1) {die "Program quit because it didn't recognize a feature filename from the input.\nYou can call FOTCL.pl with no arguments for usage details.\n";}
    #unless ($CELL2) {die "Program quit because it didn't recognize a featurelist filename from the input.\nYou can call FOTCL.pl with no arguments for usage details.\n";}
	unless($INPUT_FILENAME) {die "Program quit because it didn't recognize an input filename from the input.\n$USAGE_DESCRIPTION";}
	unless ($OUTPUT_FILENAME) {die "Program quit because it didn't recognize an output filename from the input.\n$USAGE_DESCRIPTION";}
	unless ($CELL1) {die "Program quit because it didn't recognize a cell line 1 from the input.\n$USAGE_DESCRIPTION";}
	unless ($CELL2) {die "Program quit because it didn't recognize a cell line 2 from the input.\n$USAGE_DESCRIPTION";}
}

sub printRunTime
{
	my $end_time = time();
	my $run_time = $end_time - $start_time;
	print STDERR "Done!\nThis job took $run_time seconds.\n";
}

################# END ####################
