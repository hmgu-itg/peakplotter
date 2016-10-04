package getGWAVA;

use strict;
use warnings;

# This script generates the following values:
    # gwava score,
    # gerp score
    # average gerp score
    # DNAse
    # bound motifs
    #


sub GWAVA {

    print STDERR "[Info] Now, performing GWAVA calculations.\n";

    # Hardwired paths to python folders, and gwava directory.
    my $pythonDir = '/nfs/team144/software/anaconda/bin/:/nfs/team144/software/anaconda/lib/python2.7/site-packages/';
    my $gwavaDir  = '/lustre/scratch113/teams/zeggini/users/ds26/GWAVA/gwava_release';

    # Checking samtools
    my $samtools = `which samtools`;

    # input variables have been processed:
    my %hash      = %{$_[0]};
    my @AllFields = @{$_[1]};

    unless ($samtools){
        print STDERR "[Warning] Samtools was not found in path! Use: 'module add hgi/samtools/latest'\n";
        print STDERR "[Warning] GWAVA calculations will be skipped. Returning to main program.\n";
        return (\%hash, \@AllFields);
    }

    # Generate random file names
    my $temporarySuffix = int(rand(10000));
    my $tempBedfile     = sprintf("/tmp/varAnnotation_GWAVA_input_%s.bed",  $temporarySuffix);
    my $tempGWAVAannot  = sprintf("/tmp/varAnnotation_GWAVA_annot_%s.csv",  $temporarySuffix);
    my $tempGWAVAscore  = sprintf("/tmp/varAnnotation_GWAVA_scores_%s.bed", $temporarySuffix);

    # If we can't open the output file for read, we return without doing anything.
    my $OUTBED;
    unless ( open($OUTBED, ">", $tempBedfile )){
        return (\%hash, \@AllFields);
        print STDERR  "[Error] Temporary file could not be oppened: $tempBedfile\n";
    }

    # Looping through all input entries:
    foreach my $key (keys %hash){

        # The corresponding fields will be extracted from the main hash:
        my $chr   = $hash{$key}->{"chr"};
        my $start = $hash{$key}->{"start"};
        my $end   = $hash{$key}->{"end"};
        my $input = $hash{$key}->{"input"};

        # all entries are saved in a temporary file.
        printf $OUTBED "chr%s\t%s\t%s\t%s", $chr, $start-1, $end-1, $input;
    }

    # Sorting input bed file:
    `sort -k1,1 -k2,2n $tempBedfile -o $tempBedfile`;

    # run gwava annotation script:
    `export PYTHONPATH=\$PYTHONPATH:$pythonDir export GWAVA_DIR=$gwavaDir python $gwavaDir/src/gwava_annotate.py $tempBedfile $tempGWAVAannot`;

    # Checking if the annotation file was created:
    unless (-e $tempGWAVAannot){
        print STDERR "[Error] The GWAVA annotation scrip did not run successfully. The annotation file ($tempGWAVAannot) was not created.\n";
        return (\%hash, \@AllFields);
    }

    # run gwava:
    `export PYTHONPATH=\$PYTHONPATH:$pythonDir export GWAVA_DIR=$gwavaDir python $gwavaDir/src/gwava.py tss $tempGWAVAannot $tempGWAVAscore`;

    # Checking if the GWAVA run was successful or not:
    unless (-e $tempGWAVAscore){
        print STDERR "[Error] The GWAVA scrip did not run successfully. The score file ($tempGWAVAscore) was not created.\n";
        return (\%hash, \@AllFields);
    }

    # Process output file. Update hash
    my ($hash, $AllFields) = &parse_GWAVA_File(\%hash, \@AllFields, $tempGWAVAscore);

    # remove temporary files.
    `rm $tempBedfile`;
    `rm $tempGWAVAannot`;
    `rm $tempGWAVAscore`;

    # Returning data:
    return ($hash, $AllFields);
}


sub parse_GWAVA_File{
    # input variables have been processed:
    my %hash       = %{$_[0]};
    my @AllFields  = @{$_[1]};
    my $gwava_file = $_[2];

    # The following fields will be extracted from the returned data:
    my %GWAVA_fields = ("GWAVA_score" => "NA",
                        "gerp"        => "NA",
                        "avg_gerp"    => "NA",
                        "dnase_fps"   => "NA",
                        "bound_motifs"=> "NA",
                        "DNase"       => "NA");

    push (@AllFields, keys(%GWAVA_fields));

    # Reading gwava file, fill a hash with all data:
    open(my $GWAVAFILE, "<", $gwava_file);
    while ( my $line = <$GWAVAFILE>) {
        chomp $line;
        my %indices = &getindex($line, \%GWAVA_fields) if $. == 1;
    }

    # Loop through all

}

sub getindex {

    my $

}

1;