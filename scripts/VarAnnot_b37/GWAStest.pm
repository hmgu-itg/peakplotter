package GWAStest;

=head1 Description

This module was developed as part of our custom built variant functional annotation tool.
Its role to find previously reported GWAS signals in the vicinity of the queried variant based on its position.

=head1 Function

The script uses the previously filled choromosome, start, end positions will of the data table.
These values will be used to find
GWAS signals.

=head1 Requirements

    bedtools in the path, read permission in the following folder:
    /nfs/team144/ds26/FunctionalAnnotation/20150505_new_tool/data/
    No internet connection or external libraries are needed.

=head1 SYNOPSIS

    use GWAStest;
    my %hash = ();
    # Input lines are stored in a hash: $hash{$line_no}{"input" = "input line"}
    %hash = %{GWAStest::testGWAScatalog(\%hash, $distance)} # $distance is optional! By default, it is 50kbp

=head1 Input:
    $hash{$line_no}{
        'chr'   => <X>,              # chromosome
        'start' => <position>,       # start position of the variant
        'end'   => <position>,       # end position of the variant
    }

=head1 Output2

    $hash{$line_no}{
        'chr'   => <X>,              # chromosome
        'start' => <position>,       # start position of the variant
        'end'   => <position>,       # end position of the variant
        'GWAS_trait' => <traits>     # all traits in the neighbourhood separated by a semicolon
        'GWAS_pval' => <p_value>     # lowest p value of all traits found.
    }

=head1 Output2

Besides the updated hash, the module will save all data into a separate file. Where the format will be the following:
Of course one queried SNP might have more than one line if multiple gwas signals are associated.

B<SNP_id/rsID trait signal_SNP p_value Pubmed_ID>

=head1 Generation of the test dataset

The test data was downloaded from the GWAS catalog. Then pooled together with our
anthropometric positive control list. GWAS hits overlapping at the same position were
also pooled, keeping the lowest p value.

=head1 Todo list

The test dataset needs to be cleaned further: there are way much more trait definitions
that necessary. I will manually go through all definitions to make the dataset more
homogenous. I also have to create a toolset to automatically extend the dataset with
new positive control lists. The dataset has to be further filtered to remove
signals with p values higher than 5e-8.

=cut


# Main method of the class. Loops through the hash, calls intersectbed,

use strict;
use warnings;

sub testGWAScatalog {

    my %hash      = %{$_[0]};
    my @AllFields = @{$_[2]};

    my @GWAS_fields = ("GWAS_hits");
    push (@AllFields, @GWAS_fields);

    # The distane is default +/-500kb from the variant. It can be changed by
    # submitting a second parameter.
    my $distance = 500000;
    $distance = $_[1] if $_[1];

    # This file has to be updated any time the positive control list is updated:
    my $GWAS_file = "/nfs/team144/ds26/FunctionalAnnotation/20150505_new_tool/data/gwas_catalog_20150616.bed"; # Updated gwas file with the Teslovich dataset.

    # Checking if bedtools are in the path:
    if (`which bedtools` =~ /bedtools/) {
        `bedtools  --version` =~ /^bedtools (.+)\n/;
        print STDERR "[Info] Bedtools version: $1\n";
    } else {
        warn "[Error] Bedtools did not find in path! Testing overlapping with GWAS hits skipped!.\n";
        return \%hash;
    }

    # looping through all variants:
    print STDERR "[Info] Max distance of GWAS signals from variants: $distance bp\n";
    print STDERR "[Info] GWAS datafile: $GWAS_file\n";
    print STDERR "[Info] Contains: GWAS catalog (2015.06.08), anthropometric positive controls (2015.06.10), osteoarthritis positive controls (2015.06.16)\n";
    print STDERR "[Info] Checking GWAS signals.";

    # Opening file, and writing header:
    open(my $GWAS_FILE, ">./gwas_signals.tsv") or die "GWAS output file could not be opened in the current working folder!\n";
    print $GWAS_FILE "Input\tTrait\trsID\tchr:pos\tDistance\tp-value\tPMID\n";

    foreach my $keys (keys %hash){
        print STDERR ".";

        # Creating a bed formatted line for each variant:
        my $lower = $hash{$keys}{start} - $distance; # There might be problem if the position of the variant is smaller than the window.
        $lower = 0 if $lower < 0;

        my $bed_line  = sprintf ("chr%s\t%s\t%s", $hash{$keys}{chr},  $lower, $hash{$keys}{end} + $distance);

        # Run intersect bed:
        my $intersect_GWAS =  `echo "$bed_line" | intersectBed -wb -a stdin -b $GWAS_file`;

        # Process output:
        %hash = %{&parse_GWAS_output(\%hash, $keys, $intersect_GWAS, $GWAS_FILE)};
    }

    print STDERR " done\n";
    print STDERR "[Info] Gwas signals will be saved into a separate file in the current working directory: ./gwas_signals.tsv\n";
    return (\%hash, \@AllFields);
}


# This helper routine parses the output.
sub parse_GWAS_output{
    my %hash             = %{$_[0]}; # The has with all data.
    my $key              = $_[1]; # The actual snp that we are quering.
    my $intersect_output = $_[2]; # From which we read the gwas signals.
    my $GWAS_file        = $_[3]; # File handle of the output file.

    my $variantPosition  = $hash{$key}{"start"};
    my %signals          = ();

    # If there are no GWAS signal in the neighbourhood:
    if ($intersect_output eq ""){
        $hash{$key}{"GWAS_hits"} = "-";
        return \%hash;
    }

    # There are overlapping GWAS signals:
    else {
        my @pvals  = ();
        my %trait_hash = ();

        # parse output lines:
        my @lines = split(/\n/, $intersect_output);
        foreach my $line (@lines){
            my @fields = split(/\t/, $line);

            my $chr   = $fields[3];
            my $pos   = $fields[4];
            my $rsid  = $fields[6];
            my $trait = $fields[7];
            my $pval  = $fields[8];
            my $pmid  = $fields[9];
            my $distance = abs($pos - $variantPosition);

            if ( exists $signals{$trait}) {
                $signals{$trait} = [$chr, $pos, $rsid, $pval, $pmid, $distance] if $distance < $signals{$trait}->[5];
            }
            else {
                $signals{$trait} = [$chr, $pos, $rsid, $pval, $pmid, $distance]
            }
        }

        my $trait_string = "";
        foreach my $trait (keys %signals){
            my @array = @{$signals{$trait}};
            my $string = sprintf("trait:%s rsID:%s pval:%s PMID:%s distance:%s SNPID:%s:%s", $trait, $array[2], $array[3], $array[4], $array[5], $array[0], $array[1]);
            $trait_string .= $string."|";

            # Besides creating the GWAS entry in the hash, we save the data into a separate file:
            #                 Input  Trait   rsID     chr:pos             Distance    p-value     PMID\n"
            my $input_string = $hash{$key}->{"input"};
            print $GWAS_file "$input_string\t$trait\t$array[2]\t$array[0]:$array[1]\t$array[5]\t$array[3]\t$array[4]\n";
        }

        $hash{$key}{"GWAS_hits"} = $trait_string;
    }

    return \%hash;
}

1;