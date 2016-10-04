package BasicInformation;

=head1 Description

This module was developed as the first step of our custom built variant functional annotation tool.

=head1 Function

The script processes variant id-s stored in the input field of the provided hash. Checks if the
given snpID has a known rsID. Checks the chromosome, position and alleles of variants with valid rsID.
These tests are performed via querying the REST API of Ensembl.

=head1 Requirements

    Modules: JSON, HTTP::Tiny
    Working internet connection to access the REST API of Ensembl

=head1 SYNOPSIS

    use BasicInformation;
    my %hash = ();
    # Input lines are stored in a hash: $hash{$line_no}{"input" = "input line"}
    %hash = %{BasicInformation::Input_cleaner(\%hash)}

=head1 Output

    $hash{$line_no}{
        'input' => <rsid or snp id>, # From the input
        'rsID'  => <rs####>          # rsID based on input or the snp id.
        'chr'   => <X>,              # chromosome
        'start' => <position>,       # start position of the variant
        'end'   => <position>,       # end position of the variant
        'ref'   => <residue>,        # reference allele
        'alt'   => <residue>,        # alternative allele
        'ancestral' => <residue>     # Ancestral allele, which not necessarely means the reference.
    }


=head1 Todo

Although the main script is responsible to make sure the data provided by the user
is in the requested format, but I want to add validity checks to make the script more robust.

=cut

use strict;
use warnings;
use HTTP::Tiny;
use JSON;

# Input cleaner is the primary method of the class. It call other helper functions to
# get .
sub Input_cleaner {

    print STDERR "[Info] Checking input data."; # Status update indicating each variant being checked.

    my %input_lines = %{$_[0]};

    # The following array contains the list of the processed fields, in the proper order:
    my @fields = ("input", "chr", "start", "end", "rsID", "ref", "alt", "ancestral_allele", "var_class");

    # At first we check if rsid exist, where it's missing from the input.
    foreach my $key (keys %input_lines){

        print STDERR ".";

        if ( $input_lines{$key}{"input"} =~ /rs\d+/i){
            $input_lines{$key}{"rsID"} = $input_lines{$key}{"input"}; # At this point it is a temporaty
        }
        else {
            $input_lines{$key}{"rsID"} = &rs_finder($input_lines{$key}{"input"});
        }
    }


    # The following lines are executed if we want to add chr, pos and alleles
    foreach my $key (keys %input_lines){

        print STDERR ".";

        # If the rsid field has a valid rsid, then the check out the details on Ensembl
        if ( $input_lines{$key}{"rsID"} =~ /rs\d+/i){
            $input_lines{$key} = &pos_finder_rs($input_lines{$key});
        }
        else {
            $input_lines{$key} = &pos_finder($input_lines{$key});
        }
    }

    print STDERR " done\n";
    return (\%input_lines, \@fields);
}



sub pos_finder_rs {
    # If the variant has a real rsid, then we check the ensembl for real positional data

    my %hash = %{$_[0]};

    my $rsid = $hash{"rsID"};
    my $URL = sprintf("http://grch37.rest.ensembl.org/variation/human/%s?content-type=application/json", $rsid);

    my $http = HTTP::Tiny->new();
    my $response = $http->get($URL,{
        headers => { 'Content-type' => 'application/json' }
    });

    # This short loop was added to make the code more robuts. We don't want it to crash any time,
    # there is some problem with downloading the data.
    my $fail_count = 0;
    unless ($response->{success}){
        while ($fail_count < 20) {
            sleep(5);
            print STDERR "[Warning] Downloading data from the Ensembl server failed. Trying again.\n";
            my $http     = HTTP::Tiny->new();
            my $response = $http->get($URL,{});
            $fail_count++;
        }

        # If more than 20 fails accumulated, we proceed to the next variant.
        print STDERR "[Warning] Downloading data from the Ensembl server failed. variant: $rsid. URL: $URL\n";
        next;
    }


    my $data = decode_json($response->{content});

    # Setting un-ambigous values:
    $hash{"rsID"}  = $data->{"name"};
    $hash{"ancestral_allele"} = $data->{"ancestral_allele"} || "-";
    $hash{"var_class"} = $data->{"var_class"} || "-";

    foreach my $mapping (@{$data->{"mappings"}}) {
        next if length($mapping->{"seq_region_name"}) > 2;

        $hash{"chr"}   = $mapping->{"seq_region_name"};
        $hash{"start"} = $mapping->{"start"};
        $hash{"end"}   = $mapping->{"end"};

        my $allele_string = $mapping->{"allele_string"};
        my ($ref, $alt) = ($allele_string =~ (/(\S+)\/(\S+)/));

        $hash{"alt"} = $alt;
        $hash{"ref"} = $ref;
    }

    return \%hash;
}

# If there is no rsid, then we just parse the input field
sub pos_finder {
    my %hash = %{$_[0]};
    my $input = $hash{"input"};
    $input =~ /chr(\d{1,2})\:(\d+)\-(\d+)\_(\S+)_(\S+)/i;
    my ($chr, $start, $end, $ref, $alt) = ($1, $2, $3, $4, $5);

    # Retrieving the real reference allele, that is based on the reference genome:
    my $realref = _getRefAllele($chr, $start, $end);

    # Updating $ref and $alt
    if ($realref eq $alt) { # It means that the alleles are swapped:
        $alt = $ref;
        $ref = $realref;
    }
    elsif ( $realref eq $ref ){
        # That's good, there's nothing to do.
    }
    else {
        # We only give a warning about the discrepancy if it is not observed in indels.
        print "[Warning] There migh be some problem with $input: neither alleles are matching with the real reference: $realref\n" unless ($alt eq "-") or ($ref eq "-");
    }

    my $var_class = "SNP";
    $var_class = "indel" if ($ref eq "-" or $alt eq "-");
    $var_class = "indel" if (length($ref) >1 or length($alt) > 1);

    # At this point there are no checkings, as we assume the input
    # field is properly formatted.
    $hash{"chr"}   = $chr;
    $hash{"end"}   = $end;
    $hash{"start"} = $start;
    $hash{"alt"}   = $alt;
    $hash{"ref"}   = $ref;
    $hash{"ancestral_allele"} = $ref;
    $hash{"var_class"} = $var_class;

    return \%hash;

}

sub _getRefAllele {

    my ($chr, $start, $end) = @_;

     # If we are looking at an insertion, the start-end values are a bit messed:
    my $URL = '';
    if ($end < $start) {
        $URL = sprintf("http://grch37.rest.ensembl.org/sequence/region/human/%s\:%s..%s?feature=variation", $chr, $end, $start);
    }
    else {
        $URL = sprintf("http://grch37.rest.ensembl.org/sequence/region/human/%s\:%s..%s?feature=variation", $chr, $start, $end);
    }

    my $http = HTTP::Tiny->new();

    my $response = $http->get($URL,{
        headers => { 'Content-type' => 'application/json' }
    });

    # This short loop was added to make the code more robuts. We don't want it to crash any time,
    # there is some problem with downloading the data.
    my $fail_count = 0;
    unless ($response->{success}){
        while ($fail_count < 3) {
            sleep(5);
            print STDERR "[Warning] Downloading data from the Ensembl server failed. Trying again.\n";
            my $http     = HTTP::Tiny->new();
            my $response = $http->get($URL,{});
            $fail_count++;
        }

        # If more than 20 fails accumulated, we proceed to the next variant.
        print STDERR "[Warning] Downloading data from the Ensembl server failed. URL: $URL\n";
        next;
    }

    my $hash = decode_json($response->{content});

    return $hash->{seq} || return "-";

}

sub rs_finder {
    my $input = $_[0];

    $input =~ /chr(\d{1,2})\:(\d+)\-(\d+)\_(\S+)_(\S+)/i;
    my ($chr, $start, $end, $alt, $ref) = ($1, $2, $3, $4, $5);

    # Retrieving the real reference allele, that is based on the reference genome:
    my $realref = _getRefAllele($chr, $start, $end);

    # Updating $ref and $alt
    if ($realref eq $alt) { # It means that the alleles are swapped:
        $alt = $ref;
        $ref = $realref;
    }
    elsif ( $realref eq $ref ){
        # That's good, there's nothing to do.
    }
    else {
        print STDERR "\n\t[Warning] There migh be some problem with $input: neither alleles are matching with the real reference: $realref"  unless ($alt eq "-") or ($ref eq "-");
    }


    # If we are looking at an insertion, the start-end values are a bit messed:
    my $URL = '';
    if ($end < $start) {
        $URL = sprintf("http://grch37.rest.ensembl.org/overlap/region/human/%s\:%s-%s?feature=variation", $chr, $end, $start);
    }
    else {
        $URL = sprintf("http://grch37.rest.ensembl.org/overlap/region/human/%s\:%s-%s?feature=variation", $chr, $start, $end);
    }

    my $http = HTTP::Tiny->new();

    my $response = $http->get($URL,{
        headers => { 'Content-type' => 'application/json' }
    });

    # print "\n$URL\n";

    # This short loop was added to make the code more robuts. We don't want it to crash any time,
    # there is some problem with downloading the data.
    my $fail_count = 0;
    unless ($response->{success}){
        while ($fail_count < 3) {
            sleep(5);
            print STDERR "\n\t[Warning] Downloading data from the Ensembl server failed. Trying again.";
            my $http     = HTTP::Tiny->new();
            my $response = $http->get($URL,{});
            $fail_count++;
        }

        # If more than 20 fails accumulated, we proceed to the next variant.
        print STDERR "\n\t[Warning] Downloading data from the Ensembl server failed. URL: $URL";
        next;
    }

    my $hash = decode_json($response->{content});

    # I have noticed that the overlapping variants sometimes are not accessible at the first query.
    # Therefore I don't trust the first negative answer and repeat the qurey once more.
    unless ( $hash->[0]->{id}) {
        # print STDERR "\n\t[Warning] Returned value contain no data. Repeating download.";
        for (1 .. 10){
            # print STDERR "try ...";
            $response = $http->get($URL,{
                headers => { 'Content-type' => 'application/json' }
            });
            $hash = decode_json($response->{content});
            last if $hash->[0]->{id};
        }
    }

    # If the second access still remained unsuccessful we return the empty value.
    return "-" unless $hash->[0]->{id};

    # If there are overlapping variations, each overlapping variant will be tested if the alleles are match:
    foreach my $variant (@{$hash}){

        # Reading values from the returned data:
        my $A1    = $variant->{alleles}->[0] if $variant->{alleles}->[0];
        my $A2    = $variant->{alleles}->[1] if $variant->{alleles}->[1];
        my $Start = $variant->{start} if $variant->{start};
        my $End   = $variant->{end} if $variant->{end};
        my $name  = $variant->{id};

        # If we found the alleles we are looking for:
        if ((($A1 eq $ref) and ($A2 eq $alt)) and (($start == $Start) and ($end == $End))) {
            return $name;
        }
        # print STDERR "\n[Warning] Checking overlapping variant: $name (chr$chr:$Start\_$End\_$A1/$A2) does not match with quried variant: (chr$chr:$start\_$end\_$ref/$alt)\n";
    }
    # If we don't find anything like the variant we were looking for:
    print STDERR "\n\t[Warning] There were overlapping variants, but the alleles did not match! rsID won't be considered.";
    return "-";
}
1;
