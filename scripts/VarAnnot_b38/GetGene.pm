package GetGene;

=head1 Description

This module was developed as part of our custom built variant functional annotation tool.
Its role to download information of gene a given variant is overlapping with

=head1 Function

The script uses the previously filled choromosome, start, end positions will of the data table.
These values will be used to query the Ensembl through the REST API.

=head1 Requirements

    modules: JSON, HTTP:Tiny
    Internet connection is also required to reach the REST API

=head1 SYNOPSIS

    use GetGene ".";
    my %hash = ();
    # Input lines are stored in a hash: $hash{$line_no}{"input" = "input line"}
    %hash = %{GetGene::Gene(\%hash)}

=head1 Input:
    $hash{$line_no}{
        'chr'   => <X>,              # chromosome
        'start' => <position>,       # start position of the variant
        'end'   => <position>,       # end position of the variant
    }

=head1 Output

    $hash{$line_no}{
        'chr'              => <X>,         # chromosome
        'start'            => <position>,  # start position of the variant
        'end'              => <position>,  # end position of the variant
        'Gene_name'        => "-", # External name of the gene
        'Gene_id'          => "-", # The Ensembl stable gene ID
        'Gene_description' => "-", #
        'Gene_biotype'     => "-", #
        'Gene_strand'      => "-"  # + or - direction of the transcription
    }


=head1 Generation of the test dataset

The test data was downloaded from the GWAS catalog. Then pooled together with our
anthropometric positive control list. GWAS hits overlapping at the same position were
also pooled, keeping the lowest p value.

=head1 Todo list

I have to find out what to do if a single variant is ovelapping with more than one genes.
At this time, the script does not care about multiple genes....

=cut


use strict;
use warnings;
use HTTP::Tiny;
use JSON;

our $geneBedfile    = "/nfs/team144/ds26/FunctionalAnnotation/20160907_GRCh38_copy/data/gencode.v25.annotation_sorted.bed";
our $proteinBedfile = "/nfs/team144/ds26/FunctionalAnnotation/20160907_GRCh38_copy/data/gencode.v25.annotation_sorted_protein_coding.bed";

# This routine downloads a list of overlapping genes.
# input:
    # Usual format of the input lines.
    #    hash{line#}->{
    #            "chr"   => #,
    #            "start" => #,
    #            "end"   => #}
sub Gene {
    print STDERR "[Info] Retrieving information of overlapping genes.";

    my %variants = %{$_[0]};
    my @Fields   = @{$_[1]};

    my @GeneFields = ("Gene_name","Gene_id","Gene_description","Gene_biotype","Gene_strand","Closest_gene-name", "Closest_gene-ID", "Closest_gene-distance", "Closest_protein_coding_gene-name", "Closest_protein_coding_gene-ID", "Closest_protein_coding_gene-distance");
    push (@Fields, @GeneFields);


    foreach my $variant (keys %variants){

        print STDERR ".";

        # Initializing values to be returned:
        $variants{$variant}->{"Gene_name"}        = "-";
        $variants{$variant}->{"Gene_id"}          = "-";
        $variants{$variant}->{"Gene_description"} = "-";
        $variants{$variant}->{"Gene_biotype"}     = "-";
        $variants{$variant}->{"Gene_strand"}      = "-";

        $variants{$variant}->{"Closest_gene-name"}                    = "-";
        $variants{$variant}->{"Closest_gene-ID"}                      = "-";
        $variants{$variant}->{"Closest_gene-distance"}                = "-";
        $variants{$variant}->{"Closest_protein_coding_gene-name"}     = "-";
        $variants{$variant}->{"Closest_protein_coding_gene-ID"}       = "-";
        $variants{$variant}->{"Closest_protein_coding_gene-distance"} = "-";

        # Creating a query URL string:
        my $query_string = sprintf("%s:%s-%s", $variants{$variant}->{"chr"}, $variants{$variant}->{"start"}, $variants{$variant}->{"end"});
        $query_string    = sprintf("%s:%s-%s", $variants{$variant}->{"chr"}, $variants{$variant}->{"end"}, $variants{$variant}->{"start"}) if $variants{$variant}->{"start"} > $variants{$variant}->{"end"};

        my $URL = "http://rest.ensembl.org/overlap/region/human/$query_string?feature=gene";
        my $http = HTTP::Tiny->new();

        # Processing returned data:
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
            print STDERR "[Warning] Downloading data from the Ensembl server failed. URL: $URL\n";
            next;
        }

        my $hash = decode_json($response->{content});

        # If there are overlapping gene at this point I don't care just with only one.
        # I know it is an issue to address.
        my $prot_found = "";
        my $gene_found = "";
        foreach my $gene (@{$hash}){
            # We have to check which is the overlapping gene is protein coding:
            # We also have to check what if a position is overlap with two protein coding genes...

            $gene_found ++;

            if ( $gene->{biotype} eq "protein_coding" ){

                $prot_found ++;

                $variants{$variant}->{"Gene_name"}        = $gene->{external_name} if $gene->{external_name};
                $variants{$variant}->{"Gene_description"} = $gene->{description} if $gene->{description};
                $variants{$variant}->{"Gene_biotype"}     = $gene->{biotype} if $gene->{biotype};
                $variants{$variant}->{"Gene_strand"}      = $gene->{strand} if $gene->{strand};
                $variants{$variant}->{"Gene_id"}          = $gene->{id} if $gene->{id};
            }
        }
        # The following lines were added for diagnostic purposes. Now removing them.
        # print STDERR "\n[Warning] The variant ($query_string) overlaps with multiple genes. Will considering only the protein coding.\n" if $gene_found and $gene_found > 1;
        # print STDERR "\n[Warning] The variant ($query_string) overlaps with multiple protein coding gene. Will be considering only the first.\n" if $prot_found and $prot_found > 1;

        # If none of the genes are overlapping, then we just keep the first overlapping gene:
        unless ($prot_found){
            $variants{$variant}->{"Gene_name"}        = $hash->[0]->{external_name} || "-";
            $variants{$variant}->{"Gene_description"} = $hash->[0]->{description}   || "-";
            $variants{$variant}->{"Gene_biotype"}     = $hash->[0]->{biotype}       || "-";
            $variants{$variant}->{"Gene_strand"}      = $hash->[0]->{strand}        || "-";
            $variants{$variant}->{"Gene_id"}          = $hash->[0]->{id}            || "-";
        }

        # So what happen if the variant is not overlapping with:
        my ($gene_name, $gene_id, $gene_distance, $pgene_name, $pgene_id, $pgene_distance) = &nearest_genes($variants{$variant}->{"chr"}, $variants{$variant}->{"start"}, $variants{$variant}->{"end"});
        #print join(" " &nearest_genes($variants{$variant}->{"chr"}, $variants{$variant}->{"start"}, $variants{$variant}->{"end"})),"\n";
        # Filling returned values:
        $variants{$variant}->{"Closest_gene-name"}                    = $gene_name;
        $variants{$variant}->{"Closest_gene-ID"}                      = $gene_id;
        $variants{$variant}->{"Closest_gene-distance"}                = $gene_distance;
        $variants{$variant}->{"Closest_protein_coding_gene-name"}     = $pgene_name;
        $variants{$variant}->{"Closest_protein_coding_gene-ID"}       = $pgene_id;
        $variants{$variant}->{"Closest_protein_coding_gene-distance"} = $pgene_distance;
      }

    print STDERR " done\n";
    return (\%variants, \@Fields);
}


sub nearest_genes {
    my ($chr, $pos1, $pos2) = @_;

    # Checking which position is a bigger:
    my ($start, $end) = ($pos1, $pos2);
    ($start, $end) = ($pos2, $pos1) if $pos1 > $pos2;

    # generate bed formatted string:
    my $bedString = sprintf("chr%s\\t%s\\t%s", $chr, $start, $end);

    # Query gene dataset (all genes):
    my $Closest_gene    = `echo -e \"$bedString\" | closestBed -d -a stdin -b $geneBedfile`;
    my ($gene_name, $gene_id, $gene_distance) = &bedParse($Closest_gene);

    # Query protein gene dataset:
    my $Closest_protein = `echo -e \"$bedString\" | closestBed -d -a stdin -b $proteinBedfile`;
    my ($pgene_name, $pgene_id, $pgene_distance) = &bedParse($Closest_protein);

    return ($gene_name, $gene_id, $gene_distance, $pgene_name, $pgene_id, $pgene_distance);
}

# Once the nearestBed run is completed, we extract important information from the returned line.
    # We save the distance, the name of the gene and the ensembl identifier
sub bedParse {
    my $lines = $_[0];
    my ($geneID, $gene_name,$distance) = ("-","-","-");

    # before processing the lines, we have to be ready to deal with multiple lines.
    # There are possibilities, that the variation is equally distant from multiple genes.
    # In this case all lines will be retrieved. We will select the first protein coding gene.
    my @lines = split("\n", $lines);
    my $index = 0;

    # Muliple lines!!!
    if (scalar(@lines) > 1) {
        for (my $i = 0; $i < scalar(@lines); $i++){
            $index = $i if $lines[$i] =~ /protein_coding/ig;
        }
    }

    my @fields = split("\t", $lines[$index]);

    $distance = $fields[8] if $fields[8] =~ /\d+/;

    my @gene_features = split(";",$fields[7]);

    # Capturing Ensembl id and gene name.
    $geneID    = $1 if $gene_features[0] =~ /\"(.+)\"/;
    $gene_name = $1 if $gene_features[4] =~ /\"(.+)\"/;

    # Return stuff:
    return ($gene_name, $geneID, $distance);
}

1;
