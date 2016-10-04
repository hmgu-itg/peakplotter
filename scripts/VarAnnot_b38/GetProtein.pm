package GetProtein;

=head1 Todo list

1) The more important part is the problem of UniprotKB vs. TrEMBL identifiers.
At this point I only care of the fromer one, as that is the manually curated and
updated dataset. I can imagine the TrEmbl information could also be useful.
Sometimes. Unfortunately, there could be so many entries, I think it is just
wasting time to check them all.

2) The script is not ready to handle multiple UniprotKB identifiers belonging to
a single Gene. This should also be sorted out.

=cut

use strict;
use warnings;
use HTTP::Tiny;
use JSON;



# This routine downloads a list of overlapping genes.
# input:
    # Usual format of the input lines.
    #    hash{line#}->{
    #            "Gene_id"   => #,
sub Protein {
    print STDERR "[Info] Retrieving information of overlapping protein.";

    my %variants = %{$_[0]};
    my @Fields   = @{$_[1]};

    my @ProteinFields = ("Protein_Ensembl_ID", "Transcript_Ensembl_ID", "Uniprot_ID", "Uniprot_Name",
                         "Uniprot_Function", "Uniprot_Disease", "Uniprot_Subunit", "Uniprot_Phenotype",
                         "Uniprot_localization", "Uniprot_Tissue", "Uniprot_Development");
    push(@Fields, @ProteinFields);


    foreach my $variant (keys %variants){

        print STDERR ".";

        # Initializing values to be returned:
        $variants{$variant}->{"Protein_Ensembl_ID"}   = "-";
        $variants{$variant}->{"Transcript_Ensembl_ID"}   = "-";
        $variants{$variant}->{"Uniprot_ID"}           = "-";
        $variants{$variant}->{"Uniprot_Name"}         = "-";
        $variants{$variant}->{"Uniprot_Disease"}      = "-";
        $variants{$variant}->{"Uniprot_Function"}     = "-";
        $variants{$variant}->{"Uniprot_Subunit"}      = "-";
        $variants{$variant}->{"Uniprot_Phenotype"}    = "-";
        $variants{$variant}->{"Uniprot_localization"} = "-";
        $variants{$variant}->{"Uniprot_Tissue"}       = "-";
        $variants{$variant}->{"Uniprot_Development"}  = "-";

        # Ensembl stable gene id
        my $geneID = $variants{$variant}->{"Gene_id"};

        # In same cases the variant is not overlapping with gene, then we don't care about it anymore.
        next if $geneID eq "-";

        # Step 1. Get coding sequences overlapping with the given gene.
        my $URL = sprintf("http://rest.ensembl.org/overlap/id/%s?feature=CDS", $geneID);
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
            print STDERR "[Warning] Downloading data from the Ensembl server failed. URL: $URL\n";
            next;
        }

        my $response_parsed = decode_json($response->{content});

        # Looping through the the returned values:
        my %Protein_IDs = ();
        foreach my $transcript (@{$response_parsed}){
            my $source =  $transcript->{source};
            my $ProtID = $transcript->{id};
            $Protein_IDs{$source}{$ProtID}++;
        }

        # Collecting those protein IDs that will be tested for xref to UniprotKB
        my @ids = ();

        # Is there manually annotated Protein_IDs?
        if (exists $Protein_IDs{"ensembl_havana"}) {
            @ids = keys (%{$Protein_IDs{"ensembl_havana"}});
        }
        # If not, we collect all protein IDs.
        else {
            foreach my $source (keys %Protein_IDs){
                push (@ids, keys ( %{$Protein_IDs{$source}} ));
            }
        }

        # We don't go further with variants that does not overlap with protein coding gene:
        # (We explicitly have to do this other wise the next steps will fail.)
        next unless @ids;

        # We loop through all returned protein IDs, and check which of them has
        # Cross-reference to the UniprotKB database:
        my %Ensembl_Uniprot_IDs = ();
        my %Ensembl_Transcript_IDs = ();
        foreach my $ProtID (@ids){
            my ($uniprot_ID, $transcript) = &_getUniPortID($ProtID);
            $Ensembl_Uniprot_IDs{$uniprot_ID}++;
            $Ensembl_Transcript_IDs{$transcript}++;
        }


        # At this point we are looking only for the first Protein Id, but it has to be sorted out.
        $variants{$variant}->{"Protein_Ensembl_ID"}    = $ids[0];
        $variants{$variant}->{"Transcript_Ensembl_ID"} = (keys %Ensembl_Transcript_IDs)[0];
        $variants{$variant}->{"Uniprot_ID"}            = (keys %Ensembl_Uniprot_IDs)[0];

        # Retriving
        my %uniprot_data = %{&_UniProtData($variants{$variant}->{"Uniprot_ID"})};

        # Updating dataset
        foreach my $field_name (keys %uniprot_data){
            $variants{$variant}->{$field_name} = $uniprot_data{$field_name}
        }

    }

    print STDERR " done.\n";
    return (\%variants, \@Fields);

}


sub _getUniPortID {
    my $ProtID = $_[0];

    print STDERR ".";

    # As a first step we try to find the
    my $URL = sprintf("http://rest.ensembl.org/xrefs/id/%s", $ProtID);
    my $http = HTTP::Tiny->new();
    my $UNIPROT_id = "-";

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

    # Parsing out uniprot ID from the returned text:
    foreach my $crossref (@{$hash}){
        if ($crossref->{dbname} eq "Uniprot/SWISSPROT" ) {
            $UNIPROT_id = $crossref->{primary_id};
        }
    }

    print STDERR ".";

    # Now as we know what is the Protein identifier, we check the corresponding transcriptID.
    $URL = sprintf("http://rest.ensembl.org/overlap/translation/%s", $ProtID);
    $http = HTTP::Tiny->new();
    my $Transcript_id = "-";

    $response = $http->get($URL,{
        headers => { 'Content-type' => 'application/json' }
    });

    # This short loop was added to make the code more robuts. We don't want it to crash any time,
    # there is some problem with downloading the data.
    $fail_count = 0;
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

    $hash = decode_json($response->{content});

    # Parsing out transcript ID from the returned text:
    foreach my $transcripts (@{$hash}){
        if ($transcripts->{"seq_region_name"} eq $ProtID ) {
            $Transcript_id = $transcripts->{"Parent"};
        }
    }


    return ($UNIPROT_id, $Transcript_id);
}



sub _UniProtData {

    print STDERR ".";

    my $line = $_[0];

    # Initializing returned variable:
    my %hash = (
        "Uniprot_Name"          =>  "-",
        "Uniprot_Disease"       =>  "-",
        "Uniprot_Function"      =>  "-",
        "Uniprot_Entry"         =>  "-",
        "Uniprot_Subunit"       =>  "-",
        "Uniprot_Phenotype"     =>  "-",
        "Uniprot_localization"  =>  "-",
        "Uniprot_Tissue"        =>  "-",
        "Uniprot_Development"   =>  "-"
    );

    return \%hash if $line eq "-";

    my $URL = "http://www.uniprot.org/uniprot/?query=id:".$line."&columns=id%2Ccomment%28FUNCTION%29%2Ccomment%28SUBUNIT%29%2Ccomment%28DEVELOPMENTAL%20STAGE%29%2Ccomment%28TISSUE%20SPECIFICITY%29%2Ccomment%28CATALYTIC%20ACTIVITY%29%2Ccomment%28DISRUPTION%20PHENOTYPE%29%2Ccomment%28SUBCELLULAR%20LOCATION%29%2Ccomment%28DISEASE%29%2Centry%20name&format=tab";
    my $http = HTTP::Tiny->new();

    my $response = $http->get($URL);

    my @fields = ();

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
#http://grch37.rest.ensembl.org/overlap/region/human/19:8670147-8670147?feature=variation
#http://grch37.rest.ensembl.org/overlap/region/human/19:8670147-8670147?feature=variation;content-type=application/json

    my @response_lines = split("\n",$response->{content});
    @fields = split("\t", $response_lines[1]);

    %hash = (
        "Uniprot_Name"          => $fields[9] || "-",
        "Uniprot_Disease"       => $fields[8] || "-",
        "Uniprot_Function"      => $fields[1] || "-",
        "Uniprot_Entry"         => $fields[0] || "-",
        "Uniprot_Subunit"       => $fields[2] || "-",
        "Uniprot_Phenotype"     => $fields[6] || "-",
        "Uniprot_localization"  => $fields[7] || "-",
        "Uniprot_Tissue"        => $fields[4] || "-",
        "Uniprot_Development"   => $fields[3] || "-"
    );

    # Cleaning up the downloaded data. It solely involves removing cross-references to different databases.
    foreach my $feature (keys %hash){
        $hash{$feature} = $1 if $hash{$feature} =~ /(.+?) \{.+/;
    }

    return \%hash;
}


1;
