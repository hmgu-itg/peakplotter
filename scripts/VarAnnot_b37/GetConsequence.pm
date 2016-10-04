package GetConsequence;


use strict;
use warnings;
use HTTP::Tiny;
use JSON;
use Data::Dumper;


# This routine downloads a list of overlapping genes.
# input:
    # Usual format of the input lines.
    #    hash{line#}->{
    #            "chr"   => #,
    #            "start" => #,
    #            "end"   => #,
    #            "alt"   => #
    #     }
sub Consequence {
    print STDERR "[Info] Retieving variant consequences.";

    my %variants = %{$_[0]};
    my @Fields   = @{$_[1]};

    my @ConseqenceFields = ("Variant_Consequence","Variant_Impact", "Variant_Polyphen",
                            "Variant_PolyphenScore", "Variant_SiftScore", "Variant_Codons",
                            "Variant_AminoAcids", "Variant_ProteinPosition");
    push(@Fields, @ConseqenceFields);

    foreach my $variant (keys %variants){

        print STDERR ".";

        ## Initializing values to be returned:
        $variants{$variant}->{"Variant_Consequence"}     = "-";
        $variants{$variant}->{"Variant_SiftScore"}       = "-";
        $variants{$variant}->{"Variant_PolyphenScore"}   = "-";
        $variants{$variant}->{"Variant_Impact"}          = "-";
        $variants{$variant}->{"Variant_Codons"}          = "-";
        $variants{$variant}->{"Variant_AminoAcids"}      = "-";
        $variants{$variant}->{"Variant_ProteinPosition"} = "-";
        $variants{$variant}->{"Variant_Polyphen"}        = "-";

        # Creating a query URL string:
        my $query_string = sprintf("%s:%s-%s:1/%s", $variants{$variant}->{"chr"}, $variants{$variant}->{"start"}, $variants{$variant}->{"end"}, $variants{$variant}->{"alt"});
        $query_string = sprintf("%s:%s-%s:1/%s", $variants{$variant}->{"chr"}, $variants{$variant}->{"end"}, $variants{$variant}->{"start"}, $variants{$variant}->{"alt"}) if $variants{$variant}->{"start"} > $variants{$variant}->{"end"};

        my $URL = "http://grch37.rest.ensembl.org/vep/human/region/$query_string";
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

        my $consequences = decode_json($response->{content});

        $variants{$variant}->{"Variant_Consequence"} = $consequences->[0]->{"most_severe_consequence"};
        #print "Most severe consequence: ",$consequences->{"most_severe_consequence"},"\n";
        # If there are overlapping gene at this point I don't care just with only one.
        # I know it is an issue to address.

        # Check all frequency data, but filtering for the ones we are especially interested in:
        foreach my $transcript (@{$consequences->[0]->{"transcript_consequences"}}) {

            next if $transcript->{'transcript_id'} ne $variants{$variant}->{"Transcript_Ensembl_ID"};

            # From the appropriate element, we read out the data:
            $variants{$variant}->{"Variant_SiftScore"}       = $transcript->{"sift_score"}          || "-";
            $variants{$variant}->{"Variant_PolyphenScore"}   = $transcript->{"polyphen_score"}      || "-";
            $variants{$variant}->{"Variant_Polyphen"}        = $transcript->{"polyphen_prediction"} || "-";
            $variants{$variant}->{"Variant_Impact"}          = $transcript->{"impact"}              || "-";
            $variants{$variant}->{"Variant_Codons"}          = $transcript->{"codons"}              || "-";
            $variants{$variant}->{"Variant_AminoAcids"}      = $transcript->{"amino_acids"}         || "-";
            $variants{$variant}->{"Variant_ProteinPosition"} = $transcript->{"protein_end"}         || "-";
        }
    }

    print STDERR " done.\n";
    return (\%variants,\@Fields);

}

1;