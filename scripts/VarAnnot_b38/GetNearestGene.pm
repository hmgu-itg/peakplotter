package GetMAFs;

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


# This routine retrieves the closest gene to a variant
# If it does not overlap with a gene.
# input:
    # Usual format of the input lines.
    #    hash{line#}->{
    #            "chr"   => #,
    #            "start" => #,
    #            "end"   => #,

    # Will only be executed if the $hash{$line}->{gene} is empty.
    # (we are only interested in the nearest gene, if the variant does not overlap with a gene)
sub MAFs {
    print STDERR "[Info] Finding nearest gene.";

    my %variants = %{$_[0]};
    my @Fields   = @{$_[1]};

    my @ProteinFields = ("Closest_gene_name", "Closest_gene_id", "Closest_gene_distance", "Closest_gene_type", "Closest_gene_description");
    push(@Fields, @ProteinFields);

    # Initialize database connection.
    my $registry = 'Bio::EnsEMBL::Registry';

    $registry->load_registry_from_db(
        -host => 'ensembldb.ensembl.org',
        -user => 'anonymous'
    );

    my $variation_adaptor = $registry->get_adaptor(
        'human',	# species
        'variation',	# database
        'variation'	# object type
    );

    foreach my $variant (keys %variants){

        print STDERR ".";

        # Initializing values to be returned:
        $variants{$variant}->{"Variant_Freq_CEU"} = "-";
        $variants{$variant}->{"Variant_Freq_TSI"} = "-";
        $variants{$variant}->{"Variant_Freq_FIN"} = "-";
        $variants{$variant}->{"Variant_Freq_GBR"} = "-";
        $variants{$variant}->{"Variant_Freq_IBS"} = "-";

        # We don't care variants that don't have a valid rsID:
        next if $variants{$variant}->{"rsID"} eq "-";

        # Retrieving data from Ensembl
        my $variation = $variation_adaptor->fetch_by_name($variants{$variant}->{"rsID"});

        # This database might not contain all variations, as we are using an older version.
        next unless $variation;

        my $alleles = $variation->get_all_Alleles();

        # Set which allele we are interested in:
        my $refAllele = $variants{$variant}->{"ref"};

        # looping through all alleles, but filtering those populations we are interested in:
        foreach my $allele (@{$alleles}) {

            next unless (defined $allele->population);

            my $allele_string   = $allele->allele;
            my $frequency       = $allele->frequency || '-';

            my $pop_code = substr($allele->population->name, -3);

            # Checking for population and allele:
            if (($Eur_pops =~ $pop_code) and ($allele_string eq $refAllele)){
                $variants{$variant}->{"Variant_Freq_$pop_code"}   = 1 - $frequency;
            }
        }


    }

    print STDERR " done.\n";
    return (\%variants, \@Fields);

}

1;