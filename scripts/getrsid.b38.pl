#!/usr/bin/env perl
use strict;
use warnings;
 
use HTTP::Tiny;
use Data::Dumper;
 
my $http = HTTP::Tiny->new();
 
my $server = 'http://rest.ensembl.org';
if($ARGV[0] eq "b37") {
  #b37
  print STDERR "Requested Build 37\n";
  $server='http://grch37.rest.ensembl.org';
  $ARGV[0]=$ARGV[1];
}

my $ext = '/vep/homo_sapiens/region';
open(IN, $ARGV[0]) or die "no open $ARGV[0]";
$|=1;


my @lol=<IN>;
chomp(@lol);
print STDERR "Processing RS-ids for ", scalar @lol, " variants.\n";
my @VAR;
push @VAR, [ splice @lol, 0, 200 ] while @lol;
foreach my $ar_pt (@VAR){
	my @ar=@$ar_pt;
my $loll=join("\",\"", @ar);
$loll="\"".$loll."\"";
my $string=$loll;
$string = "{ \"variants\" : [$string] }";

#print "\n\n$string\n\n";

my $response = $http->request('POST', $server.$ext, {
  headers => { 
  	'Content-type' => 'application/json',
  	'Accept' => 'application/json'
  },
 content => $string
});


#print STDERR Dumper($response) unless $response->{success};
if(! $response->{success}){
print STDERR "Received String:\n$string\nend\n";
print STDERR Dumper($response);
die("Catastrophic failure\n");
} 
 
use JSON;
use Data::Dumper;
my $hash = decode_json($response->{content});
#print Dumper $hash;
foreach my $input (@$hash){
	my $id=defined($input->{'colocated_variants'}->[0]->{'id'})?$input->{'colocated_variants'}->[0]->{'id'}:"novel";
	my $gene=defined( $input->{'transcript_consequences'}->[0]->{'gene_symbol'} )?$input->{'transcript_consequences'}->[0]->{'gene_symbol'}:"none";
	my $consequence=defined($input->{'most_severe_consequence'})?$input->{'most_severe_consequence'}:"none";
	print join ("\t", $input->{'end'}, $id ,$gene,$consequence   ), "\n";
	print STDERR join(" ", $id, $input->{'end'}), "\n";
}

}
# die();

# if(length $response->{content}) {
#   my $hash = decode_json($response->{content});
#   local $Data::Dumper::Terse = 1;
#   local $Data::Dumper::Indent = 1;
#   print Dumper $hash->[0];
#   print "\n";
# }

