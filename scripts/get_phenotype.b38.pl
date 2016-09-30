#!/usr/local/bin/perl

use strict;
use warnings;
 
use HTTP::Tiny;
 
my $http = HTTP::Tiny->new();

my $server = 'http://rest.ensembl.org';
if($ARGV[0] eq "b37") {
  #b37
  print STDERR "Requested Build 37\n";
  $server='http://grch37.rest.ensembl.org';
  $ARGV[0]=$ARGV[1];
}

 
open(IN, $ARGV[0]) or die "no open $ARGV[0]";
my @lol=<IN>;
chomp(@lol);
my $loll=join("\",\"", @lol);
$loll="\"".$loll."\"";

print STDERR "\n\n===========================\nget_phenotype.pl : Received input $loll\n=========================\n\n";




my $ext = '/variation/homo_sapiens?phenotypes=1';
my $response = $http->request('POST', $server.$ext, {
  headers => { 
  	'Content-type' => 'application/json',
  	'Accept' => 'application/json'
  },
  content => '{ "ids" : ['.$loll.' ] }',
});
 
die "Failed!\n" unless $response->{success};
 
 use JSON;
use Data::Dumper;
my $hash = decode_json($response->{content});

foreach my $key (keys(%$hash)){
  my $name=$hash->{$key}->{'name'};
  my $phenotypes="";
  my $temp= $hash->{$key}->{'phenotypes'};
  my @phs=@$temp;
  foreach my $ph (@phs){
    $phenotypes=join(",", $phenotypes, $ph->{'trait'})
  }
  $phenotypes=$phenotypes eq ""?"none":$phenotypes;
  print join ("\t", $name ,$phenotypes   ), "\n";
}

exit();
use JSON;
use Data::Dumper;
if(length $response->{content}) {
  my $hash = decode_json($response->{content});
  local $Data::Dumper::Terse = 1;
  local $Data::Dumper::Indent = 1;
  print Dumper $hash;
  print "\n";
}
