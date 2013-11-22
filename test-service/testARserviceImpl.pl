use Bio::KBase::atomic_regulons::atomic_regulonsImpl;
use strict;
use Data::Dumper;
use Bio::KBase::CDMI::CDMIClient;

my $genome_id = "kb|g.20403";
my $csO = Bio::KBase::CDMI::CDMIClient->new_for_script();
my $featuresH = $csO->genomes_to_fids([$genome_id],["CDS","rna"]);
my $sourceIdsH = $csO->get_entity_Feature($featuresH->{$genome_id},["source_id"]);
my %peg2fid;

foreach my $fid (keys %$sourceIdsH) {
    my $peg = $sourceIdsH->{$fid}->{'source_id'};
    $peg2fid{$peg} = $fid if defined $peg;
}

my $table = {};

open(MP,"mp.table");

my $line = <MP>;
chomp $line;
my @sample_names = split /\s/, $line;
$table->{"sample_names"} = \@sample_names;

while (chomp($line = <MP>)) {
    my ($peg, @values) = split /\s/, $line;
    next if ! exists $peg2fid{$peg};
    $table->{"expression_vectors"}->{$peg2fid{$peg}} = \@values;
}

close(MP);

my $impl = new Bio::KBase::atomic_regulons::atomic_regulonsImpl;
my $result = $impl->compute_atomic_regulons($table,$genome_id);
print &Dumper($result);
