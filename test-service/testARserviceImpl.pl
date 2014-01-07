use Bio::KBase::atomic_regulons::atomic_regulonsImpl;
use strict;
use Data::Dumper;
use Bio::KBase::CDMI::CDMIClient;

#my $genome_id = "kb|g.20403"; # for MP
my $genome_id = "kb|g.423"; # for BS
my $csO = Bio::KBase::CDMI::CDMIClient->new_for_script();
my $featuresH = $csO->genomes_to_fids([$genome_id],["CDS","rna"]);
my $sourceIdsH = $csO->get_entity_Feature($featuresH->{$genome_id},["source_id"]);
my (%peg2fid, %fid2peg);

foreach my $fid (keys %$sourceIdsH) {
    my $peg = $sourceIdsH->{$fid}->{'source_id'};
    if (defined $peg) {
	$peg2fid{$peg} = $fid;
	$fid2peg{$fid} = $peg;
    }
}

my $table = {};

#open(MP,"mp.table"); # for MP
open(MP,"bs.table.43"); # for BS

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
my ($atomic_regulons, $feature_calls, $ar_calls) = $impl->compute_atomic_regulons($table,$genome_id);

foreach my $ar (@$atomic_regulons) {
    my $ar_id = $ar->{"ar_id"};
    my @pegs;

    foreach my $fid (@{$ar->{"feature_ids"}}) {
	push @pegs, $fid2peg{$fid};
    }

    my $pegs = join ",", @pegs;
    print $ar_id, "\t", $pegs, "\n";
}
