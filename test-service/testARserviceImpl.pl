use Bio::KBase::atomic_regulons::atomic_regulonsImpl;
use strict;
use Data::Dumper;
use Bio::KBase::CDMI::CDMIClient;

my $BS = 0; # flag to turn on processing of B. subtilis data

my $genome_id = "kb|g.20403"; # for MP
$genome_id = "kb|g.423" if $BS;

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

open(MP,"mp.table") if ! $BS; # for MP
open(MP,"bs.table.43") if $BS;

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

my %comp;
if ($BS) {
    open (CP, "bs_comp.txt");
    while (<CP>) {
	chomp;
	my ($bs113, $bs1, $bs43, $bs49, $func) = split "\t";
	$comp{$bs43} = $bs113; # change this for other versions of BS
    }
    close CP;
}

open(AR,">ar.out");
foreach my $ar (sort { $a->{"ar_id"} <=> $b->{"ar_id"} } @$atomic_regulons) {
    my $ar_id = $ar->{"ar_id"};
    my @pegs;

    foreach my $fid (@{$ar->{"feature_ids"}}) {
	push @pegs, $BS ? $comp{$fid2peg{$fid}} : $fid2peg{$fid};
    }

    my $pegs = join ",", @pegs;
    print AR $ar_id, "\t", $pegs, "\n";
}
close(AR);

open(ARC, ">ar_calls.out");
foreach my $ac (sort { $a->{"sample_name"} cmp $b->{"sample_name"} || $a->{"ar_id"} <=> $b->{"ar_id"} } @$ar_calls) {
    print ARC $ac->{"sample_name"}, "\t", $ac->{"ar_id"}, "\t", $ac->{"on_off_unknown"}, "\n";
}
close(ARC);

open(FRC, ">feature_calls.out");
foreach my $fc (sort { $a->{"sample_name"} cmp $b->{"sample_name"} || $a->{"feature_id"} cmp $b->{"feature_id"} } @$feature_calls) {
    my $peg = $BS ? $comp{$fid2peg{$fc->{"feature_id"}}} : $fid2peg{$fc->{"feature_id"}};
    print FRC $fc->{"sample_name"}, "\t", $peg, "\t", $fc->{"on_off_unknown"}, "\n";
}
close(FRC);

