use Bio::KBase::atomic_regulons::atomic_regulonsImpl;
use strict;
use Data::Dumper;
use Bio::KBase::CDMI::CDMIClient;

my $genome_id = "kb|g.372"; # Shewanella
#my $genome_id = "kb|g.744"; # Vibrio
my $csO = Bio::KBase::CDMI::CDMIClient->new_for_script();
my $featuresH = $csO->genomes_to_fids([$genome_id],["CDS","rna"]);
my $rolesH = $csO->fids_to_roles($featuresH->{$genome_id});

my $impl = new Bio::KBase::atomic_regulons::atomic_regulonsImpl;
my ($atomic_regulons, $feature_calls, $ar_calls) = $impl->compute_atomic_regulons_CDS($genome_id);

open(AR,">ar.out");
foreach my $ar (sort { $a->{"ar_id"} <=> $b->{"ar_id"} } @$atomic_regulons) {
    foreach my $fid (@{$ar->{"feature_ids"}}) {
	print AR $ar->{"ar_id"}, "\t", $fid, "\t", (defined $rolesH->{$fid} ? (join ";", @{$rolesH->{$fid}}) : ""), "\n";
    }
}
close(AR);

open(ARC, ">ar_calls.out");
foreach my $ac (sort { $a->{"sample_name"} cmp $b->{"sample_name"} || $a->{"ar_id"} <=> $b->{"ar_id"} } @$ar_calls) {
    print ARC $ac->{"sample_name"}, "\t", $ac->{"ar_id"}, "\t", $ac->{"on_off_unknown"}, "\n";
}
close(ARC);

open(FRC, ">feature_calls.out");
foreach my $fc (sort { $a->{"sample_name"} cmp $b->{"sample_name"} || $a->{"feature_id"} cmp $b->{"feature_id"} } @$feature_calls) {
    print FRC $fc->{"sample_name"}, "\t", $fc->{"feature_id"}, "\t", $fc->{"on_off_unknown"}, "\n";
}
close(FRC);

