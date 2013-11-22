package Bio::KBase::atomic_regulons::atomic_regulonsImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

atomic_regulons

=head1 DESCRIPTION



=cut

#BEGIN_HEADER
use SeedUtils;
use gjoseqlib;
use Data::Dumper;
use Bio::KBase::CDMI::CDMIClient;
use Bio::KBase::Utilities::ScriptThing;
use File::Temp qw /tempdir/;

#-----------------------------------------------------------------------------
#  $cc = correl_coef( \@x, \@y ) ; copied from gjostat
#-----------------------------------------------------------------------------
sub correl_coef
{   
    my ($xref, $yref) = @_;
    (@$xref > 2) || return undef;
    my (@x) = @$xref;
    my (@y) = @$yref;
    my $n = @x;

    my ($xsum, $x2sum, $ysum, $y2sum, $xysum) = (0) x 5;
    my ($i, $xi, $yi);

    for ($i = 1; $i <= $n; $i++)
    {
	$xi = shift @x; $xsum += $xi; $x2sum += $xi*$xi;
	$yi = shift @y; $ysum += $yi; $y2sum += $yi*$yi;
	$xysum += $xi*$yi;
    }

    my $xsd = sqrt( ($x2sum - ($xsum*$xsum/$n)) / ($n - 1) );
    my $ysd = sqrt( ($y2sum - ($ysum*$ysum/$n)) / ($n - 1) );
    if (($xsd == 0) || ($ysd == 0)) { return undef }
    ( $xysum - $xsum * $ysum / $n ) / ( $xsd * $ysd * ( $n - 1 ) );
}

# subs from ex_get_adjacency_based_estimates.pl

sub possible_clusters {
    my($pegs_with_locs,$values) = @_;

    my $clusters = [];
    &gather_from_strand($pegs_with_locs,$values,$clusters);
    my $flipped = &flip($pegs_with_locs);
    my @pegs_with_locR = sort { ($a->[1]->[0] cmp $b->[1]->[0]) or
				    (($a->[1]->[1]+$a->[1]->[2]) <=> ($b->[1]->[1]+$b->[1]->[2]))
                              }
                         @$flipped;
    &gather_from_strand(\@pegs_with_locR,$values,$clusters);
    return $clusters;
}

sub gather_from_strand {
    my($pegs_with_locs,$values,$clusters) = @_;

    my $max = -1;
    my $i;
    while (($i = &next_to_try($max,$pegs_with_locs)) < (@$pegs_with_locs - 1))
    {
#	print STDERR "start grouping ",&Dumper($pegs_with_locs->[$i]);
	my $peg = $pegs_with_locs->[$i]->[0];
	my $cluster = [$peg];
	$max = $i;
	my $j;
	for ($j=$i+1; 
	     ($j < @$pegs_with_locs) && 
	     &ok_in_runF($pegs_with_locs->[$j-1],$pegs_with_locs->[$j],$values,$cluster); 
	     $j++)
	{
#	    print STDERR "adding ",&Dumper($pegs_with_locs->[$j]);
	    push(@$cluster,$pegs_with_locs->[$j]->[0]);
	    $max = $j;
	}
#	print STDERR "before going left ",&Dumper($cluster);
	if ((($i-1) >= 0) && 
	    (&cc_corr($pegs_with_locs->[$i-1]->[0],$cluster,$values) >= 0.4) &&
	    ($pegs_with_locs->[$i-1]->[1]->[0] eq $pegs_with_locs->[$i]->[1]->[0]) &&
	    ($pegs_with_locs->[$i-1]->[1]->[1] > $pegs_with_locs->[$i-1]->[1]->[2]))
	{
#	    print STDERR "going left with ",&Dumper($pegs_with_locs->[$i-1]);
	    push(@$cluster,$pegs_with_locs->[$i-1]->[0]);
	    for ($j=$i-2; 
	         ($j >= 0) && &ok_in_runB($pegs_with_locs->[$j+1],$pegs_with_locs->[$j],$values,$cluster); 
		 $j--)
	    {
#		print STDERR "adding $pegs_with_locs->[$j]->[0]\n";
		push(@$cluster,$pegs_with_locs->[$j]->[0]);
	    }
	}
#	print STDERR "final cluster ",&Dumper($cluster); 
	if (@$cluster > 1) 
	{ 
	    push(@$clusters,$cluster);
#	    print STDERR "keeping ",join(",",@$cluster),"\n";
	}
    }
}

sub ok_in_runB {
    my($x,$y,$values,$cluster) = @_;

    my $loc1 = $x->[1];
    my($c1,$b1,$e1) = @$loc1;
    my $loc2 = $y->[1];
    my($c2,$b2,$e2) = @$loc2;
    if ($c1 ne $c2) { return 0 }
    if ($b2 < $e2)  { return 0 }
    if (! &cc_corr($y->[0],$cluster,$values)) { return 0 }
    return (abs($e1-$b2) < 200);
}

sub ok_in_runF {
    my($x,$y,$values,$cluster) = @_;

    my $loc1 = $x->[1];
    my($c1,$b1,$e1) = @$loc1;
    my $loc2 = $y->[1];
    my($c2,$b2,$e2) = @$loc2;
    if ($c1 ne $c2) { return 0 }
    if ($b2 > $e2)  { return 0 }
    if (! &cc_corr($y->[0],$cluster,$values)) { return 0 }
    return (abs($b2-$e1) < 200);
}

sub cc_corr {
    my($peg1,$cluster,$values) = @_;

    my $sum = 0;
    foreach my $peg2 (@$cluster)
    {
	my $v = correl_coef($values->{$peg1}, $values->{$peg2});
	if ((! defined($v)) || ($v < 0.4)) { return 0 }
	$sum += $v;
    }
    return (($sum / @$cluster) >= 0.6);
}

sub next_to_try {
    my($max,$pegs_with_locs) = @_;
    my $i;
    for ($i=$max+1; ($i < @$pegs_with_locs) && 
	          ($pegs_with_locs->[$i]->[1]->[1] > $pegs_with_locs->[$i]->[1]->[2]);
	 $i++) {}
    return $i;
}

sub flip {
    my($pegs_with_locs) = @_;

    my $flipped = [];
    foreach my $x (@$pegs_with_locs)
    {
	my($peg,$loc) = @$x;
	my($contig,$beg,$end) = @$loc;
	my $loc1 = [$contig,100000000-$beg,100000000-$end];
	push(@$flipped,[$peg,$loc1]);
    }
    return $flipped;
}

# subs from ex_get_subsystem_based_estimates.pl

sub split_on_pc {
    my($pegs,$values) = @_;

    my @sets = ();
    my %used;
    my $i;
    for ($i=0; ($i < (@$pegs - 1)); $i++)
    {
	if (! $used{$pegs->[$i]})
	{
	    my @poss = ($pegs->[$i]);
	    my $j;
	    for ($j=$i+1; ($j < @$pegs); $j++)
	    {
		if (&ss_corr($pegs->[$j],\@poss,$values))
		{
		    push(@poss,$pegs->[$j]);
		    $used{$pegs->[$j]} = 1;
		}
	    }
	    push(@sets,\@poss);
	}
    }
    return @sets;
}

sub ss_corr {
    my($peg1,$cluster,$values) = @_;

    my $sum = 0;
    foreach my $peg2 (@$cluster)
    {
	my $v = correl_coef($values->{$peg1}, $values->{$peg2});
	if ((! defined($v)) || ($v < 0.4)) { return 0 }
	$sum += $v;
    }
    return (($sum / @$cluster) >= 0.7);
}

# subs from ex_make_initial_atomic_regulons.pl

sub to_pvec {
    my($exps,$pegs,$exp_peg_call) = @_;

    my $pvecs = {};
    foreach my $peg (@$pegs)
    {
	foreach my $exp (@$exps)
	{
	    push(@{$pvecs->{$peg}},$exp_peg_call->{$exp}->{$peg});
	}
    }
    return $pvecs;
}

sub same_profile {
    my($pvecsI,$pvecsJ) = @_;

    my $i;
    my $solid1_on = 0;
    my $solid1_off = 0;
    my $solid2_off = 0;
    my $solid2_on = 0;

    my $ok = 1;
    for ($i=0; ($i < @$pvecsI) && $ok; $i++) 
    {	
	my $v1 = $pvecsI->[$i];
	my $v2 = $pvecsJ->[$i];
	$ok = &ok($v1,$v2);
	if    ($v1 == 1)    { $solid1_on++ }
	elsif ($v1 == -1)   { $solid1_off++ }
	if    ($v2 == 1)    { $solid2_on++ }
	elsif ($v2 == -1)   { $solid2_off++ }
    }
    return (($i == @$pvecsI) && 
	    ($solid1_on  > (@$pvecsI * 0.2)) &&
	    ($solid1_off > (@$pvecsI * 0.2)) &&
	    ($solid2_on  > (@$pvecsI * 0.2)) &&
	    ($solid2_off > (@$pvecsJ * 0.2)));
}

sub ok {
    my($x,$y) = @_;

    return (($x == $y) || ($x == 0) || ($y == 0));
}

sub cluster_objects {
    my %to_cluster;
    my %in_cluster;

    my $nxt = 1;
    while (defined(my $input = shift @_))
    {
	foreach my $set (@$input) {
	    for (my $i = 0; $i < @$set - 1; $i++)
	    {
		my $obj1 = $set->[$i];
		my $obj2 = $set->[$i+1];
		my $in1 = $to_cluster{$obj1};
		my $in2 = $to_cluster{$obj2};

		if (defined($in1) && defined($in2) && ($in1 != $in2))
		{
		    push(@{$in_cluster{$in1}},@{$in_cluster{$in2}});
		    foreach $_ (@{$in_cluster{$in2}})
		    {
			$to_cluster{$_} = $in1;
		    }
		    delete $in_cluster{$in2};
		}
		elsif ((! defined($in1)) && defined($in2))
		{
		    push(@{$in_cluster{$in2}},$obj1);
		    $to_cluster{$obj1} = $in2;
		}
		elsif ((! defined($in2)) && defined($in1))
		{
		    push(@{$in_cluster{$in1}},$obj2);
		    $to_cluster{$obj2} = $in1;
		}
		elsif ((! defined($in1)) && (! defined($in2)))
		{   
		    $to_cluster{$obj1} = $to_cluster{$obj2} = $nxt;
		    $in_cluster{$nxt} = [$obj1,$obj2];
		    $nxt++;
		}
	    }
	}
    }

    return \%in_cluster;
}

#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR
mkdir("/mnt/tmp");
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 compute_atomic_regulons

  $atomic_regulons, $feature_calls, $ar_calls = $obj->compute_atomic_regulons($expression_values, $genome_id)

=over 4

=item Parameter and return types

=begin html

<pre>
$expression_values is an ExpressionValues
$genome_id is a string
$atomic_regulons is a reference to a list where each element is an AtomicRegulon
$feature_calls is a reference to a list where each element is a FeatureOnOffCall
$ar_calls is a reference to a list where each element is an AtomicRegulonOnOffCall
ExpressionValues is a reference to a hash where the following keys are defined:
	sample_names has a value which is a reference to a list where each element is a string
	expression_vectors has a value which is a reference to a hash where the key is a string and the value is a reference to a list where each element is a string
AtomicRegulon is a reference to a hash where the following keys are defined:
	ar_id has a value which is a string
	feature_ids has a value which is a reference to a list where each element is a string
FeatureOnOffCall is a reference to a hash where the following keys are defined:
	sample_name has a value which is a string
	feature_id has a value which is a string
	on_off_unknown has a value which is an int
AtomicRegulonOnOffCall is a reference to a hash where the following keys are defined:
	sample_name has a value which is a string
	ar_id has a value which is a string
	on_off_unknown has a value which is an int

</pre>

=end html

=begin text

$expression_values is an ExpressionValues
$genome_id is a string
$atomic_regulons is a reference to a list where each element is an AtomicRegulon
$feature_calls is a reference to a list where each element is a FeatureOnOffCall
$ar_calls is a reference to a list where each element is an AtomicRegulonOnOffCall
ExpressionValues is a reference to a hash where the following keys are defined:
	sample_names has a value which is a reference to a list where each element is a string
	expression_vectors has a value which is a reference to a hash where the key is a string and the value is a reference to a list where each element is a string
AtomicRegulon is a reference to a hash where the following keys are defined:
	ar_id has a value which is a string
	feature_ids has a value which is a reference to a list where each element is a string
FeatureOnOffCall is a reference to a hash where the following keys are defined:
	sample_name has a value which is a string
	feature_id has a value which is a string
	on_off_unknown has a value which is an int
AtomicRegulonOnOffCall is a reference to a hash where the following keys are defined:
	sample_name has a value which is a string
	ar_id has a value which is a string
	on_off_unknown has a value which is an int


=end text



=item Description

input a list of expression values for a genome, compute atomic regulons

=back

=cut

sub compute_atomic_regulons
{
    my $self = shift;
    my($expression_values, $genome_id) = @_;

    my @_bad_arguments;
    (ref($expression_values) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"expression_values\" (value was \"$expression_values\")");
    (!ref($genome_id)) or push(@_bad_arguments, "Invalid type for argument \"genome_id\" (value was \"$genome_id\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to compute_atomic_regulons:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'compute_atomic_regulons');
    }

    my $ctx = $Bio::KBase::atomic_regulons::Service::CallContext;
    my($atomic_regulons, $feature_calls, $ar_calls);
    #BEGIN compute_atomic_regulons

    # initialize return values
    $atomic_regulons = [];
    $feature_calls = [];
    $ar_calls = [];

    # create a temporary directory for intermediate results
#    my $dataD = tempdir ( DIR => "/mnt/tmp", CLEANUP => 0 );

    # get the features and roles for the genome
    my $csO = Bio::KBase::CDMI::CDMIClient->new_for_script();
    my $featuresH = $csO->genomes_to_fids([$genome_id],["CDS","rna"]);
    my $rolesH = $csO->fids_to_roles($featuresH->{$genome_id});

    # compute chromosomal clusters 
    my @fids_with_values;
    foreach my $fid (@{$featuresH->{$genome_id}}) {
	if (exists $expression_values->{"expression_vectors"}->{$fid}) {
	    push @fids_with_values, $fid;
	} else {
	    print STDERR "No expression values for $fid\n";
	}
    }
    my $locH    = $csO->fids_to_locations(\@fids_with_values);
    my @fids_with_loc = sort { ($a->[1]->[0] cmp $b->[1]->[0]) or
				   (($a->[1]->[1]+$a->[1]->[2]) <=> ($b->[1]->[1]+$b->[1]->[2]))}
                        map  { my($contig, $pos, $strand, $length) = @{$locH->{$_}->[0]};
		           [$_,[$contig,$pos,($strand eq '+') ? ($pos+($length-1)) : ($pos -($length-1))]] }
                        keys(%$locH);      
    my $chromosomal_clusters = &possible_clusters(\@fids_with_loc,$expression_values->{"expression_vectors"});

    # compute subsystem clusters
    my $genomeH = $csO->genomes_to_subsystems( [$genome_id] );
    my @subs = sort map { ($_->[0] =~ /^\*?(0|-1)$/) ? () : $_->[1] } @{$genomeH->{$genome_id}};

    my $ss_clusters = [];
    my %bad;

    foreach my $sub (@subs)
    {
	my %pegs;
	my @pegs;

	my $subH = $csO->subsystems_to_fids( [$sub], [$genome_id]);
	my $sub_entry = $subH->{$sub}->{$genome_id};
	@pegs = ();
	foreach my $pair (@$sub_entry)
	{
	    my $pegs = $pair->[1];
	    foreach my $peg (@$pegs) { 
		if (exists $expression_values->{"expression_vectors"}->{$peg}) {
		    $pegs{$peg} = 1;
		} else {
		    print STDERR "No expression values for $peg\n";
		}
	    }
	}
	@pegs = sort { &SeedUtils::by_fig_id($a,$b) } keys(%pegs);

	my @sets = grep { @$_ > 1 } split_on_pc(\@pegs,$expression_values->{"expression_vectors"});

	if (@sets > ((@pegs + 2) / 3))
	{
	    $bad{$sub} = 1;
#	print STDERR &Dumper([$sub,\@sets]);
	}
	else
	{
	    foreach my $set (@sets)
	    {
		if (@$set > 1)
		{
		    push @$ss_clusters, $set;
		}
	    }
	}
    }

    foreach $_ (keys(%bad))
    {
	print STDERR "bad subsystem\t$_\n";
    }

    # code from ex_make_on_off_calls.pl
    
    my $heredoc = <<END;
Alanyl-tRNA synthetase (EC 6.1.1.7)
Arginyl-tRNA synthetase (EC 6.1.1.19)
Asparaginyl-tRNA synthetase (EC 6.1.1.22)
Aspartyl-tRNA synthetase (EC 6.1.1.12)
Cysteinyl-tRNA synthetase (EC 6.1.1.16)
DNA-directed RNA polymerase alpha subunit (EC 2.7.7.6)
DNA-directed RNA polymerase beta subunit (EC 2.7.7.6)
DNA-directed RNA polymerase beta\' subunit (EC 2.7.7.6)
DNA-directed RNA polymerase omega subunit (EC 2.7.7.6)
Glutaminyl-tRNA synthetase (EC 6.1.1.18)
Glutamyl-tRNA synthetase (EC 6.1.1.17)
Glycyl-tRNA synthetase alpha chain (EC 6.1.1.14)
Glycyl-tRNA synthetase beta chain (EC 6.1.1.14)
Histidyl-tRNA synthetase (EC 6.1.1.21)
Isoleucyl-tRNA synthetase (EC 6.1.1.5)
Leucyl-tRNA synthetase (EC 6.1.1.4)
LSU ribosomal protein L10p (P0)
LSU ribosomal protein L11p (L12e)
LSU ribosomal protein L13p (L13Ae)
LSU ribosomal protein L14p (L23e)
LSU ribosomal protein L15p (L27Ae)
LSU ribosomal protein L16p (L10e)
LSU ribosomal protein L17p
LSU ribosomal protein L18p (L5e)
LSU ribosomal protein L19p
LSU ribosomal protein L1p (L10Ae)
LSU ribosomal protein L20p
LSU ribosomal protein L21p
LSU ribosomal protein L22p (L17e)
LSU ribosomal protein L23p (L23Ae)
LSU ribosomal protein L24p (L26e)
LSU ribosomal protein L25p
LSU ribosomal protein L27p
LSU ribosomal protein L28p
LSU ribosomal protein L29p (L35e)
LSU ribosomal protein L2p (L8e)
LSU ribosomal protein L30p (L7e)
LSU ribosomal protein L31p
LSU ribosomal protein L32p
LSU ribosomal protein L33p
LSU ribosomal protein L34p
LSU ribosomal protein L35p
LSU ribosomal protein L36p
LSU ribosomal protein L3p (L3e)
LSU ribosomal protein L4p (L1e)
LSU ribosomal protein L5p (L11e)
LSU ribosomal protein L6p (L9e)
LSU ribosomal protein L9p
Lysyl-tRNA synthetase (class II) (EC 6.1.1.6)
Methionyl-tRNA synthetase (EC 6.1.1.10)
Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)
Phenylalanyl-tRNA synthetase beta chain (EC 6.1.1.20)
Seryl-tRNA synthetase (EC 6.1.1.11)
SSU ribosomal protein S10p (S20e)
SSU ribosomal protein S11p (S14e)
SSU ribosomal protein S12p (S23e)
SSU ribosomal protein S13p (S18e)
SSU ribosomal protein S15p (S13e)
SSU ribosomal protein S16p
SSU ribosomal protein S17p (S11e)
SSU ribosomal protein S19p (S15e)
SSU ribosomal protein S1p
SSU ribosomal protein S20p
SSU ribosomal protein S21p
SSU ribosomal protein S2p (SAe)
SSU ribosomal protein S3p (S3e)
SSU ribosomal protein S4p (S9e)
SSU ribosomal protein S5p (S2e)
SSU ribosomal protein S6p
SSU ribosomal protein S7p (S5e)
SSU ribosomal protein S8p (S15Ae)
SSU ribosomal protein S9p (S16e)
Threonyl-tRNA synthetase (EC 6.1.1.3)
Tryptophanyl-tRNA synthetase (EC 6.1.1.2) ## proteobacterial type
Valyl-tRNA synthetase (EC 6.1.1.9)
END
    my %roles_always_on = map { $_ => 1 } split "\n", $heredoc;

    my %pegs_always_on;

    foreach my $fid (keys %$rolesH) {
	foreach my $role (@{$rolesH->{$fid}}) {
	    if (exists $roles_always_on{$role}) {
		$pegs_always_on{$fid} = 1;
		last;
	    }
	}
    }

    my %on_for_exp;
    my %off_for_exp;
    my %diff_by_exp;

    my $locusH = $expression_values->{"expression_vectors"};
    my $exp_counter = 0;

    foreach my $exp (@{$expression_values->{"sample_names"}})
    {
	my @on_values = sort { $a <=> $b } 
	                map { $pegs_always_on{$_} ? $locusH->{$_}->[$exp_counter] : () } keys(%$locusH);
	my $I_10 = int(@on_values * 0.1);
	my $on_threshold = $on_values[$I_10];
	$on_for_exp{$exp} = $on_threshold;
	my @not_on_values = sort { $a <=> $b } 
            	            map { ((! $pegs_always_on{$_} ) && ($locusH->{$_}->[$exp_counter] < $on_threshold)) ? $locusH->{$_}->[$exp_counter] : () } 
	                    keys(%$locusH);
	my $I_80  = int(@not_on_values * 0.8);
	my $off_threshold = $not_on_values[$I_80];

	$off_for_exp{$exp} = $off_threshold;
	$diff_by_exp{$exp} = $on_threshold - $off_threshold;

	$exp_counter++;
    }
    my @diffs = sort { $a <=> $b } map { $diff_by_exp{$_} } keys(%diff_by_exp);
    my $I_25  = int(@diffs * 0.25);
    my $diff_threshold = $diffs[$I_25];
    foreach my $exp (keys(%off_for_exp))
    {
	if (($on_for_exp{$exp} - $off_for_exp{$exp}) < $diff_threshold)
	{
	    $off_for_exp{$exp} = $on_for_exp{$exp} - $diff_threshold;
	}
    }

    my %peg_on_off_by_exp;

    $exp_counter = 0;

    foreach my $exp (@{$expression_values->{"sample_names"}})
    {
	my $on_threshold = $on_for_exp{$exp};
	my $off_threshold = $off_for_exp{$exp};

	foreach my $peg (keys(%$locusH))
	{
	    my $v = $locusH->{$peg}->[$exp_counter];
	    if ($v >= $on_threshold)
	    {
		$peg_on_off_by_exp{$exp}->{$peg} = 1;
	    }
	    elsif ($v <= $off_threshold)
	    {
		$peg_on_off_by_exp{$exp}->{$peg} = -1;
	    }
	    else
	    {
		$peg_on_off_by_exp{$exp}->{$peg} = 0;
	    }
	}

	$exp_counter++;
    }
    
    # code from ex_make_initial_atomic_regulons.pl

    my @exp = sort keys(%peg_on_off_by_exp);
    my @pegs = keys %{$expression_values->{"expression_vectors"}};
    my $pvecs = &to_pvec(\@exp,\@pegs,\%peg_on_off_by_exp);

    # compute fids with same on/off profiles
    my $same_profile_clusters = [];

    for (my $i=0; ($i < $#pegs); $i++)
    {
	my $pvecsI = $pvecs->{$pegs[$i]};

	for (my $j=$i+1; ($j < @pegs); $j++)
	{
	    my $pvecsJ = $pvecs->{$pegs[$j]};
	    if (&same_profile($pvecsI,$pvecsJ))
	    {
		push @$same_profile_clusters, [$pegs[$i],$pegs[$j]];
	    }
	}
    }

    my $merged_clusters = &cluster_objects($chromosomal_clusters, $ss_clusters, $same_profile_clusters);
    print Dumper($merged_clusters);

    #END compute_atomic_regulons
    my @_bad_returns;
    (ref($atomic_regulons) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"atomic_regulons\" (value was \"$atomic_regulons\")");
    (ref($feature_calls) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"feature_calls\" (value was \"$feature_calls\")");
    (ref($ar_calls) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"ar_calls\" (value was \"$ar_calls\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to compute_atomic_regulons:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'compute_atomic_regulons');
    }
    return($atomic_regulons, $feature_calls, $ar_calls);
}




=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}

=head1 TYPES



=head2 genome_id

=over 4



=item Description

genome id needs to be KBase interpretable


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 ExpressionValues

=over 4



=item Description

table of expression levels for features in samples


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
sample_names has a value which is a reference to a list where each element is a string
expression_vectors has a value which is a reference to a hash where the key is a string and the value is a reference to a list where each element is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
sample_names has a value which is a reference to a list where each element is a string
expression_vectors has a value which is a reference to a hash where the key is a string and the value is a reference to a list where each element is a string


=end text

=back



=head2 AtomicRegulon

=over 4



=item Description

atomic regulon has an id and a set of features


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
ar_id has a value which is a string
feature_ids has a value which is a reference to a list where each element is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
ar_id has a value which is a string
feature_ids has a value which is a reference to a list where each element is a string


=end text

=back



=head2 FeatureOnOffCall

=over 4



=item Description

on/off/unknown call for a feature in one sample


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
sample_name has a value which is a string
feature_id has a value which is a string
on_off_unknown has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
sample_name has a value which is a string
feature_id has a value which is a string
on_off_unknown has a value which is an int


=end text

=back



=head2 AtomicRegulonOnOffCall

=over 4



=item Description

on/off/unknown call for an AtomicRegulon in one sample


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
sample_name has a value which is a string
ar_id has a value which is a string
on_off_unknown has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
sample_name has a value which is a string
ar_id has a value which is a string
on_off_unknown has a value which is an int


=end text

=back



=cut

1;
