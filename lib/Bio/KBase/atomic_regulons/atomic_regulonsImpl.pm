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
	    (&corr($pegs_with_locs->[$i-1]->[0],$cluster,$values) >= 0.4) &&
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
    if (! &corr($y->[0],$cluster,$values)) { return 0 }
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
    if (! &corr($y->[0],$cluster,$values)) { return 0 }
    return (abs($b2-$e1) < 200);
}

sub corr {
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
    my $featuresH = $csO->genomes_to_fids([$genome_id],["CDS"]);
    my $rolesH = $csO->fids_to_roles($featuresH->{$genome_id});

    # compute chromosomal clusters 
    my @fids_with_values;
    foreach my $fid (@{$featuresH->{$genome_id}}) {
	if (exists $expression_values->{"expression_vectors"}->{$fid}) {
	    push @fids_with_values, $fid;
	} else {
	    print STDERR "No values for $fid\n";
	}
    }
    my $locH    = $csO->fids_to_locations(\@fids_with_values);
    my @fids_with_loc = sort { ($a->[1]->[0] cmp $b->[1]->[0]) or
				   (($a->[1]->[1]+$a->[1]->[2]) <=> ($b->[1]->[1]+$b->[1]->[2]))}
                        map  { my($contig, $pos, $strand, $length) = @{$locH->{$_}->[0]};
		           [$_,[$contig,$pos,($strand eq '+') ? ($pos+($length-1)) : ($pos -($length-1))]] }
                        keys(%$locH);      
    my $clusters = &possible_clusters(\@fids_with_loc,$expression_values->{"expression_vectors"});

    foreach my $cluster (@$clusters)
    {
	if (@$cluster > 1)
	{
	    my @pegs = sort { &SeedUtils::by_fig_id($a,$b) } @$cluster;
#	    join(",",@pegs),"\tClusterOnChromosome:$pegs[0],$pegs[$#pegs]\n";
	}
    }

    # compute subsystem clusters
    my $genomeH = $sapO->genomes_to_subsystems( -ids => [$genome_id] );
    my @subs = map { ($_->[1] =~ /^\*?(0|-1)$/) ? () : $_->[0] } @{$genomeH->{$genome_id}};

    my $subH = $csO->ids_in_subsystems( -subsystems => \@subs,
					-genome     => $genome_id);
    my @subs = sort  keys %$subH;

    my %bad;

    foreach my $sub (@subs)
    {
	my %pegs;
	my @pegs;

	my $sub_entry = $subH->{$sub};
	@pegs = ();
	foreach my $role (keys(%$sub_entry))
	{
	    my $pegs = $sub_entry->{$role};
	    foreach $_ (@$pegs) { $pegs{$_} = 1 }
	}
	@pegs = sort { &SeedUtils::by_fig_id($a,$b) } keys(%pegs);
	
	my @sets = grep { @$_ > 1 } split_on_pc(\@pegs,$corrH);

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
		    print join(",",@$set),"\tInSubsystem:$sub\n";
		}
	    }
	}
    }
foreach $_ (keys(%bad))
{
    print STDERR "bad subsystem\t$_\n";
}

sub split_on_pc {
    my($pegs,$corrH) = @_;

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
		if (&corr($pegs->[$j],\@poss,$corrH))
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

sub corr {
    my($peg1,$cluster,$corrH) = @_;

    my $sum = 0;
    foreach my $peg2 (@$cluster)
    {
	my $v = $corrH->{$peg1}->{$peg2};
	if ((! defined($v)) || ($v < 0.4)) { return 0 }
	$sum += $v;
    }
    return (($sum / @$cluster) >= 0.7);
}
    


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
