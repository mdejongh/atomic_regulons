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
#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR
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
