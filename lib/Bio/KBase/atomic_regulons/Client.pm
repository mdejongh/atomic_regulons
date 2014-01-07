package Bio::KBase::atomic_regulons::Client;

use JSON::RPC::Client;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

Bio::KBase::atomic_regulons::Client

=head1 DESCRIPTION





=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => Bio::KBase::atomic_regulons::Client::RpcClient->new,
	url => $url,
    };


    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




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
	expression_vectors has a value which is a reference to a hash where the key is a string and the value is a reference to a list where each element is a float
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
	expression_vectors has a value which is a reference to a hash where the key is a string and the value is a reference to a list where each element is a float
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
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 2)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function compute_atomic_regulons (received $n, expecting 2)");
    }
    {
	my($expression_values, $genome_id) = @args;

	my @_bad_arguments;
        (ref($expression_values) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"expression_values\" (value was \"$expression_values\")");
        (!ref($genome_id)) or push(@_bad_arguments, "Invalid type for argument 2 \"genome_id\" (value was \"$genome_id\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to compute_atomic_regulons:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'compute_atomic_regulons');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "atomic_regulons.compute_atomic_regulons",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'compute_atomic_regulons',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method compute_atomic_regulons",
					    status_line => $self->{client}->status_line,
					    method_name => 'compute_atomic_regulons',
				       );
    }
}



=head2 compute_atomic_regulons_CDS

  $atomic_regulons, $feature_calls, $ar_calls = $obj->compute_atomic_regulons_CDS($genome_id)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_id is a string
$atomic_regulons is a reference to a list where each element is an AtomicRegulon
$feature_calls is a reference to a list where each element is a FeatureOnOffCall
$ar_calls is a reference to a list where each element is an AtomicRegulonOnOffCall
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

$genome_id is a string
$atomic_regulons is a reference to a list where each element is an AtomicRegulon
$feature_calls is a reference to a list where each element is a FeatureOnOffCall
$ar_calls is a reference to a list where each element is an AtomicRegulonOnOffCall
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

compute atomic regulons for a genome from expression values in the CDS

=back

=cut

sub compute_atomic_regulons_CDS
{
    my($self, @args) = @_;

# Authentication: none

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function compute_atomic_regulons_CDS (received $n, expecting 1)");
    }
    {
	my($genome_id) = @args;

	my @_bad_arguments;
        (!ref($genome_id)) or push(@_bad_arguments, "Invalid type for argument 1 \"genome_id\" (value was \"$genome_id\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to compute_atomic_regulons_CDS:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'compute_atomic_regulons_CDS');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "atomic_regulons.compute_atomic_regulons_CDS",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'compute_atomic_regulons_CDS',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method compute_atomic_regulons_CDS",
					    status_line => $self->{client}->status_line,
					    method_name => 'compute_atomic_regulons_CDS',
				       );
    }
}



sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, {
        method => "atomic_regulons.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'compute_atomic_regulons_CDS',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method compute_atomic_regulons_CDS",
            status_line => $self->{client}->status_line,
            method_name => 'compute_atomic_regulons_CDS',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for Bio::KBase::atomic_regulons::Client\n";
    }
    if ($sMajor == 0) {
        warn "Bio::KBase::atomic_regulons::Client version is $svr_version. API subject to change.\n";
    }
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

table of expression levels for features in samples;
the order of the expression_levels should match the order of the sample_names
and there should be no missing values


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
sample_names has a value which is a reference to a list where each element is a string
expression_vectors has a value which is a reference to a hash where the key is a string and the value is a reference to a list where each element is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
sample_names has a value which is a reference to a list where each element is a string
expression_vectors has a value which is a reference to a hash where the key is a string and the value is a reference to a list where each element is a float


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

on (1) / off (-1) / unknown (0) call for a feature in one sample


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

on (1) / off (-1) / unknown (0) call for an AtomicRegulon in one sample


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

package Bio::KBase::atomic_regulons::Client::RpcClient;
use base 'JSON::RPC::Client';

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $obj) = @_;
    my $result;

    if ($uri =~ /\?/) {
       $result = $self->_get($uri);
    }
    else {
        Carp::croak "not hashref." unless (ref $obj eq 'HASH');
        $result = $self->_post($uri, $obj);
    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
