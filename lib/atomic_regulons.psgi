use Bio::KBase::atomic_regulons::atomic_regulonsImpl;

use Bio::KBase::atomic_regulons::Service;
use Plack::Middleware::CrossOrigin;



my @dispatch;

{
    my $obj = Bio::KBase::atomic_regulons::atomic_regulonsImpl->new;
    push(@dispatch, 'atomic_regulons' => $obj);
}


my $server = Bio::KBase::atomic_regulons::Service->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub { $server->handle_input(@_) };

$handler = Plack::Middleware::CrossOrigin->wrap( $handler, origins => "*", headers => "*");
