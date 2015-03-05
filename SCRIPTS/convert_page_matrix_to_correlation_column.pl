my $pagedir ;
BEGIN{
    if ((!$ENV{PAGEDIR}) || ($ENV{PAGEDIR} eq '')) {
	$pagedir="./" ;
	print "The PAGEDIR environment variable is not set. It is set to default.\n";
    }
    else{
	$pagedir = $ENV{PAGEDIR};
    }
}

use lib "$pagedir/SCRIPTS";
use lib "$pagedir/SCRIPTS/PostScript-Simple-0.07/lib";

my $programdir = $pagedir."/PROGRAMS" ;
my $scriptdir  = $pagedir."/SCRIPTS" ;

use Table;
use Sets;
use Getopt::Long;
use PostScript::Simple;
use AggloClust;
use strict;
use Data::Dumper ;

my $pvaluematrixfile = shift(@ARGV) ;

print "Reading matrix ... ";

my $ta = Table->new;

#
#  read in the matrix file
#
$ta->loadFile($pvaluematrixfile);

# get an 2D array
my $a_ref_M      = $ta->getArray();

# header
my $a_ref_H      = shift @$a_ref_M; shift @$a_ref_H;
print "Done.\n";

my @ref ;
foreach my $c (@$a_ref_H) {
    my ($i1, $i2) = $c =~ /\[(.+?)\ (.+?)\]/;    
    push(@ref, ($i1+$i2)/2) ;
}



my $correl ;
my $cats ;
for (my $i=0; $i<@$a_ref_M; $i++) {
    my $r  = $a_ref_M->[$i];
    my $go = shift @$r;
    $cats->[$i] = $go ;
    my @pv ;
    for (my $j=0; $j<@$r; $j++) {
	# get log pvalues for over and under rep
	my ($lpo,$lpu) = $r->[$j] =~ /^(.+?)\/(.+)$/;

	# get a single value, pos if over-rep, neg if under-rep
	my $lp = undef;
	if (abs($lpo) > abs($lpu)) {
	    $lp = -$lpo;
	} else {
	    $lp = $lpu;
	}

	push(@pv, $lp) ;
    }
    $correl->[$i] = Sets::pearson(\@ref, \@pv) ;
}

open O, "> $pvaluematrixfile.correl" ;
for (my $i=0; $i<@$a_ref_M; $i++) {
    print O $cats->[$i], "\t", $correl->[$i], "\n" ;
}
