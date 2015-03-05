use lib "$ENV{PAGEDIR}/SCRIPTS";

use Getopt::Long;
use Table;
use Sets;

GetOptions("expfile=s" => \$expfile,
	   "pathway=s" => \$pathway,
	   "bin=s"     => \$bin,
	   "species=s" => \$species);

die "Please specify --species=STR (e.g. human_go) \n" if !defined($species);

#
# load pathways
#
my $speciesdir  = "$ENV{PAGEDIR}/PAGE_DATA/ANNOTATIONS/$species";

# load annotation
my $goindexfile = "$speciesdir/$species\_index.txt";
my $gonamesfile = "$speciesdir/$species\_names.txt";
my $descfile    = "$speciesdir/$species\_genedesc.txt";

# get matching pathway
my $ta = Table->new;
$ta->loadFile($gonamesfile);
my $a_ref = $ta->getArray();
my %matchinggo = ();
foreach my $r (@$a_ref) {
  if (($r->[1] =~ /$pathway/) || ($r->[0] eq $pathway)) {
    $matchinggo{ $r->[0] } = 1;
  }
}

# get genes in the corresp pathway
$ta->loadFile($goindexfile);
$a_ref = $ta->getArray();

my %genes = ();
foreach my $r (@$a_ref) {
  my $g = shift @$r;
  
  foreach my $c (@$r) {

    if (defined($matchinggo{ $c })) {
      $genes{ $g } = 1; 
    }	

  }
  
}

# load $descfile
$ta->loadFile($descfile);
my $h_ref_desc = $ta->getIndex(0);


# now go thru the expfile 
$ta->loadFile($expfile);
$a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  if (defined($genes{$r->[0]})) {
    if (!defined($bin) || (defined($bin) && ($r->[1] == $bin))) {
      print "$r->[0]\t$h_ref_desc->{$r->[0]}->[1]\t$h_ref_desc->{$r->[0]}->[2]\n";
    }
  }
}
