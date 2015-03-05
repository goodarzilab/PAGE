BEGIN {
	if ( ( !$ENV{PAGEDIR} ) || ( $ENV{PAGEDIR} eq '' ) ) {
		print
"The PAGEDIR environment variable is not set. Please set it using export or setenv (see PAGE tutorial online).\n";
		exit;
	}
}

my $pagedir = $ENV{PAGEDIR};
use lib "$ENV{PAGEDIR}/SCRIPTS/PostScript-Simple-0.07/lib";
my $programdir = $pagedir . "/PROGRAMS";
my $scriptdir  = $pagedir . "/SCRIPTS";

use lib "$ENV{PAGEDIR}/SCRIPTS";

use strict;
use Sets;
use PBS;
use Table;
use Getopt::Long;
use AggloClust;
use Data::Dumper;
use PostScript::Simple;
use FileHandle;

if ( @ARGV == 0 ) {
	die
"Usage: perl prmg.pl --expfile=FILE --firefile=FILE --pagefile=FILE --goindexfile=FILE --max_p=F (--tu_enabled=0/1 --tu_file=FILE --printmode=0/1) --dorna=0/1\n";
}

my $expfile     = undef;
my $firefile    = undef;
my $pagefile    = undef;
my $goindexfile = undef;
my $gonamesfile = undef;
my $species     = undef;
my $gonamesfile = undef;
my $max_p       = 0.05;
my $submit      = 0;
my $printmode   = 0;
my $tu_enabled  = 0;
my $tu_file     = undef;
my $dorna       = 1;

my $holdjob0 = undef;
my $holdjob1 = undef;

GetOptions(
	'expfile=s'     => \$expfile,
	'firefile=s'    => \$firefile,
	'goindexfile=s' => \$goindexfile,
	'gonamesfile=s' => \$gonamesfile,
	'species=s'     => \$species,
	'max_p=s'       => \$max_p,
	'submit=s'      => \$submit,
	'tu_enabled=s'  => \$tu_enabled,
	'tu_file=s'     => \$tu_file,
	'printmode=s'   => \$printmode,
	'dorna=s'       => \$dorna,
	'holdjob0=s'    => \$holdjob0,
	'holdjob1=s'    => \$holdjob1
);

if ( !( -d "$pagedir/TEMP" ) ) {
	mkdir "$pagedir/TEMP";
}
if ( ( $expfile =~ /\*/ ) or ( $submit == 1 ) ) {
	my $files = Sets::getFiles($expfile);

	my $walltime = "20:00:00";
	my $platform = undef;

	foreach my $file (@$files) {
		my $f = Sets::filename($file);
		mkdir "$file\_PAGE/" if !( -d "$file\_PAGE/" );
		my $expfile_nodups_prmg = "$file\_PAGE/$f";

		my $pwd = `pwd`;
		$pwd =~ s/\n//;
		my $time = Sets::getNiceDateTime(1);

		my $pbs = PBS->new;
		$pbs->setPlatform($platform) if ( defined($platform) );
		$pbs->setWallTime($walltime);
		$pbs->addCmd("cd $pwd");

		$pbs->setScriptName("$expfile_nodups_prmg.script");

		$pbs->addCmd("date");
		$pbs->addCmd("echo \"Running PRMG\"");
		$pbs->addCmd("export PAGEDIR=$pagedir");

		my $firefile = $file . "_FIRE/DNA_RNA/$f.summary";
		my $pagefile = $file . "_PAGE/pvmatrix.txt";

		my $cmd =
"perl $pagedir/prmg.pl --expfile=$file --firefile=$firefile --pagefile=$pagefile --goindexfile=$goindexfile --max_p=$max_p --dorna=$dorna --species=$species";
		$pbs->addCmd($cmd);

		my $page_jobid;
		if ( $submit == 0 ) {
			$pbs->execute;
		}
		elsif ( $submit == 1 ) {

			if ( defined($holdjob0) && length( $holdjob0 > 0 ) ) {
				$pbs->addDepJob($holdjob0);
			}
			if ( defined($holdjob1) && length( $holdjob1 > 0 ) ) {
				$pbs->addDepJob($holdjob1);
			}

			$page_jobid = $pbs->submit;
			print "Submitted job $page_jobid.\n";
		}
	}
	exit(1);
}

if ( !defined($firefile) and $dorna == 1 ) {
	my $fn = Sets::filename($expfile);
	$firefile = "$expfile\_FIRE/DNA_RNA/$fn.summary";
}
if ( !defined($firefile) and $dorna == 0 ) {
	my $fn = Sets::filename($expfile);
	$firefile = "$expfile\_FIRE/DNA/$fn.summary";
}
if ( !defined($pagefile) ) {
	$pagefile = "$expfile\_PAGE/pvmatrix.txt";
}

if ( defined $species and $species ne '' and $goindexfile == undef and $gonamesfile == undef ) {
	$goindexfile =
	  "$pagedir/PAGE_DATA/ANNOTATIONS/$species/$species\_go_index.txt";
	$gonamesfile =
	  "$pagedir/PAGE_DATA/ANNOTATIONS/$species/$species\_go_names.txt";
}

my %TU;
if ($tu_enabled) {
	open I, "< $tu_file" or die;
	while (<I>) {
		s/\s+$//;
		my ( $tu, @a ) = split( /\t/, $_ );
		$TU{$tu} = \@a;
	}
}

#
# Loading expfile
#
print "Loading expfile ...";
my $ta = Table->new;
$ta->loadFile($expfile);
my $a_expfile = $ta->getArray();
shift @$a_expfile;
print "Done\n";

#
# Loading motifs
#
$ta->loadFile("$firefile");
my @motif = @{ $ta->getColumn(0) };

my $firebinfile = substr( $firefile, 0, rindex( $firefile, "DNA" ) );
my $filename = Sets::filename($firefile);
print $firebinfile ; <STDIN> ;
$filename =~ s/\.summary/\.profiles/;
$ta->loadFile("$firebinfile/DNA/$filename");
my $a_dna_pro = $ta->getArray();
my %motif_gene;
my %gene_motif;
print "Loading fire profile for DNA...";

foreach my $r (@$a_dna_pro) {
	push( @{ $motif_gene{ $r->[0] } }, $r->[1] );
	push( @{ $gene_motif{ $r->[1] } }, $r->[0] );

	#push(@motif, $r->[0]) if ! (grep {$_ eq $r->[0]} (@motif)) ;
}
print "Done\n";
$ta->loadFile("$firebinfile/RNA/$filename");
my $a_rna_pro = $ta->getArray();
my %isRNA;
if ( $dorna == 1 ) {
	print "Loading fire profile for RNA...";
	foreach my $r (@$a_rna_pro) {
		push( @{ $motif_gene{ $r->[0] } }, $r->[1] );
		push( @{ $gene_motif{ $r->[1] } }, $r->[0] );
		if ( !defined( $isRNA{ $r->[0] } ) ) {
			$isRNA{ $r->[0] } = 1;
		}

		#push(@motif, $r->[0]) if ! (grep {$_ eq $r->[0]} (@motif)) ;
	}
	print "Done\n";
}

print "Loading motif names for DNA...";
my %motif_names;
$filename =~ s/profiles$/motifnames/;
$ta->loadFile("$firebinfile/DNA/$filename");
my $a_motif_name = $ta->getArray();
foreach my $r (@$a_motif_name) {
	my $motif = $r->[0];
	my $name  = $r->[1];
	$name =~ s/^J\_//;
	$name =~ s/^M\S(\d+)\_//;
	$name =~ s/\.txt$//;
	$motif_names{$motif} = $name;
}
print "Done\n";
if ( $dorna == 1 ) {
	print "Loading motif names for RNA...";
	$ta->loadFile("$firebinfile/RNA/$filename");
	$a_motif_name = $ta->getArray();
	foreach my $r (@$a_motif_name) {
		my $motif = $r->[0];
		my $name  = $r->[1];
		$name =~ s/^J\_//;
		$name =~ s/^M\S(\d+)\_//;
		$name =~ s/\.txt$//;
		$motif_names{$motif} = $name;
	}
	print "Done\n";
}

my $GOINDEX = FileHandle->new("< $goindexfile");
my %go_gene;
my %gene_go;
while (<$GOINDEX>) {
	s/\s+$//;
	my ( $gene, @a ) = split( /\t/, $_ );

	foreach my $i ( 0 .. $#a ) {
		push( @{ $go_gene{ $a[$i] } }, $gene );
		push( @{ $gene_go{$gene} }, $a[$i] );
	}
}


print "Running...\n";
my %h_pv;
my $count = 0;
open M, "> $pagedir/TEMP/temp";

foreach my $motif (@motif) {
	print $motif , " (", $count + 1, "/", $#motif + 1, ")";
	my $a_mo_pro;
	foreach my $r (@$a_expfile) {
		my $g = $r->[0];
		if ( grep { $motif eq $_ } ( @{ $gene_motif{$g} } ) ) {
			if ($tu_enabled) {
				next if ( !defined $TU{$g} );
				my @genes = @{ $TU{$g} };
				foreach my $b (@genes) {
					push( @$a_mo_pro, [ $b, 1 ] );
				}
			}
			else {
				push( @$a_mo_pro, [ $g, 1 ] );
			}
		}
		else {
			if ($tu_enabled) {
				next if ( !defined $TU{$g} );
				my @genes = @{ $TU{$g} };
				foreach my $b (@genes) {
					push( @$a_mo_pro, [ $b, 0 ] );
				}
			}
			else {
				push( @$a_mo_pro, [ $g, 0 ] );
			}
		}
	}
	&saveTable( $a_mo_pro, "$pagedir/TEMP/motif_profile.in" );

	my $cmd =
"perl $pagedir/page.pl --expfile=$pagedir/TEMP/motif_profile.in --species=human --exptype=discrete --max_p=0.005";
	system($cmd) ;
	open MI, "< $pagedir/TEMP/motif_profile.in_PAGE/pvmatrix.txt" or die;
	while (<MI>) {
		chomp;
		my ( $go, $pvalue ) = split( /\t/, $_ );
		$h_pv{$go}{$motif} = $pvalue;
		print M "$motif\t$go\t$pvalue\n" ;
	}
	print "\n";
	$count++;
}
close (M) ;
sub transposeTable {
	my $ref_M  = shift @_;
	my $a_ta_M = &copyTable($ref_M);
	my $a_ta_H = shift(@$a_ta_M);
	my $dummy  = shift(@$a_ta_H);
	my $a_ta_R;
	foreach my $r (@$a_ta_M) {
		push( @$a_ta_R, $r->[0] );
		shift @$r;
	}

	my $new_M;
	for ( my $i = 0 ; $i < @$a_ta_H ; $i++ ) {
		for ( my $j = 0 ; $j < @$a_ta_R ; $j++ ) {
			$new_M->[$i]->[$j] = $a_ta_M->[$j]->[$i];
		}
	}
	for ( my $i = 0 ; $i < @$new_M ; $i++ ) {
		my $r = $new_M->[$i];
		unshift( @$r, $a_ta_H->[$i] );
	}
	unshift( @$a_ta_R, $dummy );
	unshift( @$new_M,  $a_ta_R );

	return $new_M;
}

sub copyTable {
	my $a_ta_M = shift @_;
	my $new_M;
	for ( my $i = 0 ; $i < @$a_ta_M ; $i++ ) {
		my $r = $a_ta_M->[$i];
		for ( my $j = 0 ; $j < @$r ; $j++ ) {
			$new_M->[$i]->[$j] = $r->[$j];
		}
	}
	return $new_M;
}

sub saveTable {
	my ( $a_ta_M, $fn, $dummy ) = @_;
	open O, "> $fn" or die "couldn't open $fn...";

	foreach my $r (@$a_ta_M) {
		for ( my $i = 1 ; $i < @$r ; $i++ ) {
			if ( !defined $r->[$i] ) {
				$r->[$i] = $dummy;
			}
		}
		print O join( "\t", @$r ) . "\n";
	}
	close O;
}
