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
my $species     = undef;
my $gonamesfile = undef;
my $max_p       = 0.05;
my $submit      = 0;
my $printmode   = 0;
my $tu_enabled  = 0;
my $tu_file     = undef;
my $dorna       = 0;

my $holdjob0 = undef;
my $holdjob1 = undef;

my $annochoice = 0;
my $protein = 0;
my $homologiesfile=undef;
my $genelist=undef;
my $nodups = 0;

GetOptions(
	'expfile=s'     => \$expfile,
	'firefile=s'    => \$firefile,
	'pagefile=s'    => \$pagefile,
	'goindexfile=s' => \$goindexfile,
	'species=s'     => \$species,
	'max_p=s'       => \$max_p,
	'submit=s'      => \$submit,
	'tu_enabled=s'  => \$tu_enabled,
	'tu_file=s'     => \$tu_file,
	'printmode=s'   => \$printmode,
	'dorna=s'       => \$dorna,
	'holdjob0=s'    => \$holdjob0,
	'holdjob1=s'    => \$holdjob1,
	'annochoice=i' => \$annochoice,
	'protein=i' => \$protein,
	'homologiesfile=s' => \$homologiesfile,
	'genelist=s' => \$genelist,
	'nodups=i' => \$nodups,
);

# protein add
if (defined $species && $protein == 2 && ($species eq 'human' || $species eq 'mouse' || $species eq 'ecoli')) {
	$species .= 'p';
}
# protein end

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

		#my $firefile = $file . "_FIREPRO/$f.final.motifs";
		my $pagefile = $file . "_PAGE/pvmatrix.txt";

#		my $cmd = "perl $pagedir/prmg.pl --expfile=$file --firefile=$firefile --pagefile=$pagefile --goindexfile=$goindexfile --max_p=$max_p --dorna=$dorna --species=$species";
		my $cmd = "perl $pagedir/prmg_protein.pl --expfile=$file --pagefile=$pagefile --goindexfile=$goindexfile --max_p=$max_p --dorna=$dorna --species=$species  --homologiesfile=$homologiesfile --genelist=$genelist";

		$pbs->addCmd($cmd);

		my $generatehtmlcmd =
		  "perl $pagedir/SCRIPTS/motif_cat_draw_html_protein.pl --expfile=$file";
		$pbs->addCmd($generatehtmlcmd);

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

## done with job submission; below is actual processing

# add remove duplicates
## can only happen if species is known or if a homologies, gene list files are supplied ($homologiesfile, $genelist)
my $removedups = 0;
if ($nodups == 0) {
	if (defined $species || (defined $homologiesfile && defined $genelist)) {
	
		if (defined $species && !defined $homologiesfile && !defined $genelist) {
			$homologiesfile = "$pagedir/PAGE_DATA/ANNOTATIONS/$species/$species\.homologies";
			$genelist = "$pagedir/PAGE_DATA/ANNOTATIONS/seqnames/$species\.txt";
		}
		
		my $nodupsexptype=0;

		my $fileout = "$expfile\.nodups";
		my $removedupcmd = "perl $scriptdir/remove_homologous_sequences_withseed.pl -expfile $expfile -quantized $nodupsexptype -genelist $genelist -dupfile $homologiesfile -outfile $fileout";

		print "REMOVING DUPLICATES\n";
		print $removedupcmd,"\n\n";
		system($removedupcmd);
		$removedups = 1;
		system("mv $expfile $expfile\.original");
		system("cp $expfile\.nodups $expfile");
	}
}
# end add remove duplicates


if ( $printmode == 1 ) {
	if ( !defined($pagefile) ) {
		$pagefile = "$expfile\_PAGE/pvmatrix.txt";
	}

	my $dir = Sets::dirname($pagefile);

	my $ta = Table->new;
	$ta->loadFile("$dir/motif_cat.cdt");
	my $a_ta_M = $ta->getArray();

#$my ($a_ta_M, $fn, $vline, $hline, $vcol, $min, $mid, $max, $h, $w, $xplus, $yplus, $header_motif, $row_motif, $s_y, $s_h, $s_w, $sep, $res, $scalefont,$uppertext, $lowertext,$lcol, $hcol, $mcol) = @_ ;
	&printHeatmapBars(
		$a_ta_M, "$dir/motif_cat.pdf",
		"true",  "true",
		[ 0, 0, 0 ], -3,
		0,      3,
		35,     30,
		560,    50,
		"true", "false",
		0,      2,
		20,     8,
		30,     12,
		"Pos",  "Neg",
		[ 0, 0, 255 ], [ 255, 0, 0 ],
		[ 255, 255, 255 ]
	);

	die;
}

# build FIRE-PRO file paths; replaces firefile due to name inconsistencies
my $FP_summary = undef;
my $FP_profile = undef;
my $FP_names = undef;

if ( !defined($FP_summary)) {
	my $fn = Sets::filename($expfile);
	$FP_summary = $expfile.'_FIREPRO/'.$fn.'.final.motifs';
}

if ( !defined($FP_profile)) {
	my $fn = Sets::filename($expfile);
	$FP_profile = $expfile.'_FIREPRO/'.$fn.'-motif_profiles.txt';
}

if ( !defined($FP_names)) {
	my $fn = Sets::filename($expfile);
	$FP_names = $expfile.'_FIREPRO/'.$fn.'.names';
}

if ( !defined($pagefile) ) {
	$pagefile = "$expfile\_PAGE/pvmatrix.txt";
}

if ( defined $species && $species ne '' && !defined $goindexfile ) {
	if($annochoice == 2) {
		$goindexfile = "$pagedir/PAGE_DATA/ANNOTATIONS/$species/$species\_index_k.txt" ;
		if ((!-e $goindexfile)) {
			die ("FAILURE: KEGG annotations not available.");
		}
	} elsif ($annochoice == 1) {
		$goindexfile = "$pagedir/PAGE_DATA/ANNOTATIONS/$species/$species\_index_gk.txt" ;
		if ((!-e $goindexfile)) {
			$goindexfile = "$pagedir/PAGE_DATA/ANNOTATIONS/$species/$species\_index.txt" ;
		}
	} else {
		$goindexfile = "$pagedir/PAGE_DATA/ANNOTATIONS/$species/$species\_index.txt" ;
	}
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
# Loading pvmatrix file
#
print "Loading pagefile ...";
$ta->loadFile($pagefile);
my $a_pv_page = $ta->getArray();
shift @$a_pv_page;
my $a_catlist;
my $c = 0;
foreach my $r (@$a_pv_page) {
	$r->[0] =~ s/,//;
	my (@f) = split( /\s/, $r->[0] );
	my $go  = shift(@f);
	my $cat = join( " ",   @f );
	$a_catlist->[$c]->[0] = $go;
	$a_catlist->[$c]->[1] = $cat;
	if ( $cat eq "" ) {
		$a_catlist->[$c]->[1] = $go;
	}
	$a_catlist->[$c]->[2] = "P";
	$c++;
}
my $tmpfile0001 = Sets::getTempFile("$pagedir/TEMP/go_names.in");
&saveTable( $a_catlist, $tmpfile0001 );
print "Done\n";

#
# Loading motifs
#
print "Loading fire summary file...\n";
print $FP_summary,"\n";
$ta->loadFile($FP_summary);
my @motif = @{ $ta->getColumn(0) };
print "Done\n";

$ta->loadFile($FP_profile);
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

print "Loading motif names for DNA...";
my %motif_names;
my $a_motif_name;

if (-s $FP_names) {
	$ta->loadFile($FP_names);
	$a_motif_name = $ta->getArray();
	foreach my $r (@$a_motif_name) {
		my $motif = $r->[1];
		my $name  = $r->[3];
		$name =~ s/^J\_//;
		$name =~ s/^M\S(\d+)\_//;
		$name =~ s/\.txt$//;
		$motif_names{$motif} = $name;
	}
}
print "Done\n";

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
print scalar(@motif);
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
	my $tmpfile0002 = Sets::getTempFile("$pagedir/TEMP/motif_profile.in");
	&saveTable( $a_mo_pro, $tmpfile0002 );

	my $tmpfile0003 = Sets::getTempFile("$pagedir/TEMP/temp");

	my $cmd =
"$pagedir/PROGRAMS/mi_go_motif_calculator -expfile $tmpfile0002 -gonamesfile $tmpfile0001 -goindexfile $goindexfile -P 1 -quantized 1 -independence 0 -max_p $max_p -pvaluematrixfile $tmpfile0003";
	system($cmd) ;
	open MI, "< $tmpfile0003" or die;
	while (<MI>) {
		chomp;
		my ( $go, $pvalue ) = split( /\t/, $_ );
		$h_pv{$go}{$motif} = $pvalue;
	}
	print "\n";
	$count++;
	system("rm $tmpfile0002");
	system("rm $tmpfile0003");
}

my $tmpfile0003 = Sets::getTempFile("$pagedir/TEMP/temp");

open M, "> $tmpfile0003";
foreach my $motif (@motif) {
	my $m = $motif;
#	if ( $isRNA{$motif} == 1 ) {
#		$m =~ s/T/U/gi;
#		print M "\t$m/3' ";
#	}
#	else {
		print M "\t$m/ ";
#	}
	if ( defined $motif_names{$motif} and $motif_names{$motif} ne "-" ) {
		print M $motif_names{$motif};
	}
}
print M "\n";
foreach my $r (@$a_catlist) {
	print M $r->[1], " ", $r->[0];
	print $r->[1], " ", $r->[0], "\n";

	foreach my $motif (@motif) {
		print M "\t", $h_pv{ $r->[1] }{$motif}
		  if ( defined $h_pv{ $r->[1] }{$motif} );
		print M "\t0" if !( defined $h_pv{ $r->[1] }{$motif} );
	}
	print M "\n";
}
close M;

$ta->loadFile($tmpfile0003);
my $a_ta_M = $ta->getArray();
my $a_n_M;

my $cluster;

push( @$a_n_M, $a_ta_M->[0] );
for ( my $i = 1 ; $i < @$a_ta_M ; $i++ ) {
	my $r   = $a_ta_M->[$i];
	my $sum = 0;
	for ( my $j = 1 ; $j < @$r ; $j++ ) {
		$sum += $r->[$j];
	}
	if ( $sum != 0 ) {
		push( @$a_n_M, $a_ta_M->[$i] );
	}
}
$a_ta_M = $a_n_M;
$a_n_M  = undef;
if (scalar(@$a_ta_M) == 1) {
	print "Table is empty. No results to be returned.","\n";
	exit;
}
$a_ta_M = &transposeTable($a_ta_M);
push( @$a_n_M, $a_ta_M->[0] );
for ( my $i = 1 ; $i < @$a_ta_M ; $i++ ) {
	my $r   = $a_ta_M->[$i];
	my $sum = 0;
	for ( my $j = 1 ; $j < @$r ; $j++ ) {
		$sum += $r->[$j];
	}
	if ( $sum != 0 ) {
		push( @$a_n_M, $a_ta_M->[$i] );
	}
}
$a_ta_M = $a_n_M;
$a_ta_M = &transposeTable($a_ta_M);

if ( @$a_ta_M > 3 ) {
	$cluster = int( sqrt( scalar(@$a_ta_M) - 1 ) );
	print "Clustering rows into $cluster clusters...";

	my $ac = AggloClust->new;

	my $H    = shift @$a_ta_M;
	my @dist = ();
	my $n    = @$a_ta_M;
	for ( my $i = 0 ; $i < $n - 1 ; $i++ ) {
		$dist[$i][$i] = 0;
		for ( my $j = $i + 1 ; $j < $n ; $j++ ) {
			my @a1 = @{ $a_ta_M->[$i] };
			shift @a1;
			my @a2 = @{ $a_ta_M->[$j] };
			shift @a2;
			$dist[$i][$j] = 1 - Sets::pearson( \@a1, \@a2 );
			$dist[$j][$i] = $dist[$i][$j];
		}
	}

	$ac->setDistanceMatrix( \@dist );
	$ac->setMaxNbClusters($cluster);
	my $a_ref_c = $ac->agglomerate_using_avg_linkage();

	my @NEWMAT = ();
	foreach my $c (@$a_ref_c) {
		print join( " ", @$c );
		print "\n";
		foreach my $i (@$c) {
			push @NEWMAT, $a_ta_M->[$i];
		}
	}
	$a_ta_M = \@NEWMAT;
	unshift( @$a_ta_M, $H );
	print "Done.\n";
}
$a_ta_M = &transposeTable($a_ta_M);
if ( @$a_ta_M > 3 ) {
	$cluster = int( sqrt( scalar(@$a_ta_M) - 1 ) );
	print "Clustering columns into $cluster clusters...";

	my $ac = AggloClust->new;

	my $H    = shift @$a_ta_M;
	my @dist = ();
	my $n    = @$a_ta_M;
	for ( my $i = 0 ; $i < $n - 1 ; $i++ ) {
		$dist[$i][$i] = 0;
		for ( my $j = $i + 1 ; $j < $n ; $j++ ) {
			my @a1 = @{ $a_ta_M->[$i] };
			shift @a1;
			my @a2 = @{ $a_ta_M->[$j] };
			shift @a2;
			$dist[$i][$j] = 1 - Sets::pearson( \@a1, \@a2 );
			$dist[$j][$i] = $dist[$i][$j];
		}
	}

	$ac->setDistanceMatrix( \@dist );
	$ac->setMaxNbClusters($cluster);
	my $a_ref_c = $ac->agglomerate_using_avg_linkage();

	my @NEWMAT = ();
	foreach my $c (@$a_ref_c) {
		print join( " ", @$c );
		print "\n";
		foreach my $i (@$c) {
			push @NEWMAT, $a_ta_M->[$i];
		}
	}
	$a_ta_M = \@NEWMAT;
	unshift( @$a_ta_M, $H );
	print "Done.\n";
}
$a_ta_M = &transposeTable($a_ta_M);

my $dir = Sets::dirname($pagefile);
&saveTable( $a_ta_M, "$dir/motif_cat.cdt" );

my $ta = Table->new;
$ta->loadFile("$dir/motif_cat.cdt");
my $a_ta_M = $ta->getArray();

#$my ($a_ta_M, $fn, $vline, $hline, $vcol, $min, $mid, $max, $h, $w, $xplus, $yplus, $header_motif, $row_motif, $s_y, $s_h, $s_w, $sep, $res, $scalefont,$uppertext, $lowertext,$lcol, $hcol, $mcol) = @_ ;
&printHeatmapBars(
	$a_ta_M, "$dir/motif_cat.pdf",
	"true",  "true",
	[ 0, 0, 0 ], -3,
	0,      3,
	35,     30,
	560,    50,
	"true", "false",
	0,      2,
	20,     8,
	30,     12,
	"Pos",  "Neg",
	[ 0, 0, 255 ], [ 255, 0, 0 ],
	[ 255, 255, 255 ]
);

system("rm $tmpfile0001");
system("rm $tmpfile0003");

sub printHeatmapBars {
	my (
		$a_ta_M,    $fn,        $vline,        $hline,     $vcol,
		$min,       $mid,       $max,          $h,         $w,
		$xplus,     $yplus,     $header_motif, $row_motif, $s_y,
		$s_h,       $s_w,       $sep,          $res,       $scalefont,
		$uppertext, $lowertext, $lcol,         $hcol,      $mcol
	  )
	  = @_;

	my $a_ta_H = shift(@$a_ta_M);
	my $dummy  = shift(@$a_ta_H);

	my $xbase = 100;
	my $ybase = 200;
	my $ysize = $ybase + $h * scalar(@$a_ta_M) + $yplus;
	my $xsize = $xbase + $w * scalar(@$a_ta_H) + $xplus;

	my $p = new PostScript::Simple(
		xsize  => $xsize,
		ysize  => $ysize,
		colour => 1,
		eps    => 1,
		units  => "pt"
	);
	$p->setlinewidth(0.5);

	if ( $header_motif eq "true" or $row_motif eq "true" ) {
		system("mkdir ./Temp") if !( -d "./Temp" );
	}

	for ( my $j = 0 ; $j < @$a_ta_H ; $j++ ) {
		if ( $header_motif eq "true" ) {
			my $motif = $a_ta_H->[$j];
			$motif =~ /(\S+)\/(.+)$/;
			$motif = $1 if defined $1;
			my $name = $2;
			$motif =~ s/\/$//;
			print $motif, "\t", $name, "\n";
			my $mo = Sets::myre2wm($motif);
			open OUT, "> ./Temp/$motif.txt" or die "cannot open\n";
			print OUT $mo;
			close OUT;
			system(
"$pagedir/SCRIPTS/weblogo/seqlogo -f ./Temp/$motif.txt -F EPS -k 0 -a -c -M -n -Y -w 5 -h 3 > ./Temp/$motif.eps"
			);

			my $e  = new PostScript::Simple::EPS( file => "./Temp/$motif.eps" );
			my $eh = $e->height;
			my $ew = $e->width;

			# height must be $h, so scale down to $h = k * LO
			$e->scale( $w / $eh );
			$e->rotate(90);
			my $ew_new = int( 0.5 + $ew * $h / $eh );
			$p->_add_eps(
				$e,
				$xbase + ( $j + 1.1 ) * $w,
				$ysize - ( $ybase - 3 )
			);
			if ( defined $name or $name ne "" ) {
				$p->setcolour("black");
				$p->setfont( "TimesBold", 10 );
				$p->text(
					{ rotate => 90, align => "left" },
					$xbase + $j * $w + 3 * $w / 4,
					$ysize - ( $ybase - 5 - $w / $eh * $ew ), $name
				);
			}
			else {
				$p->setcolour("black");
				$p->setfont( "TimesBold", 10 );
				$p->text(
					{ rotate => 90, align => "left" },
					$xbase + $j * $w + 3 * $w / 4,
					$ysize - ( $ybase - 5 - $w / $eh * $ew ), "-"
				);
			}
		}
		else {
			$p->setcolour("black");
			$p->setfont( "Garamond", 6 );
			$p->text(
				{ rotate => 90 },
				$xbase + $j * $w + 3 * $w / 4,
				$ysize - ( $ybase - 3 ),
				$a_ta_H->[$j]
			);
		}
	}

	my @col = ();
	for ( my $i = 0 ; $i < @$a_ta_M ; $i++ ) {
		my $r    = $a_ta_M->[$i];
		my $name = shift @$r;
		print $name, "\n";
		for ( my $j = 0 ; $j < @$r ; $j++ ) {
			my $v = $r->[$j];

			#defining the color
			if ( defined $mcol ) {
				if ( $v > $mid ) {
					@col = Sets::interp_general( $v, $mcol, $hcol, $mid, $max );
				}
				else {
					@col = Sets::interp_general( $v, $lcol, $mcol, $min, $mid );
				}
			}
			else {
				@col = Sets::interp_general( $v, $lcol, $hcol, $min, $max );
			}

			$p->setcolour(@col);
			$p->box(
				{ filled => 1 },
				$xbase + $j * $w,
				$ysize - ( $ybase + $i * $h ),
				$xbase + $j * $w + $w,
				$ysize - ( $ybase + ( $i * $h + $h ) )
			);

			$p->setcolour(@$vcol);
			$p->line(
				$xbase + $j * $w,
				$ysize - $ybase,
				$xbase + $j * $w,
				$ysize - ( $ybase + ( $i + 1 ) * $h )
			  )
			  if $vline eq "true";
			$p->line(
				$xbase,
				$ysize - ( $ybase + $i * $h ),
				$xbase + ( $j + 1 ) * $w,
				$ysize - ( $ybase + $i * $h )
			  )
			  if $hline eq "true";
		}

		if ( $row_motif eq "true" ) {
			my $motif = $name;
			$motif =~ /(\S+)\/(\S+)/;
			$motif = $1 if defined $1;
			my $m_name = $2;
			$motif =~ s/\/$//;
			print $motif, "\t", $m_name, "\n";

			my $mo = Sets::myre2wm($motif);
			open OUT, "> ./Temp/$motif.txt" or die "cannot open\n";
			print OUT $mo;
			close OUT;
			system(
"$pagedir/SCRIPTS/weblogo/seqlogo -f ./Temp/$motif.txt -F EPS -k 0 -a -c -M -n -Y -w 5 -h 3 > ./Temp/$motif.eps"
			);

			my $e  = new PostScript::Simple::EPS( file => "./Temp/$motif.eps" );
			my $eh = $e->height;
			my $ew = $e->width;

			# height must be $h, so scale down to $h = k * LO
			$e->scale( $h / $eh );
			my $ew_new = int( 0.5 + $ew * $h / $eh );
			$p->_add_eps(
				$e,
				$xbase + @$r * $w + 10,
				$ysize - ( $ybase + $i * $h + $h )
			);
			if ( defined $name or $name ne "" ) {
				$p->setcolour("black");
				$p->setfont( "TimesBold", 10 );
				$p->text(
					{ align => "left" },
					$xbase + @$r * $w + 30 + $ew * $h / $eh,
					$ysize - ( $ybase + $i * $h + $h / 2 ), $m_name
				);
			}
			else {
				$p->setcolour("black");
				$p->setfont( "TimesBold", 10 );
				$p->text(
					{ align => "left" },
					$xbase + @$r * $w + 30 + $ew * $h / $eh,
					$ysize - ( $ybase + $i * $h + $h / 2 ), "-"
				);
			}
		}
		else {
			$p->setfont( "Garamond", 15 );
			$p->setcolour("black");
			my ( $p1, $p2 ) = split( /\s/, $name );
			if ( $p1 eq $p2 ) {
				$name = $p1;
			}
			$p->text(
				{ align => "left", rotate => 0 },
				$xbase + @$r * $w + 10,
				$ysize - ( $ybase + $i * $h + $h / 2 + 4 ), $name
			);
		}
	}

	$p->setcolour(@$vcol);
	$p->line(
		$xbase + @$a_ta_H * $w,
		$ysize - $ybase,
		$xbase + @$a_ta_H * $w,
		$ysize - ( $ybase + @$a_ta_M * $h )
	);
	$p->line( $xbase, $ysize - $ybase, $xbase + @$a_ta_H * $w,
		$ysize - $ybase );
	$p->line(
		$xbase,
		$ysize - ( $ybase + @$a_ta_M * $h ),
		$xbase + @$a_ta_H * $w,
		$ysize - ( $ybase + @$a_ta_M * $h )
	);
	$p->line( $xbase, $ysize - $ybase,
		$xbase, $ysize - ( $ybase + @$a_ta_M * $h ) );

	&drawScale(
		$xbase / 4, $ybase + $s_y, $min, $mid,       $max,
		$sep,       $res,          $p,   $xsize,     $ysize,
		$scalefont, $s_h,          $s_w, $uppertext, $lowertext,
		$lcol,      $hcol,         $mcol
	);

	my $outpdf = $fn;
	my $outeps = $outpdf;
	$outeps =~ s/pdf$/eps/;
	$p->output("$outeps");
	system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");

	system("rm -r ./Temp");
		
	unshift( @$a_ta_H, $dummy );
	unshift( @$a_ta_M, $a_ta_H );
}

sub transposeTable {
	my $ref_M  = shift @_;
	my $a_ta_M = &copyTable($ref_M);
	print Dumper($a_ta_M),"\n";
	my $a_ta_H = shift(@$a_ta_M);
	my $dummy  = shift(@$a_ta_H);
	my $a_ta_R;
	foreach my $r (@$a_ta_M) {
		push( @$a_ta_R, $r->[0] );
		shift @$r;
	}

	my $new_M;
	#print Dumper($a_ta_M),"\n";
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
	#print Dumper($new_M);
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

sub drawScale {
	my ( $x, $y, $min, $mid, $max, $sep, $res, $p, $xsize, $ysize, $scalefont,
		$h, $w, $uppertext, $lowertext, $lcol, $hcol, $mcol )
	  = @_;

	$p->setcolour("black");
	$p->setfont( "Courier", $scalefont );
	$p->text(
		{ align => "center" },
		$x + $w / 2 + 0,
		$ysize - ( $y - 3 ), $max
	);
	$p->text(
		{ align => "center" },
		$x + $w / 2 + 0,
		$ysize - ( $y - 13 ), $uppertext
	);

	my $t = $max;
	for ( my $i = 0 ; $i < $res / 2 ; $i++ ) {
		my @col = ();

		if ( defined $mcol ) {
			if ( $t > $mid ) {
				@col = Sets::interp_general( $t, $mcol, $hcol, $mid, $max );
			}
			else {
				@col = Sets::interp_general( $t, $lcol, $mcol, $min, $mid );
			}
		}
		else {
			@col = Sets::interp_general( $t, $lcol, $hcol, $min, $max );
			$t -= ( $max - $min ) / $res;
		}

		$p->setcolour(@col);
		$p->box(
			{ filled => 1 },
			$x, $ysize - ( $y + $i * $h ),
			$x + $w, $ysize - ( $y + $i * $h + $h )
		);
		$t -= ( $max - $min ) / $res;
	}

	$p->setcolour("black");
	$p->setfont( "Courier", $scalefont );

	$p->text(
		{ align => "center" },
		$x + $w / 2 + 0,
		$ysize - ( $y + ( $res / 2 ) * $h + $sep * 3 / 4 ), $mid
	);

	if ( defined $mcol ) {
		for ( my $i = $res - 1 ; $i >= $res / 2 ; $i-- ) {
			my @col = ();

			@col = Sets::interp_general( $t, $mcol, $lcol, $min, $mid );

			$p->setcolour(@col);
			$p->box(
				{ filled => 1 },
				$x, $ysize - ( $y + $sep + $i * $h ),
				$x + $w, $ysize - ( $y + $sep + $i * $h + $h )
			);
			$t -= ( $max - $min ) / $res;
		}
		$p->setcolour('black');
		$p->text(
			{ align => "center" },
			$x + $w / 2 + 1,
			$ysize - ( $y + $sep + $res * $h + 20 ), $lowertext
		);
		$p->text(
			{ align => "center" },
			$x + $w / 2 + 0,
			$ysize - ( $y + $sep + $res * $h + 10 ), $min
		);
	}
}
