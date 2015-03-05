my $pagedir;

BEGIN {
	if ( ( !$ENV{PAGEDIR} ) || ( $ENV{PAGEDIR} eq '' ) ) {
		$pagedir = "./";
		print
"The PAGEDIR environment variable is not set. It is set to default.\n";
	}
	else {
		$pagedir = $ENV{PAGEDIR};
	}
}

use lib "$pagedir/SCRIPTS";
use lib "$pagedir/SCRIPTS/PostScript-Simple-0.07/lib";

my $programdir = $pagedir . "/PROGRAMS";
my $scriptdir  = $pagedir . "/SCRIPTS";

use Table;
use Sets;
use Getopt::Long;
use PostScript::Simple;
use AggloClust;
use strict;
use Data::Dumper;

my $pvaluematrixfile    = undef;
my $expfile             = undef;
my $datafile            = undef;
my $colmap              = "$scriptdir/HEATMAPS/cmap_dens.txt";
my $cluster             = 5;
my $sortrowsbyphase     = 1;
my $max_p               = undef;
my $minmax_lp           = 3;
my $quantized           = 1;
my $min                 = undef;
my $max                 = undef;
my $xsize               = undef;
my $xscale              = 5;
my $yscale              = 75;
my $scalefont           = 9;
my $h                   = undef;
my $w                   = undef;
my $draw_sample_heatmap = "false";
my $order               = 1;

GetOptions(
	'pvaluematrixfile=s'    => \$pvaluematrixfile,
	'cluster=s'             => \$cluster,
	'expfile=s'             => \$expfile,
	'datafile=s'            => \$datafile,
	'quantized=s'           => \$quantized,
	'colmap=s'              => \$colmap,
	'minmax_lp=s'           => \$minmax_lp,
	'min=s'                 => \$min,
	'max=s'                 => \$max,
	'xsize=s'               => \$xsize,
	'xscale=s'              => \$xscale,
	'yscale=s'              => \$yscale,
	'scalefont=s'           => \$scalefont,
	'h=s'                   => \$h,
	'w=s'                   => \$w,
	'max_p=s'               => \$max_p,
	'draw_sample_heatmap=s' => \$draw_sample_heatmap,
	'order=s'               => \$order
);

my $file    = Sets::filename($expfile);
my $dir     = $expfile . "_PAGE/";
my $outfile = "$expfile" . "\_PAGE/$file" . ".summary.html";

my $ta = Table->new;
$ta->loadFile($pvaluematrixfile);

# get an 2D array
my $a_ref_M = $ta->getArray();
if(scalar @$a_ref_M == 1) {
die "nothing to output for html page";
}

# header
my $a_ref_H = shift @$a_ref_M;
shift @$a_ref_H;
if ( !defined($max_p) ) {
	$max_p = 0.05 / @$a_ref_H;
}
print "Done.\n";
if ($cluster>@$a_ref_M){
    $cluster = scalar(@$a_ref_M) ;
}
if ( defined($cluster) && ( $cluster > 2 ) && ( @$a_ref_M > 2 ) ) {

	print "Cluster rows .. ";

	my $ac = AggloClust->new;

	my @dist = ();
	my $n    = @$a_ref_M;
	for ( my $i = 0 ; $i < $n - 1 ; $i++ ) {
		$dist[$i][$i] = 0;
		for ( my $j = $i + 1 ; $j < $n ; $j++ ) {
			my @a1 = @{ $a_ref_M->[$i] };
			shift @a1;
			my @a2 = @{ $a_ref_M->[$j] };
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
			push @NEWMAT, $a_ref_M->[$i];
		}
	}
	$a_ref_M = \@NEWMAT;

	print "Done.";
}

if ( $order == 1 ) {
	my @NEWMAT = ();
	my @enr;
	for ( my $i = 0 ; $i < @$a_ref_M ; $i++ ) {
		my $r   = $a_ref_M->[$i];
		my $min = 0;
		my $pos = 1;
		for ( my $j = 1 ; $j < @$r ; $j++ ) {
			my $lpo = -1*$r->[$j] ;
	    	if ($lpo<$min and $r->[$j]>0){
				$min = $lpo;
				$pos = $j;
			}
		}
		$enr[$i]->{ind} = $i;
		$enr[$i]->{pos} = $pos;
	}
	my @enr = sort { $b->{pos} <=> $a->{pos} } (@enr);
	foreach my $e (@enr) {
		my $i = $e->{ind};
		push @NEWMAT, $a_ref_M->[$i];
	}
	$a_ref_M = \@NEWMAT;
}

my $A_REF_COLMAP = undef;
if ( defined($colmap) ) {

	$ta->setDelim(" ");
	$ta->loadFile($colmap);
	$A_REF_COLMAP = $ta->getArray();
	$ta->setDelim('\t');
}

# set min max p-values
$max = $minmax_lp  if !( defined $max );
$min = -$minmax_lp if !( defined $min );

my $xsize = 600;
my $xbase = 120;
my $ybase = 130;
my $w     = int( 0.5 + $xsize / @$a_ref_H );
if ( $w < 10 ) {
	$w     = 10;
	$xsize = $w * @$a_ref_H + 500;
}
else {
	$xsize = $w * @$a_ref_H + 500;
}

my $h = 30;

open O, "> $outfile";
print O (
	"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" 
   \"http://www.w3.org/TR/html4/strict.dtd\">\n"
);
print O ("<head>\n");
print O ("<title>iPAGE results</title>\n");
print O ("<style type=\"text/css\">");
print O ("
* {
padding: 0;
margin: 0;
}
body {
font-family: Helvetica, Arial, sans-serif;
}
h1, p, p.head {
width:960px; 
padding: 10px;
}
td.grid {
width: $w\px;
height: $h\px;
}
td.gridheader {
height: auto;
width: $w\px;
border: 1px solid #FFF;
background-color: #000;
text-align: center;
}
td.gridheaderQ {
height: auto;
width: $w\px;
border: 1px solid #FFF;
text-align: center;
}
td.gridlabel {
white-space: nowrap;
height: $h\px;
padding-left: 5px;
}
table.legend {
margin-left: 10px;
top: 0px;
left: 0px;
position: relative;
}
table.legend td {
text-align: center;
}
table.legend tr {
font-size:0px;
height: 2px;
}
table.main {
position: absolute;
top: $ybase\px;
left: $xbase\px;
padding: 10px;
}
");
print O ("</style>\n");
print O ("</head>\n");
print O ("<body>\n\n");


print O ("<h1><abbr title=\"Information-theoretic Pathway Analysis of Gene Expression\">iPAGE</abbr> results</h1>");

print O ("<p class=\"head\">These results are also available as <a href=\"./$file.summary.pdf\"><abbr title=\"Portable Document Format\">PDF</abbr></a> and <a href=\"./$file.summary.eps\"><abbr title=\"Encapsulated PostScript\">EPS</abbr></a> documents.</p>\n");

print O ("<p>Depending on your display resolution, scrolling or zooming may be necessary.</p>");

my $dname = 'imgsrc';
my $d     = "$dir/$dname";
if ( !-e $d ) {
	mkdir $d;
}

####### Scale ########
my $scalex = $xbase / 3;

my $hscale = $h;
my $W      = $hscale;

print O ("<table class=\"legend\" width=\"$hscale\px\" cellpadding=\"0\" cellspacing=\"0\">\n");

my $p = new PostScript::Simple(
	xsize  => $W / 2,
	ysize  => $h * 2,
	colour => 1,
	eps    => 1,
	units  => "pt"
);
$p->setcolour("black");
$p->setfont( "Courier", 10 );
$p->text( { rotate => 90 }, 5, 5, "over-rep" );
$p->output("$d/over.eps");
system("convert -density 100 $d/over.eps $d/over.png");
print O
"<tr><td style=\"border:1px #FFFFFF solid;\">";
print O
"<img src=\"./$dname/over.png\" alt=\"over-representation\" /></td></tr>\n";

print O "<tr><td style=\"font-size: 12px; height: 14px\">$max</td></tr>\n";
my $t   = $max;
my $res = 70;
for ( my $i = 0 ; $i <= $res ; $i++ ) {
	my @col = ();
	if ( !defined($colmap) ) {
		if ( $i > $res / 2 ) {
			@col =
			  Sets::interp_general( $t, [ 0, 0, 0 ], [ 255, 0, 0 ], $min, 0 );
		}
		else {
			@col = Sets::interp_general(
				$t,
				[ 255, 0,   0 ],
				[ 255, 255, 0 ],
				0, $max
			);
		}
	}
	else {
		@col =
		  Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max );
	}

	if ( $i == $res / 2 ) {
		print O "<tr>\n";
		print O
		  "<td style=\"background-color:#FFFFFF; height: 14px; font-size: 12px;\">0</td>\n";
		print O "</tr>\n";
	}
	else {
		my $color = &RGBDecToHex( $col[0], $col[1], $col[2] );
		print O "<tr>\n";
		print O "<td style=\"background-color:$color\">&nbsp;</td>\n";
		print O "</tr>\n";
		$t -= ( $max - $min ) / $res;
	}
}

print O "<tr style=\"font-size: 12px; height: 14px\"><td>$min</td></tr>\n";
my $p = new PostScript::Simple(
	xsize  => $W / 2,
	ysize  => $h * 2,
	colour => 1,
	eps    => 1,
	units  => "pt"
);
$p->setcolour("black");
$p->setfont( "Courier", 10 );
$p->text( { rotate => 90 }, 5, 5, "under-rep" );
$p->output("$d/under.eps");
system("convert -density 100 $d/under.eps $d/under.png");
print O
"<tr><td style=\"border:1px #FFFFFF solid;\">";
print O
"<img src=\"./$dname/under.png\" alt=\"under-representation\" /></td></tr>\n";

print O ("\n</table>\n");

my $min_i1 = 1000000;
my $max_i2 = -1000000;
my @bins   = ();
if ( $quantized == 0 ) {
	foreach my $c (@$a_ref_H) {
		my ( $i1, $i2 ) = $c =~ /\[(.+?)\ (.+?)\]/;
		$min_i1 = $i1 if ( $i1 < $min_i1 );
		$max_i2 = $i2 if ( $i2 > $max_i2 );
		my @a_tmp = ( $i1, $i2 );
		push @bins, \@a_tmp;
	}
	my $l = $xbase - 20;
	my $t = $ybase + 12;
	$l = sprintf( "%2.2f", $l );
	print O
"<div style=\"font-size:9px; top:$t\px; left:$l\px; position:absolute;\">$max_i2</div>";
	my $t = $ybase + $h * 1.5 + 8;
	$l = sprintf( "%2.2f", $l );
	print O
"<div style=\"font-size:9px; top:$t\px; left:$l\px; position:absolute;\">$min_i1</div>";

}

print O (
"<table class=\"main\" cellpadding=\"0\" cellspacing=\"0\">\n"
);
if ( $quantized == 1 ) {
	print O "<tr>\n";
	for ( my $j = 0 ; $j < @$a_ref_H ; $j++ ) {
	  my $p = new PostScript::Simple(
					 xsize  => $w ,
					 ysize  => $h * 2,
					 colour => 1,
					 eps    => 1,
					 units  => "pt"
					);
	  $p->setcolour("black");
	  my $fsize = 8 ;
	  $fsize = $w/2 if ($w/2>8) ;
	  $fsize = 16 if ($fsize>16) ;

	  $p->setfont( "Courier", $w/2 );
	  $p->text( { rotate => 90 }, 5, 5, $a_ref_H->[$j]);
	  my $fn = "C".$a_ref_H->[$j] ;
	  $p->output("$d/$fn.eps");
	  system("convert -density 100 $d/$fn.eps $d/$fn.png");

		print O
"<td class=\"gridheaderQ\" >";
		print O "<img src=\"./$dname/$fn.png\" width=\"$w\px\" alt=\"$fn\" />";
		print O "</td>\n";
	}
	print O
"<td style=\"border:1px #FFFFFF solid;\" >";

	print O "&nbsp;</td>\n";
	print O "</tr>\n";
} else {
	my $th = $max_i2 - $min_i1;
	$th = 1e-3 if ($th < 1e-3) ;

	my $H = $h * 1.5;
	print O "<tr style=\"height:$H\px\">\n";
	$H -= 3;
	my $j = 0;
	foreach my $c (@$a_ref_H) {
		my $h1 = $H * ( $bins[$j]->[0] - $min_i1 ) / $th;
		my $h2 = $H * ( $bins[$j]->[1] - $min_i1 ) / $th;

		my $s =
		    "lower bound = "
		  . $bins[$j]->[0]
		  . " and upper bound = "
		  . $bins[$j]->[1];
		print O
"<td class=\"gridheader\" onclick=\"alert('$s')\" title=\"$s\">";
		my $dw = $w;
		my $dh = $h2 - $h1;
		$dh = 1 if ( $dh < 1 );
		my $y = ( $h1 + $h2 ) / 2 - ( $H - 2 ) / 2;
		print "$y\t$H\t$h1\n";
		print O
"<div style=\"width:100\%; height:$dh\px; background-color:#F00 ; left:0\px ; bottom:$y\px; position:relative\">";
		print O "</div>\n";

		print O "</td>\n";
		$j++;
	}

	print O
	  "<td style=\"border:1px #FFFFFF solid;\">";
	print O "&nbsp;</td>\n";
	print O "</tr>\n";
}

for ( my $i = 0 ; $i < @$a_ref_M ; $i++ ) {
	print O "<tr>\n";

	# get row
	my $r = $a_ref_M->[$i];

	# get GO description
	my $go = shift @$r;

	# go through entries in current row
	for ( my $j = 0 ; $j < @$r ; $j++ ) {

		# get a single value, pos if over-rep, neg if under-rep
		my $lp = $r->[$j] ;
		my $pv;
		my $s;
		if ( $lp>0 ) {
			$pv = -1*$lp;
			$s  = "over-representation";
		}
		else {
			$pv = $lp;
			$s  = "under-representation";
		}
		$pv = 10**$pv;
		$pv = sprintf( "%s: log(p) = %1.4e", $s, $pv );

		# create appropriate color
		my @col = ();    #$colmap = undef;
		if ( !defined($colmap) ) {
			if ( $lp < 0 ) {
				@col = Sets::interp_general(
					$lp,
					[ 0, 255, 0 ],
					[ 0, 0,   0 ],
					$min, 0
				);
			}
			else {
				@col = Sets::interp_general(
					$lp,
					[ 0,   0, 0 ],
					[ 255, 0, 0 ],
					0, $max
				);
			}
		}
		else {
			@col =
			  Sets::interp_from_matlab_colormap( $lp, $A_REF_COLMAP, $min,
				$max );
		}

		my $color = &RGBDecToHex( $col[0], $col[1], $col[2] );
		print O "<td class=\"grid\" style=\"background-color:$color\" onclick=\"alert('$pv')\" title=\"$pv\">&nbsp;</td>\n";
	}
	print O "<td class=\"gridlabel\">$go</td>\n";
	print O "</tr>\n";
}
print O ("\n</table>\n");
print O ("\n</body>\n");
print O ("</html>\n");

sub RGBDecToHex {
	my ( $red, $green, $blue ) = @_;
	return sprintf( "#%02lx%02lx%02lx", $red, $green, $blue );
}
