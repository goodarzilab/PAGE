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

my $expfile   = undef;
my $min       = -3;
my $max       = 3;
my $xsize     = undef;
my $xscale    = 5;
my $yscale    = 75;
my $scalefont = 9;
my $h         = undef;
my $w         = undef;

if ( @ARGV == 0 ) {
	die "Usage: perl motif_cat_draw_html.pl  --expfile=FILE\n";
}

GetOptions(
	'expfile=s'   => \$expfile,
	'min=s'       => \$min,
	'max=s'       => \$max,
	'xsize=s'     => \$xsize,
	'xscale=s'    => \$xscale,
	'yscale=s'    => \$yscale,
	'scalefont=s' => \$scalefont,
	'h=s'         => \$h,
	'w=s'         => \$w,
);

#
# creating the summary file
#
my $fn      = Sets::filename($expfile);
my $dir     = $expfile . "_PAGE";
my $cdtfile = "$dir/motif_cat.cdt";
my $outfile = "$dir/motif_cat.html";

my $ta = Table->new;

print "Reading matrix ... ";
$ta->loadFile($cdtfile);
my $a_ref_M = $ta->getArray();
my $a_ref_H = shift @$a_ref_M;
shift @$a_ref_H;
print "Done.\n";


my $xbase = 50;
my $ybase = 130;
my $w     = 35;
my $xsize = $w * @$a_ref_H + 500;
my $h     = $w;

open O, "> $outfile";
print O (
	"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" 
   \"http://www.w3.org/TR/html4/strict.dtd\">\n"
);
print O ("<head>\n");
print O ("<title>Pathway regulatory interaction map (PRMG output)</title>\n");
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
height: $w\px;
border: 1px solid #000;
}
td.gridheader {
height: auto;
width: $w\px;
border: 1px solid #FFF;
text-align: center;
}
td.gridlabel {
white-space: nowrap;
height: $w\px;
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

print O ("<h1>Pathway-regulatory interaction map</h1>\n");

print O ("<p class=\"head\">These results are also available as <a href=\"./motif_cat.pdf\"><abbr title=\"Portable Document Format\">PDF</abbr></a> and <a href=\"./motif_cat.eps\"><abbr title=\"Encapsulated PostScript\">EPS</abbr></a> documents.</p>\n");

print O ("<p>Depending on your display resolution, scrolling or zooming may be necessary.</p>\n");

my $dname = 'imgsrc';
my $d     = "$dir/$dname";
if ( !-e $d ) {
	mkdir $d;
}

############ Scale ############
my $scalex = $xbase / 3;
my $hscale = $h * 2 / 3;
my $W      = $hscale;
print O (
"<table class=\"legend\" width=\"$hscale\px\" cellpadding=\"0\" cellspacing=\"0\">\n"
);
my $p = new PostScript::Simple(
	xsize  => $W / 2,
	ysize  => $h + 40,
	colour => 1,
	eps    => 1,
	units  => "pt"
);
$p->setcolour("black");
$p->setfont( "courier", 8 );
$p->text( { rotate => 90 }, 7, 5, "over-rep" );
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
	if ( $i > $res / 2 ) {
		@col =
		  Sets::interp_general( $t, [ 0, 0, 255 ], [ 255, 255, 255 ], $min, 0 );
	}
	else {
		@col =
		  Sets::interp_general( $t, [ 255, 255, 255 ], [ 255, 0, 0 ], 0, $max );
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
		print O "<td style=\"background-color: $color\">&nbsp;</td>\n";
		print O "</tr>\n";
		$t -= ( $max - $min ) / $res;
	}
}

print O "<tr style=\"font-size: 12px; height: 14px\"><td>$min</td></tr>\n";
my $p = new PostScript::Simple(
	xsize  => $W / 2,
	ysize  => $h + 20,
	colour => 1,
	eps    => 1,
	units  => "pt"
);
$p->setcolour("black");
$p->setfont( "courier", 8 );
$p->text( { rotate => 90 }, 7, 5, "under-rep" );
$p->output("$d/under.eps");
system("convert -density 100 $d/under.eps $d/under.png");
print O
"<tr><td style=\"border:1px #FFFFFF solid;\">";
print O
"<img src=\"./$dname/under.png\" alt=\"under-representation\" /></td></tr>\n";

print O ("\n</table>\n");

print O (
"<table class=\"main\" cellpadding=\"0\" cellspacing=\"0\">\n"
);

print O ("<tr>");
for ( my $j = 0 ; $j < @$a_ref_H ; $j++ ) {
	my $motif = $a_ref_H->[$j];
	$motif =~ /(\S+)\/(.+)$/;
	$motif = $1 if defined $1;
	my $name = $2;
	$name =~ s/\s+$//;
	my @namesplit = split( /\//, $name );
	@namesplit = split(/\;/, $namesplit[0]);
	@namesplit = split(/\,/, $namesplit[0]);
	$name = $namesplit[0];
	my $p = new PostScript::Simple(
		xsize  => $w / 2,
		ysize  => 100,
		colour => 1,
		eps    => 1,
		units  => "pt"
	);
	$p->setcolour("black");
	$p->setfont( "Garamond", 10 );
	$p->text( { rotate => 90 }, 10, 0, "$name" );

	$name =~ s/\ /_/gi;
	$name =~ s/'//gi;
	my $outname    = $d . '/' . $name . '.eps';
	my $outnamepng = $d . '/' . $name . '.png';
	$p->output($outname);
	system("convert -density 100 $outname $outnamepng");
	print O
"<td class=\"gridheader\">";
	print O
"<img src=\"./$dname/$name.png\" alt=\"$name\" title=\"$name\"/></td>\n";
}
print O ("</tr>\n");

print O ("<tr>");
for ( my $j = 0 ; $j < @$a_ref_H ; $j++ ) {
	my $motif = $a_ref_H->[$j];
	$motif =~ /(\S+)\/(.+)$/;
	$motif = $1 if defined $1;
	my $name = $2;
	$motif =~ s/\/$//;
	print $motif, "\t", $name, "\n";
	my $mo = Sets::myre2wm($motif);
	open OUT, "> $d/$motif.txt" or die "cannot open\n";
	print OUT $mo;
	close OUT;
	system(
"$pagedir/SCRIPTS/weblogo/seqlogo -f $d/$motif.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $d/a$j.eps"
	);

	my $e  = new PostScript::Simple::EPS( file => "$d/a$j.eps" );
	my $eh = $e->height;
	my $ew = $e->width;
	$e->rotate(90);
	my $p = new PostScript::Simple(
		xsize  => $eh,
		ysize  => $ew,
		colour => 1,
		eps    => 1,
		units  => "pt"
	);
	$p->_add_eps( $e, $ew * 2 / 3, 0 );
	my $outname    = $d . '/a' . $j . '.eps';
	my $outnamepng = $d . '/a' . $j . '.png';
	$p->output("$outname");

	system("convert -density 100 $outname $outnamepng");

	print O
"<td class=\"gridheader\">";
	print O
"<img src=\"./$dname/a$j.png\" width=\"$w\px\" alt=\"$motif\" title=\"$motif\" /></td>\n";
}
print O ("</tr>");

for ( my $i = 0 ; $i < @$a_ref_M ; $i++ ) {
	print O
	  "<tr>\n";
	my $r    = $a_ref_M->[$i];
	my $name = shift @$r;
	print $name, "\n";
	for ( my $j = 0 ; $j < @$r ; $j++ ) {
		my $v = $r->[$j];

		#defining the color
		my @col = ();
		if ( $v < 0 ) {
			@col = Sets::interp_general(
				$v,
				[ 0,   0,   255 ],
				[ 255, 255, 255 ],
				$min, 0
			);
		}
		else {
			@col = Sets::interp_general(
				$v,
				[ 255, 255, 255 ],
				[ 255, 0,   0 ],
				0, $max
			);
		}
		my $color = &RGBDecToHex( $col[0], $col[1], $col[2] );
		my $pv;
		if ( $v > 0 ) {
			$pv = 10**( -1 * $v );
		}
		else {
			$pv = 10**$v;
		}
		$pv = sprintf( "log(p) = %1.4e", $pv );

		if ( abs($v) > 0 ) {
			if ( $v > 0 ) {
				print O
"<td class=\"grid\" style=\"background-color:$color\" title=\"over-representation: $pv\" onclick=\"alert('over-representation: $pv')\">&nbsp;</td>\n"
				  ;
			}
			else {
				print O
"<td class=\"grid\" style=\"background-color:$color\" title=\"under-representation: $pv\" onclick=\"alert('under-representation: $pv')\">&nbsp;</td>\n";
			}
		}
		else {
			print O
"<td class=\"grid\" title=\"not significant\" >&nbsp;</td>\n";
		}
	}
	print O
"<td class=\"gridlabel\">"
	  . $name . "</td>";
	print O "\n</tr>\n";
}
print O ("\n</table>\n");
print O ("\n</body>\n");
print O ("</html>\n");

sub RGBDecToHex {
	my ( $red, $green, $blue ) = @_;
	return sprintf( "#%02lx%02lx%02lx", $red, $green, $blue );
}
