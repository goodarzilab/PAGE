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


my $pvaluematrixfile  = undef;
my $expfile           = undef;
my $datafile          = undef;
my $colmap            = "$scriptdir/HEATMAPS/cmap_dens.txt";
my $cluster           = undef;
my $sortrowsbyphase   = 1;
my $max_p             = undef;
my $minmax_lp         = 3;
my $quantized         = 1;
my $min               = undef;
my $max               = undef;
my $xsize             = undef;
my $xscale            = 5 ;
my $yscale            = 75 ;
my $scalefont         = 9 ;
my $h                 = undef ;
my $w                 = undef ;
my $draw_sample_heatmap = "false" ;
my $order             = 0 ;

  if (@ARGV == 0) {
 die "Usage: perl mi_go_draw_matrix.pl  --pvaluematrixfile=FILE --expfile=FILE --max_p=P\n";
}

GetOptions ('pvaluematrixfile=s'  => \$pvaluematrixfile,
            'cluster=s'           => \$cluster,
			'expfile=s'           => \$expfile,
			'datafile=s'          => \$datafile,
			'quantized=s'         => \$quantized,
			'colmap=s'            => \$colmap,
	    'minmax_lp=s'         => \$minmax_lp,
	    'min=s'               => \$min,
	    'max=s'               => \$max,
	    'xsize=s'             => \$xsize,
	    'xscale=s'            => \$xscale,
	    'yscale=s'            => \$yscale,
	    'scalefont=s'         => \$scalefont,
	    'h=s'                 => \$h,
	    'w=s'                 => \$w,
	    'max_p=s'             => \$max_p,
	    'draw_sample_heatmap=s' => \$draw_sample_heatmap,
	    'order=s'             => \$order);

#
# creating the summary file
#
my $file = substr($expfile, rindex($expfile, '/')) ;
my $dir = $expfile."_PAGE/" ;
my $ta = Table->new;

print "Reading matrix ... ";

#
#  read in the matrix file
#
$ta->loadFile($pvaluematrixfile);

# get an 2D array
my $a_ref_M      = $ta->getArray();

# header
my $a_ref_H      = shift @$a_ref_M; shift @$a_ref_H;
if (!defined($max_p)) {
  $max_p = 0.05 / @$a_ref_H;
}
print "Done.\n";

if ($cluster>@$a_ref_M){
    $cluster = scalar(@$a_ref_M) ;
}

if (defined($cluster) && (@$a_ref_M > 2)) {

 print "Cluster rows .. ";

 my $ac = AggloClust->new;

 my @dist = ();
 my $n    = @$a_ref_M;
 for (my $i=0; $i<$n-1; $i++) {
   $dist[$i][$i] = 0;
   for (my $j=$i+1; $j<$n; $j++) {
     my @a1 = @{ $a_ref_M->[$i] }; shift @a1;
     my @a2 = @{ $a_ref_M->[$j] }; shift @a2;
     $dist[$i][$j] = 1 - Sets::pearson(\@a1, \@a2);
     $dist[$j][$i] = $dist[$i][$j];
   }
 }

 $ac->setDistanceMatrix(\@dist);
 $ac->setMaxNbClusters($cluster);
 my $a_ref_c = $ac->agglomerate_using_avg_linkage();

 my @NEWMAT = ();
 foreach my $c (@$a_ref_c) {
   print join(" ", @$c); print "\n";
   foreach my $i (@$c) {
     push @NEWMAT, $a_ref_M->[$i];
   }
 }
 $a_ref_M = \@NEWMAT;

 print "Done.";
}

#
# load color map
#
my $A_REF_COLMAP = undef;
if (defined($colmap)) {

 $ta->setDelim(" ");
 $ta->loadFile($colmap);
 $A_REF_COLMAP = $ta->getArray();
 $ta->setDelim('\t');
}


#
#  START DRAWING
#
#
my %CENTROIDS ;
my @samples ;
my @clusters ;
my $nbsamples ;
my $minimum ;
my $maximum ;
my %ARRAYS ;

#
# left and top margins
#
my $xbase        = 60;
my $ybase        = 40;
$ybase           = 100 if ($quantized==0) ;

# height of each entry in the matrix
$h            = 30 if ! (defined $h) ;

# size of the EPS image to be generated
my $ysize        = $ybase + $h * scalar(@$a_ref_M) + 100;
my $h_h = 10 ;
if ($quantized == 1 and defined $datafile and $draw_sample_heatmap eq "true")
{
    my $todo = "perl $scriptdir/draw_clusters_and_gocats.pl --expfile=$expfile --data=$datafile --pvaluematrixfile=$pvaluematrixfile" ;
    system("$todo") == 0 or die "system failed: $?";
    
    open CL, "< $dir/sample_centroids.txt" or die;
    print "$dir/sample_centroids.txt" ;
    
    chomp ;
    
    $minimum = 1000;
    $maximum = 0;
    my $l = <CL> ;
    chomp ($l) ;
    @samples = split(/\t/, $l) ;
    shift (@samples) ;
    while(<CL>)
    {
	chomp ;
	my ($cluster, @a) = split(/\t/, $_) ;
	$nbsamples = $#a+1 ;
	for (my $i=0; $i<$nbsamples; $i++) 
	{
	    $CENTROIDS{$cluster}{$samples[$i]} = $a[$i] ;
	    push(@{$ARRAYS{$samples[$i]}}, $a[$i]) ;
	}
	my $m1 = Sets::minInArray(\@a) ;
	my $m2 = Sets::maxInArray(\@a) ;
	$minimum = $m1 if ($m1<$minimum) ;
	$maximum = $m2 if ($m2>$maximum) ;
	push (@clusters, $cluster) ;
    }
    $ysize += $h_h * $nbsamples ;

    for (my $j=0; $j<$nbsamples; $j++){
	my $average = Sets::average($ARRAYS{$samples[$j]}) ;
	my $std = Sets::stddev($ARRAYS{$samples[$j]}) ;
	for (my $i=0 ; $i<=$#clusters ; $i++){
	    $CENTROIDS{$clusters[$i]}{$samples[$j]} = ($CENTROIDS{$clusters[$i]}{$samples[$j]}-$average)/$std;
	}
    }

    my @CEN ;
    my @NEWMAT = ();
    my @HEADER = () ;
    my @enr ;
    for (my $i=0 ; $i<@$a_ref_H ; $i++){
	my $c = $a_ref_H->[$i] ;
	my $delta = $CENTROIDS{$c}{$samples[1]}-$CENTROIDS{$c}{$samples[0]} ;
	$enr[$i]->{delta} = $delta ;
	$enr[$i]->{index} = $i+1 ;
    }
    my @enr = sort {$b->{delta} <=> $a->{delta}} (@enr) ;

    for (my $i=0 ; $i<@$a_ref_M ; $i++){
	$NEWMAT[$i][0] = $a_ref_M->[$i][0] ;
    }
    my $cnt = 1 ;
    foreach my $e (@enr){
        my $j = $e->{index} ;
	print $j, "\t" ;
	for (my $i=0 ; $i<@$a_ref_M ; $i++){
	    $NEWMAT[$i][$cnt] = $a_ref_M->[$i][$j] ;
	}
	$HEADER[$cnt-1] = $a_ref_H->[$j-1] ;
	$cnt++ ;
    }
    $a_ref_M = \@NEWMAT ;
    $a_ref_H = \@HEADER ;
}

if ($order ==1){
    my @NEWMAT = ();
    my @enr ;
    for (my $i=0 ; $i<@$a_ref_M ; $i++){
	my $r = $a_ref_M->[$i] ;
	my $min = 0 ;
	my $pos = 1 ;
	for (my $j=1 ; $j<@$r ; $j++){
	    my $lpo = -1*$r->[$j] ;
	    if ($lpo<$min and $r->[$j]>0){
		$min = $lpo ;
		$pos = $j ;
	    }
	}
	$enr[$i]->{ind} = $i ;
	$enr[$i]->{pos} = $pos ;
    }
    my @enr = sort {$b->{pos} <=> $a->{pos}} (@enr) ;
    foreach my $e (@enr){
	my $i = $e->{ind} ;
	push @NEWMAT, $a_ref_M->[$i];
    }
    $a_ref_M = \@NEWMAT ;
}


$xsize           = 1000 if ! (defined $xsize);

# width of each entry in the matrix
$w = int( 0.5 + ((1 * $xsize / 2.0) + 5)  / @$a_ref_H )-1 if (! defined($w));
print "Start drawing\n";

my $p = new PostScript::Simple(xsize     => $xsize,
			       ysize     => $ysize,
                              colour    => 1,
                              eps       => 1,
                              units     => "pt");

$p->setlinewidth(0.5);

my $f = Sets::filename($expfile) ;
#
# show header (cluster indices)
#
if ($quantized == 1) {

    if (defined $datafile and $draw_sample_heatmap eq "true") {
	my @color = ();
	#my $cmap = "$scriptdir/HEATMAPS/c.txt" ;
	my $cmap = undef ;
	my $REF_COLMAP = undef;
	for (my $j=0; $j<$nbsamples; $j++) 
	{
	    for (my $i=0 ; $i<@$a_ref_H ; $i++)
	    {
		if (defined($cmap)) 
		{
		    $ta->setDelim(" ");
		    $ta->loadFile($cmap);
		    $REF_COLMAP = $ta->getArray();
		    $ta->setDelim('\t');
		}
		my @color = () ;
		if ($CENTROIDS{$a_ref_H->[$i]}{$samples[$j]}>0){
		    @color = Sets::interp_general( $CENTROIDS{$a_ref_H->[$i]}{$samples[$j]},[255, 0, 0], [255, 255, 0], 0, 2);
		}
		else{
		    @color = Sets::interp_general( $CENTROIDS{$a_ref_H->[$i]}{$samples[$j]},[0, 0, 0], [255, 0, 0], -2, 0)
		}
		#my @color = Sets::interp_general( $CENTROIDS{$a_ref_H->[$i]}{$samples[$j]},[0, 255, 0], [255, 0, 0], -2, 2);
		#my @color = Sets::interp_from_matlab_colormap( $CENTROIDS{$a_ref_H->[$i]}{$samples[$j]},$REF_COLMAP, -2, 2);
		
		$p->setcolour(@color);
		$p->box({filled => 1},
			$xbase + $i * $w,      $ysize - ($ybase + $j*$h_h) ,
			$xbase + $i * $w + $w, $ysize - ($ybase + ($j*$h_h+$h_h)));
		
	    }
	    $p->setfont("Arial", 8);
	    $p->setcolour("black");
	    $p->text({align => "left", rotate => 0}, $xbase + scalar(@$a_ref_H) * $w + 10, $ysize - ($ybase + $j*$h_h+$h_h/2+4), $samples[$j]);
	    print $samples[$j], "\n" ;
	}

	drawScale($xscale+20, $ybase , -2, 2, 50, $p, $xsize, $ysize, $cmap, $REF_COLMAP, "", "", 8, 0.3, 20);
	$ybase = $ybase - ($nbsamples*$h_h) + $h_h*6 + 10;
    }

    # drawing header
    for (my $j=0; $j<@$a_ref_H; $j++) {
	$p->setcolour("black");
	$p->setfont("Arial", $w);
	$p->text( { rotate => 90 }, $xbase + $j * $w + 3*$w/4, $ysize - ($ybase-3), $a_ref_H->[$j]);
    }
    
} else {
    
  # graphical header for continuous data
  
  my $min_i1 =  1000000;
  my $max_i2 = -1000000;
  my @bins   = ();
  foreach my $c (@$a_ref_H) {
    my ($i1, $i2) = $c =~ /\[(.+?)\ (.+?)\]/;    
    $min_i1 = $i1 if ($i1 < $min_i1);
    $max_i2 = $i2 if ($i2 > $max_i2);
    my @a_tmp = ($i1, $i2); push @bins, \@a_tmp;
  }
  
  my $th = $max_i2 - $min_i1;
  my $hi = $h * 1.5;

  my $j = 0;
  foreach my $c (@$a_ref_H) {

    my $h1 = $hi * ($bins[$j]->[0] - $min_i1 ) / $th ;
    my $h2 = $hi * ($bins[$j]->[1] - $min_i1 ) / $th ;

    $p->setcolour("black");    
    $p->box({filled => 1}, 
	    $xbase + $j * $w,      $ysize - ($ybase - 5 - $hi) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 + 0 ));

   
    $p->setcolour("white");
    $p->line( $xbase + $j * $w + $w,      $ysize - ($ybase - 5 - $hi) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 + 0 ) );
    
    $p->setcolour("red");    
    $p->box({filled => 1}, 
	    $xbase + $j * $w,      $ysize - ($ybase - 5 - $h2) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 - $h1 ));
    
    my $cc = $c; $cc =~ s/^.+\]//; #print "$c $cc\n";
    
    $p->setcolour("black");
    $p->setfont("Arial", 6);
    $p->text( { rotate => 90 }, $xbase + $j * $w + 3*$w/4, $ysize - ($ybase-5 - $hi), $cc);  


    $j ++;
  }

  $p->setlinewidth(0.5);

  $j = 0;
   foreach my $c (@$a_ref_H) {
    
    my $h1 = $hi * ($bins[$j]->[0] - $min_i1 ) / $th ;
    my $h2 = $hi * ($bins[$j]->[1] - $min_i1 ) / $th ;
   
    $p->setcolour("white");
    $p->line( $xbase + $j * $w + $w,      $ysize - ($ybase - 5 - $hi) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 + 0 ) );
    
    $j ++;
  }

  

  print "$max_i2\t$min_i1\n";

  $p->setfont("Arial", 8);

  $p->setcolour("black");    
  $p->text( {align => 'right'}, $xbase - 1 , $ysize - ($ybase - 5 - $hi + 4), $max_i2);
  $p->text( {align => 'right'}, $xbase - 1, $ysize - ($ybase - 5 - 0   + 1), $min_i1);

}



$p->setcolour("black");
$p->setfont("Arial", 8);


#
# draw (i,j) p-value matrix itself
#

# set min max p-values
$max =  $minmax_lp if !(defined $max);
$min = -$minmax_lp if !(defined $min);
my @col = ();
my @go ;
for (my $i=0; $i<@$a_ref_M; $i++) {
 # get row
 my $r  = $a_ref_M->[$i];

 # get GO description
 my $go = shift @$r;

 $go =~ s/^(.+?)\ //;
 $go = "$go, $1" if ($go ne $1) ;
 
 print "$go\n";
 push(@go, $go) ;

 # go through entries in current row
 for (my $j=0; $j<@$r; $j++) {

   # get log pvalues for over and under rep
   my $lp = $r->[$j] ;

   # get a single value, pos if over-rep, neg if under-rep

   # create appropriate color
   my @col = ();  #$colmap = undef;
   if (!defined($colmap)) {
       if ($lp<0)
       {
	   @col = Sets::interp_general( $lp, [0, 255, 0], [0, 0, 0], $min, 0);
       }
       else
       {
	   @col = Sets::interp_general( $lp, [0, 0, 0], [255, 0, 0], 0, $max);
       }
   } else {
     @col = Sets::interp_from_matlab_colormap( $lp, $A_REF_COLMAP, $min, $max);
   }

   # draw the matrix entry with appropriate color
   $p->setcolour(@col);
   $p->box({filled => 1},
           $xbase + $j * $w,      $ysize - ($ybase + $i*$h) ,
           $xbase + $j * $w + $w, $ysize - ($ybase + ($i*$h+$h)));

 }

 $p->setfont("Arial", 14);
 $p->setcolour("black");
 $p->text({align => "left", rotate => 0}, $xbase + @$r * $w + 10, $ysize - ($ybase + $i*$h+$h/2+4), $go);


}



#
# draw scale bar
#
drawScale($xscale, $ysize / 2 - $yscale , $min, $max, 50, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, "Over-representation", "Under-representation", $scalefont, 2, 20);

my $outeps = "$expfile"."\_PAGE/$f".".summary.eps" ;
my $outpdf = "$expfile"."\_PAGE/$f".".summary.pdf" ;
my $ps2pdf = 1;


# output EPS file
print "Outputing EPS file $outeps\n";
$p->output("$outeps");

# convert to PDF
print "Convert to PDF $outpdf\n";
system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");

print "Finished.\n";
exit(0);


sub drawScale {
 my ($x, $y, $min, $max, $res, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, $upper_text, $lower_text, $scalefont, $h, $w) = @_;

 my $sep = 0;

 $p->setcolour("black");
 $p->setfont("Arial", 10);
 $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y - 3), $max);

 $p->setfont("Arial", $scalefont);

 $p->text({align => "left", rotate => 90}, $x+$w/2+5, $ysize - ($y - 16), $upper_text);


 my $t = $max;


 for (my $i=0; $i<=$res; $i++) {

   my @col = () ;
   if (!defined($colmap)) {
       if ($i>$res/2)
       {
	   @col = Sets::interp_general( $t, [0, 0, 0], [255, 0, 0], $min, 0);
       }
       else
       {
	   @col = Sets::interp_general( $t, [255, 0, 0], [255, 255, 0], 0, $max);
       }
   } else {
     @col = Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max);
   }

   $p->setcolour( @col );
   $p->box({filled => 1}, $x, $ysize - ($y + $sep + $i*$h) , $x+$w, $ysize - ($y + $sep + $i*$h + $h));
   $t -= ($max - $min) / $res;
 }

 $p->setcolour( 'black' );
 $p->setfont("Arial", 10);
 $p->text({align => "center"}, $x+$w/2+1, $ysize - ($y + $sep + $res*$h + 11), $min);
 $p->setfont("Arial", $scalefont);
 $p->text({align => "right", rotate => 90}, $x+$w/2+5, $ysize - ($y + $sep + $res*$h + 11 + 10), $lower_text);


}





sub interp {
 my ($r, $min, $max) = @_;

 if ($r < $min) {
   $r = $min;
 }

 if ($r > $max) {
   $r = $max;
 }

 #  1 --> red           => "0.8  0    0",
 #  0 --> blue          => "0    0    0.8",

 #  0.8 = $max . a + b
 #  0   = $min . a + b

 my $a1 = 0.8 / ($max - $min);
 my $b1 = - $a1 * $min;

 #  0   = $max . a + b
 #  0.8 = $min . a + b

 my $a3 = -0.8 / ($max - $min);
 my $b3 = -$a3 * $max;

 my $c1 = int( 0.5 + 256 * ($a1 * $r + $b1) );
 my $c2 = 0;
 my $c3 = int( 0.5 + 256 * ($a3 * $r + $b3) );


 return ($c1, $c2, $c3);

}
