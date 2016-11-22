#!/usr/bin/perl -w
use POSIX qw(ceil floor);
use Cwd;

my $savept = 100;  # save point
my $repno0 = 1;    # replica index 1
my $repnof = 80;   # replica index final
my $repskip =1;   # replica index final

# define derectory
# my $loaddir = "../data/"
my $dir = getcwd;
my $title = "randcopoly";
my $ind = index($dir, $title);
my $folder = substr $dir, $ind, -6;
my $loaddir = "../../../sim-".$folder."data/";

my $ratviz=0.5; # Ratio of visualization
my $cut = 0;    # Whether show cross-section

my $xlbox=20; # Box edge length
my $ylbox=20; # Box edge length
my $zlbox=20; # Box edge length
my $xshift=0;  # shift of x periodicity
my $yshift=0;  # shift of x periodicity
my $zshift=0;  # shift of x periodicity
my $nbead=40*2000;

my $radius=1;           # Radius of chain
my $ratio=1;            # Resize ratio
my $snapcount=1;        # File count for output

my $nbpbead=1;         # Number of bp per bead
my $nbpturn=10;         # Number of bp per turn
my $phidna=135;         # Rotation angle for the double helix
my $pi=3.14159265359;   # Value of pi
my $gamma=2*$pi*$nbpbead/$nbpturn; # Twist angle per bead

my $filein1;            # File with bead coordinates
my $fileout1;           # Output file for pdb

print "\n Start loading data from: ".$loaddir."...\n ";
for (my $filecount=$repno0; $filecount<=$repnof; $filecount+=$repskip) {
$filein1 =sprintf($loaddir."r%dv%d", $savept, $filecount);
$fileout1=sprintf(">snap%03d.pdb",$filecount);

open(COORD1, $filein1) || die('cannot open file:'. $!);
open(PDB1, $fileout1)  || die('cannot open file:'. $!);

my $count=1;

# Array to hold positions of backbone
   my @atomx1;
   my @atomy1;
   my @atomz1;
   my @atomx4;
   my @atomy4;
   my @atomz4;
   my @ab;
   my @viz;

   while(<COORD1>){
      my @info = split;
      $atomx1[$count] = $ratio*$info[0]+$xshift;
      $atomx1[$count] = $atomx1[$count]-floor($atomx1[$count]/$xlbox)*$xlbox;
      $atomy1[$count] = $ratio*$info[1]+$yshift;
      $atomy1[$count] = $atomy1[$count]-floor($atomy1[$count]/$ylbox)*$ylbox;
      $atomz1[$count] = $ratio*$info[2]+$zshift;
      $atomz1[$count] = $atomz1[$count]-floor($atomz1[$count]/$zlbox)*$zlbox;

      $ab[$count]=$info[3];
      
      $viz[$count]=0;
      if (($atomx1[$count] > $ratviz*$xlbox) || ($atomx1[$count] < (1-$ratviz)*$xlbox))
      {	 
	  $viz[$count]=1;
      }
      if (($atomy1[$count] > $ratviz*$ylbox) || ($atomy1[$count] < (1-$ratviz)*$ylbox))
      {
	  $viz[$count]=1;
      }
      if (($atomz1[$count] > $ratviz*$zlbox) || ($atomz1[$count] < (1-$ratviz)*$zlbox))
      { 
	  $viz[$count]=1;
      } 

      if ($cut == 1){
	  if (-3/4*$atomx1[$count] + $atomy1[$count] > 3.0 )
	  {
	      $viz[$count]=1;
	  }
	  
	  if (-3/4*$atomx1[$count] + $atomy1[$count] > 5 )
	  {
	      $viz[$count]=0;
	  }
      }

      $atomx1[$count] = $atomx1[$count]-floor($atomx1[$count]/$xlbox)*$xlbox;
      $atomy1[$count] = $atomy1[$count]-floor($atomy1[$count]/$ylbox)*$ylbox;
      $atomz1[$count] = $atomz1[$count]-floor($atomz1[$count]/$zlbox)*$zlbox;
      
      $count++;
	
   }
#   my $nbead = $count-1;          # Index of last element in atomx1

   ##############################################################
   # Assemble single PDB file
   ##############################################################

   my $atomname1 = "A1";           # Chain atom type
   my $atomname2 = "A2";           # Ribbon atom type
   my $atomname3 = "A3";           # Extra atom type
   my $atomname4 = "A4";           # Extra atom type
   my $atomname5 = "A5";           # Extra atom type
   my $atomname6 = "A6";           # Extra atom type
   my $atomname7 = "A7";           # Extra atom type
   my $atomname8 = "A8";           # Extra atom type
   my $resname = "SSN";           # Type of residue (UNKnown/Single Stranded Nucleotide)
   my $chain = "A";               # Chain identifier
   my $resnum = "1";
   my $numresidues = $nbead;
   my $descrip = "Pseudo atom representation of DNA";
   my $chemicalname = "Body and ribbon spatial coordinates";
   
   # Het Header info
   printf PDB1 "HET    %3s  %1s%4d   %5d     %-38s\n",$resname,$chain,$resnum, $numresidues,$descrip ;
   printf PDB1 "HETNAM     %3s %-50s\n",$resname, $chemicalname;
   printf PDB1 "FORMUL  1   %3s    C20 N20 P21\n",$resname;

   $count=1;

   $x=0;
   $y=0;
   $z=0;
   printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n" ,
   ,$count,$atomname3,$resname,$chain,
   ,$x,$y,$z,
   ,1.00,1.00;
   $count++;
   $x=$xlbox;
   $y=0;
   $z=0;
   printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n" ,
   ,$count,$atomname3,$resname,$chain,
   ,$x,$y,$z,
   ,1.00,1.00;
   $count++;
   $x=0;
   $y=$ylbox;
   $z=0;
   printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n" ,
   ,$count,$atomname3,$resname,$chain,
   ,$x,$y,$z,
   ,1.00,1.00;
   $count++;
   $x=0;
   $y=0;
   $z=$zlbox;
   printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n" ,
   ,$count,$atomname3,$resname,$chain,
   ,$x,$y,$z,
   ,1.00,1.00;
   $count++;
   $x=$xlbox;
   $y=$ylbox;
   $z=0;
   printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n" ,
   ,$count,$atomname3,$resname,$chain,
   ,$x,$y,$z,
   ,1.00,1.00;
   $count++;
   $x=$xlbox;
   $y=0;
   $z=$zlbox;
   printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n" ,
   ,$count,$atomname3,$resname,$chain,
   ,$x,$y,$z,
   ,1.00,1.00;
   $count++;
   $x=0;
   $y=$ylbox;
   $z=$zlbox;
   printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n" ,
   ,$count,$atomname3,$resname,$chain,
   ,$x,$y,$z,
   ,1.00,1.00;
   $count++;
   $x=$xlbox;
   $y=$ylbox;
   $z=$zlbox;
   printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n" ,
   ,$count,$atomname3,$resname,$chain,
   ,$x,$y,$z,
   ,1.00,1.00;
   $count++;

   printf PDB1 "CONECT%5d%5d%5d%5d\n",1,2,3,4;
   printf PDB1 "CONECT%5d%5d%5d%5d\n",2,1,5,6;
   printf PDB1 "CONECT%5d%5d%5d%5d\n",3,1,5,7;
   printf PDB1 "CONECT%5d%5d%5d%5d\n",4,1,6,7;
   printf PDB1 "CONECT%5d%5d%5d%5d\n",5,2,3,8;
   printf PDB1 "CONECT%5d%5d%5d%5d\n",6,2,4,8;
   printf PDB1 "CONECT%5d%5d%5d%5d\n",7,3,4,8;
   printf PDB1 "CONECT%5d%5d%5d%5d\n",8,5,6,7;   

   $ii=1;
   while ($ii <= $nbead) {
       if ($ab[$ii]==1 && $viz[$ii]==1) {
	   printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n" ,
	       ,$count,$atomname1,$resname,$chain,
	       ,$atomx1[$ii],
	       ,$atomy1[$ii],
	       ,$atomz1[$ii],1.00,1.00;
#       printf PDB1 "CONECT%10d%10d\n",$ii+8,$ii+9;
       }
       if ($ab[$ii]==0 && $viz[$ii]==1) {
	   printf PDB1 "ATOM%7d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f           C\n",
	       $count,$atomname2,$resname,$chain,
	       ,$atomx1[$ii],
	       ,$atomy1[$ii],
	       ,$atomz1[$ii],1.00,1.00;
#       printf PDB1 "CONECT%10d%10d\n",$ii+8,$ii+9;
       }
#       printf PDB1 "CONECT%10d%10d\n",$ii-1+8,$ii+8;

       $ii++;
       $count++;
   }

   printf PDB1 "END";
   # Clean up and close files

   close(PDB1);
   close(COORD1);
}
