#!perl
#Generated at Dec 19, 2025
#updated at Dec 20, 2025: enable the statistics for copy number of repetitive contigs.
package graphStat;

use strict;
use warnings;
#use List::MoreUtils qw(any uniq);
use List::Util qw(max min sum);
#use List::Util 'min';
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";
use GSAT::graphIO;

sub mmgStat{
  my ($gfa,$outp,$est_stat)=@_;
  $outp=$gfa if not defined $outp;
  open(gfaFile,"$gfa") || die "Error: Cannot open the file: $gfa \n";
  chomp(my @gfa_content=<gfaFile>);
  close gfaFile;
  my ($gfa_S,undef,$gfa_L)=graphIO::readGfa(\@gfa_content,'SL');

  #记录node相连的前后信息
  my %node_lk=graphIO::nodeNeighbors($gfa_L);
  
  my @summary=(0) x 3;
  my @line;
  my @cn_r;
  my @cl_r;
  open outFile,">${outp}.mmgStat";
  my $h_str="Node\tNum\tNodeLength";
  do {
    $h_str .= "\tCopyNum\tTotalLength";
    @summary=(0) x 5;
    warn "Please note that this is an experimental feature and the resulting data should be used with caution.\n";
  } if $est_stat > 0;

  print outFile $h_str."\n";
  foreach my $ctg(sort keys %node_lk){
    @line=($ctg,'1',${$gfa_S}{$ctg}{'len'});

    if($est_stat > 0){
      my $lc=$node_lk{$ctg}{'+'}{'left'};
      my $rc=$node_lk{$ctg}{'+'}{'right'};
      my $c1=@{$lc};
      my $c2=@{$rc};
      
      if($c1 == 1 || $c2 == 1){
        my $lc_n=$lc->[0];
        my $lc_o=substr($lc_n,-1,1,'');
        my $rc_n=$rc->[0];
        my $rc_o=substr($rc_n,-1,1,'');
        die "Error: The estimation for copies of repetitive contigs is only available for merged graph!\n" if ($lc_n ne $ctg && $c1 == 1 && @{$node_lk{$lc_n}{$lc_o}{'right'}} == 1) || ($ctg ne $rc_n && $c2 == 1 && @{$node_lk{$rc_n}{$rc_o}{'left'}} == 1);       
      }

      my $cn_c = $c1 * $c2;

      push @line,$cn_c;
      push @line,$cn_c * ${$gfa_S}{$ctg}{'len'};
      if($cn_c > 1){
        push @cn_r,$cn_c;
        push @cl_r,$line[-1];
      }
    }

    my @tmp_sum=map {$summary[$_-1]+$line[$_]} 1..$#line;
    @summary=@tmp_sum;

    print outFile join("\t",@line)."\n";
  }

  print outFile join("\t",('all',@summary))."\n";
  my @sort_cn=sort {$a <=> $b} @cn_r;
  my $cn_r_sum=sum(@sort_cn);
  my $cn_r_num=@sort_cn;
  my $cl_r_sum=sum(@cl_r);
  print outFile "\n#\tRepContigNum\tminCN\tmaxCN\ttotalCN\ttotalCopyLen\n#\t".join("\t",($cn_r_num,$sort_cn[0],$sort_cn[-1],$cn_r_sum,$cl_r_sum))."\n" if $est_stat > 0;
  close outFile;
}

1;