#Generated at Dec 4, 2025

package graphMerge;

use strict;
use warnings;
use List::MoreUtils qw(any uniq);
use List::Util 'max';
use List::Util 'min';
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";
use GSAT::graphIO;

sub rev_comp{
  my $seq=shift;
  my $out = reverse($seq);
  $out =~ tr/atcgATCG/tagcTAGC/;
  return $out;
}

#my $gfa_S;
#my $gfa_L;
sub getSeq{
  my ($seq,$ori)=@_;
  my $out;
  if($ori eq '+'){
    return $seq;
  }
  else{
    return rev_comp($seq);
  }
}

sub gmerge{
  my ($f1,$outp)=@_;
  $outp='gmerge.gfa' if not defined $outp;
  open(gfaFile,"$f1") || die "Error: Cannot open the file: $f1 \n";
  chomp(my @gfa_content=<gfaFile>);
  close gfaFile;
  
  my ($gfa_S,undef,$gfa_L)=graphIO::readGfa(\@gfa_content,'SL');
  my $gfa_lnk_k=${$gfa_L}[0]->[-1];

  #记录node相连的前后信息
  my %node_lk;
  my %ori_r=('+','-','-','+');
  foreach my $ii(0..$#{$gfa_L}){
    my $lk=$gfa_L->[$ii];
    my $qo_rev=$ori_r{$lk->[1]};
    my $so_rev=$ori_r{$lk->[3]};
    if(exists $node_lk{$lk->[0]}{$lk->[1]}{'right'}){
      push @{$node_lk{$lk->[0]}{$lk->[1]}{'right'}},$lk->[2].$lk->[3].":$ii";
    }
    else{
      $node_lk{$lk->[0]}{$lk->[1]}{'right'}=[$lk->[2].$lk->[3].":$ii"];
    }
    if(exists $node_lk{$lk->[0]}{$qo_rev}{'left'}){
      push @{$node_lk{$lk->[0]}{$qo_rev}{'left'}},$lk->[2].$so_rev.":$ii";
    }
    else{
      $node_lk{$lk->[0]}{$qo_rev}{'left'}=[$lk->[2].$so_rev.":$ii"];
    }
    if(exists $node_lk{$lk->[2]}{$lk->[3]}{'left'}){
      push @{$node_lk{$lk->[2]}{$lk->[3]}{'left'}},$lk->[0].$lk->[1].":$ii";
    }
    else{
      $node_lk{$lk->[2]}{$lk->[3]}{'left'}=[$lk->[0].$lk->[1].":$ii"];
    }
    if(exists $node_lk{$lk->[2]}{$so_rev}{'right'}){
      push @{$node_lk{$lk->[2]}{$so_rev}{'right'}},$lk->[0].$qo_rev.":$ii";
    }
    else{
      $node_lk{$lk->[2]}{$so_rev}{'right'}=[$lk->[0].$qo_rev.":$ii"];
    }
  }

  my $c1=0;
  my @m_lkn;
  my %lk2lks;
  foreach my $i(0..$#{$gfa_L}){
    my ($n1,$r1,$n2,$r2,undef)=@{$gfa_L->[$i]};
    if(exists $node_lk{$n1}{$r1}{'right'} && exists $node_lk{$n2}{$r2}{'left'} && @{$node_lk{$n1}{$r1}{'right'}} == 1 && @{$node_lk{$n2}{$r2}{'left'}} == 1 && $n1 ne $n2){
      push @m_lkn,$i;
      $c1 ++;

      my @ns_l;
      if(exists $node_lk{$n1}{$r1}{'left'}){
        @ns_l=map {(split(/:/))[-1]} @{$node_lk{$n1}{$r1}{'left'}};
      }
      my @ns_r;
      if(exists $node_lk{$n2}{$r2}{'right'}){
        @ns_r=map {(split(/:/))[-1]} @{$node_lk{$n2}{$r2}{'right'}};
      }
      $lk2lks{$i}=[uniq(@ns_l,@ns_r)];
    }
  }

  my @del_lkn;#待删除no
  #my @m_lk_rn;#每个lkn对应的相关lk的no
  foreach my $i(@m_lkn){
    my ($n1,$r1,$n2,$r2,undef)=@{$gfa_L->[$i]};
    next if $n1 eq $n2;    
    push @del_lkn,$i;
    #执行序列合并，同时检查ovl区域是否一致，不一致则警告信息
    my $s1=getSeq(${$gfa_S}{$n1}{'seq'},$r1);
    my $s2=getSeq(${$gfa_S}{$n2}{'seq'},$r2);

    my $s1o=substr($s1,${$gfa_S}{$n1}{'len'} - $gfa_lnk_k);
    my $s2o=substr($s2,0,$gfa_lnk_k);

    warn "Warnning: The overlaping sequences for this link are not consistent! [ $n1 $r1 $n2 $r2 ]\n" if $s1o ne $s2o;
    my $s3=$s1.substr($s2,$gfa_lnk_k);
    my $s3n=$n1.'_'.$n2;
    ${$gfa_S}{$s3n}{'seq'}=$s3;
    ${$gfa_S}{$s3n}{'len'}=length($s3);
    ${$gfa_S}{$s3n}{'nstat'}=${$gfa_S}{$n1}{'nstat'} + ${$gfa_S}{$n2}{'nstat'};
    ${$gfa_S}{$s3n}{'nstat'} = 1 if ${$gfa_S}{$s3n}{'nstat'} > 1;
    ${$gfa_S}{$s3n}{'dep'}= ((${$gfa_S}{$n1}{'dep'} * ${$gfa_S}{$n1}{'len'}) + (${$gfa_S}{$n2}{'dep'} * ${$gfa_S}{$n2}{'len'}))/(${$gfa_S}{$n1}{'len'} + ${$gfa_S}{$n2}{'len'});

    delete ${$gfa_S}{$n1};
    delete ${$gfa_S}{$n2};

    #foreach my $x(@{$lk2lks{$i}}){
    foreach my $x(@{$gfa_L}){
      if($x->[0] eq $n1){
        $x->[0] = $s3n;
        $x->[1] = $ori_r{ $x->[1]} if $r1 eq '-';
      }
      if($x->[2] eq $n1){
        $x->[2] = $s3n;
        $x->[3] = $ori_r{ $x->[3]} if $r1 eq '-';
      }
      if($x->[0] eq $n2){
        $x->[0] = $s3n;
        $x->[1] = $ori_r{ $x->[1]} if $r2 eq '-';
      }
      if($x->[2] eq $n2){
        $x->[2] = $s3n;
        $x->[3] = $ori_r{ $x->[3]} if $r2 eq '-';
      }
    }
  }

  my @del_lkn_sort = sort {$b <=> $a} @del_lkn;
  foreach my $y(@del_lkn_sort){
    splice(@{$gfa_L},$y,1);
  }

  my @S_content=map {"S\t".$_."\t".${$gfa_S}{$_}{'seq'}."\tDP:f:".${$gfa_S}{$_}{'dep'}} sort keys %{$gfa_S};
  my @L_content=map {"L\t".join("\t",@{$_}).'M'} @{$gfa_L};
  open(OUTF,'>'.$outp.'.merge.gfa') or die $!;
  print OUTF join("\n",(@S_content,@L_content));
  close OUTF;
  warn "[info] $c1 links have been merged!\n";
}

1;