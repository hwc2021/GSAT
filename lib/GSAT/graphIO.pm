#Updated on Apr 10, 2023
#Updated on Mar 19, 2024
#Updated on Feb 15, 2025
#updated on Feb 26, 2025: fixed a bug for gfa_L
#updated on Dec 23, 2025

package graphIO;

use strict;
use warnings;
use List::MoreUtils qw(any uniq);
use List::Util 'max';
use List::Util 'min';

sub readFas{#为提高运算效率，使用引用传递数组/也可能完全用不到这个模块
  my $fas=shift @_;
  my $count=-1;
  my @fas_name;
  my @fas_seq;
  foreach (@{$fas}){
	chomp;
    if(/^>/){
      $count ++;
      $fas_name[$count]=$_;
      $fas_name[$count] =~ s/^>//;
    }
    else{
      $fas_seq[$count] .= $_;
    }
  }
  return (\@fas_name,\@fas_seq);
}

sub readGfa{#为提高运算效率，使用引用传递数组
  my $gfa=shift @_;
  my $get=shift @_;#SPL or S or SL, etc.
  my %S_info;#name-seq/depth;
  my @P_info;#[array]
  my @L_info;#[node1,or1,node2,or2,cov_length]
  chomp(@{$gfa});

  if($get =~ /S/){
    foreach (@{$gfa}){
      if(/^S/){
        my @temp_info=split(/\t/,$_,4);
        $S_info{$temp_info[1]}{'seq'}=$temp_info[2];
        if($temp_info[2] =~ /NNNNNNNNNN|nnnnnnnnnn/){
          $S_info{$temp_info[1]}{'nstat'}=1;
        }
        else{
          $S_info{$temp_info[1]}{'nstat'}=0;
        }
        
        my $temp_length=length($temp_info[2]);
        $S_info{$temp_info[1]}{'len'}=$temp_length;

        if($temp_info[3] =~ /KC:i:([0-9]+)|RC:i:([0-9]+)/){
         $S_info{$temp_info[1]}{'dep'} = $1 / $temp_length;
        }
        elsif($temp_info[3] =~ /DP:f:([0-9\.]+)/){
          $S_info{$temp_info[1]}{'dep'} = $1;
        }
        else{
          die"Error: The gfa file seems not in correct format! \n";
        }
      }
    }
  }

  if($get =~ /P/){
    @P_info=map {my @x=split(/\t/);shift @x;[@x]} (grep {/^P/} @{$gfa});
  }

  if($get =~ /L/){
    foreach (@{$gfa}){
      if(/^L/){
        my @L_inf=split(/\t/);
        shift @L_inf;
        $L_inf[-1] =~ s/M\s*$//;
        push @L_info,[@L_inf];
      }
    }
  }

  return (\%S_info,\@P_info,\@L_info);
}

sub nodeNeighbors{
  my $links=shift;
  my %node_lk;
  my %ori_r=('+','-','-','+');
  foreach my $lk(@{$links}){
    my $qo_rev=$ori_r{$lk->[1]};
    my $so_rev=$ori_r{$lk->[3]};
    if(exists $node_lk{$lk->[0]}{$lk->[1]}{'right'}){
      push @{$node_lk{$lk->[0]}{$lk->[1]}{'right'}},$lk->[2].$lk->[3];
    }
    else{
      $node_lk{$lk->[0]}{$lk->[1]}{'right'}=[$lk->[2].$lk->[3]];
    }
    if(exists $node_lk{$lk->[0]}{$qo_rev}{'left'}){
      push @{$node_lk{$lk->[0]}{$qo_rev}{'left'}},$lk->[2].$so_rev;
    }
    else{
      $node_lk{$lk->[0]}{$qo_rev}{'left'}=[$lk->[2].$so_rev];
    }
    if(exists $node_lk{$lk->[2]}{$lk->[3]}{'left'}){
      push @{$node_lk{$lk->[2]}{$lk->[3]}{'left'}},$lk->[0].$lk->[1];
    }
    else{
      $node_lk{$lk->[2]}{$lk->[3]}{'left'}=[$lk->[0].$lk->[1]];
    }
    if(exists $node_lk{$lk->[2]}{$so_rev}{'right'}){
      push @{$node_lk{$lk->[2]}{$so_rev}{'right'}},$lk->[0].$qo_rev;
    }
    else{
      $node_lk{$lk->[2]}{$so_rev}{'right'}=[$lk->[0].$qo_rev];
    }
  }
  return %node_lk;
}

1;