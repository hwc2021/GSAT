#Updated on Apr 10, 2023

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
        
        my $temp_length=length($temp_info[2]);
        $S_info{$temp_info[1]}{'len'}=$temp_length;

        if($temp_info[3] =~ /KC:i:([0-9]+)/){
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
        $L_inf[-1] =~ s/M$//;
        push @L_info,[@L_inf];
      }
    }
  }

  return (\%S_info,\@P_info,\@L_info);
}

1;