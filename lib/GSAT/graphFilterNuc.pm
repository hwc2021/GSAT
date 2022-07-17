#Written by He W. (nongke2@163.com) at AGIS.
#updated on Jul 8.

package graphFilterNuc;

use strict;
use warnings;
use List::MoreUtils qw(any uniq);
use List::Util 'max';
use List::Util 'min';

my @seq_name;
my @seq_length;
my @seq_depth;
my @delete_nos1;
my @seqs;
my @links;

sub filterNuc{
  my $infile1=shift;
  my $depth500=shift;
  my $depth1000=shift;
  my $rm_sep=shift;
  my $outP=shift;

  open (infile1,$infile1) || die "Cannot open the file: $!\n";

  my @content=<infile1>;
  chomp(@content);
  close infile1;

  @seqs=grep {/^S/} @content;
  @links=grep {/^L/} @content;

  my $temp_count = -1;
  foreach (@seqs){
    my (undef,$name,$the_seq,$temp)=split(/\t/,$_,4);
    my $length=length($the_seq);
    my $temp_dep;
    if($temp =~ /DP:f:([0-9\.]+)/){
      $temp_dep=$1;
    }
    elsif($temp =~ /KC:i:([0-9]+)/){
      $temp_dep=$1/$length;
    }
    else{
      die"Error: Unspecified depth information in the gfa file!\n";
    }
  
    push @seq_name,$name;
    push @seq_length,$length;
    push @seq_depth,$temp_dep;
  
    $temp_count++;
    if(($length>500 && $temp_dep < $depth500) || ($length>1000 && $temp_dep < $depth1000)){
      push @delete_nos1,$temp_count;
    }
  }

  #处理残留的nuc相关nodes
  remove_nodes(@delete_nos1);

  #处理cp相关nodes
  if($rm_sep eq 'T'){
    my @delete_nos2;
    foreach my $this_seq(0..$#seqs){
      my $before=0;
      my $after=0;
      my $temp_target=$seq_name[$this_seq];
      my @temp_hits=grep {/\b${temp_target}\b/} @links;

      foreach my $each(@temp_hits){
        my (undef,$node1,$oritation1,$node2,$oritation2,undef)=split(/\t/,$each);
  	    if (($node1 eq $temp_target && $oritation1 eq '+') || ($node2 eq $temp_target && $oritation2 eq '-')){
          $after++;
  	    }
  	    elsif(($node1 eq $temp_target && $oritation1 eq '-') || ($node2 eq $temp_target && $oritation2 eq '+')){
          $before++;
        }
        else{
          die "Fatal error detected!\n";
	      }
      }
	
      if($after + $before == 0){
        push @delete_nos2,$this_seq;
        print "Will be deleted (seperate node): ".$seq_name[$this_seq]."\n";
	      next;
      }
      else{
        ;
	    }
    }
    remove_nodes(@delete_nos2);#注意这里接受的是数组的下标
  }

  $infile1 =~ s/\.gfa$//;
  $outP=$infile1 if not defined($outP) || length($outP) == 0;
  open outfile1,">${outP}.filtered.gfa";

  print outfile1 join("\n",@seqs);
  print outfile1 "\n";
  print outfile1 join("\n",@links);
  close outfile1;

  print "All completed!\n";
 
}

#只接受数组下标的列表
sub remove_nodes{
  my @sorted_nos=sort {$b <=> $a} @_;
	
  foreach my $this_no(@sorted_nos){
	 my $temp_offset=$this_no;
	 my $temp_name=$seq_name[$this_no];
	 
   print "Deleting the node:".$seq_name[$this_no]."\n";
   splice(@seqs,$temp_offset,1);
	 splice(@seq_name,$temp_offset,1);
	 splice(@seq_length,$temp_offset,1);
	 splice(@seq_depth,$temp_offset,1);
	 
	 my @new_links=grep {$_ !~ /\b${temp_name}\b/} @links;
	 @links=@new_links;
  }
}

1;