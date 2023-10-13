#updated on Jul 8, 2022
#updated on Jul 29, 2022
#updated on Mar 29, 2023
package graphSimply;

use strict;
use warnings;
use List::MoreUtils qw(any uniq);
use List::Util 'max';
use List::Util 'min';
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";

use GSAT::pathFilter;

my $infile1;
my $infile2;
my $min_path_no;
my $min_end_length;

my @seqs;
my @links;

my @seq_cp_no;
my @seq_mt_no;
my @seq_long_mt_name;

my @seq_name;
my @seq_length;
my @seq_depth;
my %ori_rev=('+'=>'-','-'=>'+');
my @bandage;
my $out_pre;

my %seqn_uniq;

sub filterGraph{
  ($infile1,$infile2,$min_path_no,$min_end_length,$out_pre)=@_;

  open (infile1,$infile1) || die "Cannot open the file: $infile1\n";
  open (infile2,$infile2) || die "Cannot open the file: $infile2\n";

  my @content=<infile1>;
  chomp(@content);
  close infile1;

  @seqs=grep {/^S/} @content;
  @links=grep {/^L/} @content;

  my @delete_nos1;

  foreach my $each_info(@seqs){
    $each_info =~ /^S\t([^\t]+)\t([^\t]+)\t/;
    my $name=$1;
    my $length=length($2);
    my $temp_dep;
    if($each_info =~ /KC:i:([0-9]+)/){
      $temp_dep=$1/$length;
      }
      elsif($each_info =~ /DP:f:([0-9\.]+)/){
      $temp_dep=$1;
      }
      else{
      die "Error in the gfa file!\n";
      }
  
    push @seq_name,$name;
    $seqn_uniq{$name}=$name;
    push @seq_length,$length;
    push @seq_depth,$temp_dep;
  }

  my $file_head=<infile2>;
  chomp(my @file2=<infile2>);
  @bandage=pathFilter::filterPath(\@file2,\@seq_name,\@seq_length,$min_end_length);

  #处理重复序列
  my $add_num1=-1;
  $out_pre=$infile1 if length($out_pre) == 0;
  $out_pre =~ s/\.gfa$//;
  open out_support,">${out_pre}.supporting";
  while ($add_num1 != 0){
    $add_num1=0;
    foreach my $this_seq(0..$#seqs){
      my $temp_target=$seq_name[$this_seq];
      my $temp_t_length=$seq_length[$this_seq];
      my $temp_t_depth=$seq_depth[$this_seq];
      my @link_hits=grep {/\t${temp_target}\t/} @links;
  
      my @before_names;
      my @before_ori;
      my @after_names;
      my @after_ori;
      my $next_stat=0;
      if(@link_hits < 4){
        $next_stat ++;
      }
  
      foreach my $this_link(@link_hits){
        my (undef,$node1,$oritation1,$node2,$oritation2,undef)=split(/\t/,$this_link);
	      if ($node1 eq $temp_target && $oritation1 eq '+'){
          push @after_names,$node2;
          push @after_ori,$oritation2;
	      }
	      elsif($node2 eq $temp_target && $oritation2 eq '-'){
          push @after_names,$node1;
          push @after_ori,$ori_rev{$oritation1};
	      }
	      elsif($node1 eq $temp_target && $oritation1 eq '-'){
          push @before_names,$node2;
          push @before_ori,$ori_rev{$oritation2};
          }
        elsif($node2 eq $temp_target && $oritation2 eq '+'){
          push @before_names,$node1;
          push @before_ori,$oritation1;
	      }
        else{
          die "Fatal error detected again!\n";
	      }
      }

      if(@before_names < 2 || @after_names <2){
        $next_stat ++;
      }
      
      #遍历连接对，方法为两个foreach
      my %true_link;
      #my %true_link_num;
    
      foreach my $l_no(0..$#before_names){
        my $left=$before_names[$l_no];
        my $l_ori=$before_ori[$l_no];
        foreach my $r_no(0..$#after_names){
          my $right=$after_names[$r_no];
          my $r_ori=$after_ori[$r_no];
          #my $path_code1='\b'."${left}".'\b[^,]*,[^,]*\b'."${temp_target}".'\b[^,]*,[^,]*\b'."${right}".'\b';#兼容GraphMapper
          my $path_code1='\b'.$seqn_uniq{$left}."\\$l_ori".','.$seqn_uniq{$temp_target}.'\+,'.$seqn_uniq{$right}."\\$r_ori";

	        #my $path_code2='\b'."${right}".'\b[^,]*,[^,]*\b'."${temp_target}".'\b[^,]*,[^,]*\b'."${left}".'\b';#兼容GraphMapper
          my $path_code2='\b'.$seqn_uniq{$right}."\\".$ori_rev{$r_ori}.','.$seqn_uniq{$temp_target}.'\-,'.$seqn_uniq{$left}."\\".$ori_rev{$l_ori};

          my @paths=grep {/${path_code1}|${path_code2}/} @bandage;
	        #print "@bandage \n";#test
          my $paths_no=@paths;
          print out_support "${temp_target}\t${temp_t_length}\t${temp_t_depth}\t${left}${l_ori}\t${right}${r_ori}\t$paths_no\n";
          next if $next_stat > 0;
          print "---Checking repeat node: $path_code1 / $path_code2 --- \n Paths: $paths_no \n";

	        if(@paths >= $min_path_no){
            if(exists $true_link{'left'}{$left.$l_ori}){
              push @{$true_link{'left'}{$left.$l_ori}},$right.$r_ori;
            }
            else{
              $true_link{'left'}{$left.$l_ori}=[$right.$r_ori];
            }
          
            if(exists $true_link{'right'}{$right.$r_ori}){
              push @{$true_link{'right'}{$right.$r_ori}},$left.$l_ori;
            }
            else{
              $true_link{'right'}{$right.$r_ori}=[$left.$l_ori];
            }
	        }
	      }
      }

      next if $next_stat > 0;
  
      foreach (0..($#before_names-1)){
        my $this_info=$before_names[$_].$before_ori[$_];
        next if not exists $true_link{'left'}{$this_info};
        my $right_info=$true_link{'left'}{$this_info}->[0];
        my $left_info=$true_link{'right'}{$right_info}->[0];
        if(($this_info eq $left_info) && (@{$true_link{'left'}{$left_info}} == 1) && (@{$true_link{'right'}{$right_info}} == 1)){
          add_copy($left_info,$this_seq,$right_info);
          print "Node has been copied: $temp_target \n";
          $add_num1=1;
        }
      }
    }
    $add_num1 =0 if $add_num1 <1;
  }

  open outfile1,">${out_pre}.simplified.gfa";

  print outfile1 join("\n",@seqs);
  print outfile1 "\n";
  print outfile1 join("\n",@links);
  close outfile1;

  print "All completed!\n";
}

sub add_copy{#before_info,copy_no,after_info
  #my ($be_name,$be_ori,$c_no,$af_name,$af_ori)=@_;
  my ($be_name,$c_no,$af_name)=@_;
  my $be_ori=substr($be_name,-1,1,'');
  my $af_ori=substr($af_name,-1,1,'');
  my $sd_name=$seq_name[$c_no];
  my $new_name=$sd_name;
  
  my @used_names=grep {$new_name eq $_} @seq_name;
  while(@used_names > 0){
    $new_name .= '_copy';
    @used_names=grep {$new_name eq $_} @seq_name;
  }
  
  push @seq_name,$new_name;
  $seqn_uniq{$new_name}=$seqn_uniq{$sd_name};
  push @seq_length,$seq_length[$c_no];
  push @seq_depth,$seq_depth[$c_no];
  
  push @seqs,$seqs[$c_no];
  $seqs[-1] =~ s/${sd_name}/${new_name}/;
  
  #对links进行改名拆解
  my $code1='\t'.$sd_name.'\t\+\t'.$af_name.'\t'."\\$af_ori";
  my $code2='\t'.$af_name.'\t'."\\".$ori_rev{$af_ori}.'\t'.$sd_name.'\t\-';
  my $code3='\t'.$be_name.'\t'."\\$be_ori".'\t'.$sd_name.'\t\+';
  my $code4='\t'.$sd_name.'\t\-\t'.$be_name.'\t'."\\".$ori_rev{$be_ori};


  my @need_rename_no=grep {$links[$_] =~ /$code1|$code2|$code3|$code4/} 0..$#links;
  die "Error: not detected enough links for adding copies! [ $code1|$code2|$code3|$code4 ]\n" if @need_rename_no < 2;
  foreach (@need_rename_no){
    $links[$_] =~ s/\b${sd_name}\b/${new_name}/g;
    }
 
  my @threads_depths=map {$seq_depth[$_]} (grep {$seq_length[$_] >= 350} 0..$#seq_length);
 
  my $thread=zhongweishu(@threads_depths)*2.2;
  print "Thread was changed into: $thread \n";
  splice @seq_mt_no;
  splice @seq_long_mt_name;
  splice @seq_cp_no;

  foreach (0..$#seqs){
  if ($seq_depth[$_] <= $thread){
    push @seq_mt_no,$_;
	
	push @seq_long_mt_name,$seq_name[$_] if $seq_length[$_] > 1000;
    }
	else{
      push @seq_cp_no,$_;
	  }
  }
}

sub zhongweishu{
  my @sorted=sort {$a <=> $b} @_;
  my $mid_num;
  if(@sorted%2){
    $mid_num=$sorted[int(@sorted/2)+1];
  }
  else{
    $mid_num=($sorted[int(@sorted/2)]+$sorted[int(@sorted/2)+1])/2;
  }
  return $mid_num;
}

1;