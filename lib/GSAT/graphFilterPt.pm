#!perl
#updated on Jul 15, 2022
#fixed several bugs in filterPt and rm_bb_cp functions on Mar 27, 2023
#fixed a bug when remove the edge-bubble; updated the messages. Oct 13, 2023

package graphFilterPt;

use strict;
use warnings;
use File::Basename;
use List::MoreUtils qw(any uniq);
use List::Util 'max';
use List::Util 'min';

my @seqs;
my @links;

my @seq_cp_no;
my @seq_mt_no;
my @seq_long_mt_name;

my @seq_name;
my @seq_length;
my @seq_depth;

sub filterPt{
  our ($infile1,$infile2,$min_path_no,$min_end_length,$rm_bubble_cp,$remain_bubble,$outP)=@_;
  open (infile1,$infile1) || die "Cannot open the file: $infile1\n";
  open (infile2,$infile2) || die "Cannot open the file: $infile2\n";

  my @content=<infile1>;
  chomp(@content);
  close infile1;

  @seqs=grep {/^S/} @content;
  @links=grep {/^L/} @content;
  my %ori_rev=('+'=>'-','-'=>'+');

  my $temp_count = -1;
  foreach my $each_info(@seqs){
    #my (undef,$name,undef,$length,$temp)=split(/\t/);
    $each_info =~ /^S\t([^\t]+)\t([^\t]+)\t/;
    my $name=$1;
    my $length=length($2);
    my $temp_dep;
    if($each_info =~ /[KR]C:i:([0-9]+)/){
      $temp_dep=$1/$length;
    }
    elsif($each_info =~ /DP:f:([0-9\.]+)/){
      $temp_dep=$1;
    }
    else{
      die "Error in the gfa file!\n";
    }
  
    push @seq_name,$name;
    push @seq_length,$length;
    push @seq_depth,$temp_dep;
  
    $temp_count++;
  }

  remove_nodes();
  my @bandage;#最终结果赋值给@bandage，在此处直接对末端进行过滤即可。
  my $file_head=<infile2>;
  foreach (<infile2>){
    chomp;
    my @p_info=split(/\t/);
    my $this_path=$p_info[4];
    my $f_stat=-1;
    my $l_stat=-1;
    $this_path =~ s/^\([^\(\)]+\)|\([^\(\)]+\)$//g;
    #rm_edge_bubble #循环去除边界bubble，如何处理边界比对结果：打满比对范围即可
    while($this_path =~ m/^\{/){
      $this_path =~ s/^\{[^,\}]+\},*//;
      $f_stat=1;
    }
    while($this_path =~ m/\}$/){
      $this_path =~ s/,*\{[^,\}]+\}$//;
      $l_stat=1;
    }

    #rm_bubble_cp
    if($rm_bubble_cp eq 'T'){
      $this_path = rm_bb_cp($this_path);
    }
    else{
      ;
    }

    #remain_bubble
    next if $remain_bubble ne 'T' && $this_path =~ /\{[^@]/;
  
    #去除匹配过短的末端
    $this_path =~ s/\{\@[^\}]+\},*|\([^\)]+\),*|\[[^\]]+\],*//g;#去除距离等附加信息
	if(length($this_path) == 0){
	  warn "Warning: A path was fitered due to zero length after removing distance and other information. \n";
	  next;
	}
	
    my @this_paths=split(/,/,$this_path);
    if($f_stat == -1){
      my $f_ctg=$this_paths[0];
      $f_ctg =~ s/[\+\-]$//;
      die "Error: wrong ctg names were detected: first ctg($f_ctg)\n" if length($f_ctg) ==0;
      my ($f_length)=map {$seq_length[$_]} (grep {$seq_name[$_] eq $f_ctg} 0..$#seq_length);
      if(length($f_length) ==0){
        warn "Warning: Cannot find the first ctg($f_ctg) of paths from the gfa file! This path will be skipped!\n";
        next;
      }
      my ($fs,$fe)=split(/\-/,$p_info[5]);
      my $f_ali_length=$fe-$fs+1;
      if(($f_ali_length < $f_length) && ($f_ali_length < $min_end_length)){
        $f_stat=0;
      }
    }

    if($l_stat == -1){
      my $l_ctg=$this_paths[-1];
      $l_ctg =~ s/[\+\-]$//;
      die "Error: wrong ctg names were detected: last ctg($l_ctg)\n" if length($l_ctg)==0;
      my ($l_length)=map {$seq_length[$_]} (grep {$seq_name[$_] eq $l_ctg} 0..$#seq_length);
      if(length($l_length==0)){
        warn "Warning: Cannot find the last ctg($l_ctg) of paths from the gfa file! This path will be skipped!\n";
        next;
      }
      my ($ls,$le)=split(/\-/,$p_info[6]);
      my $l_ali_length=$le-$ls+1;

      if(($l_ali_length < $l_length) && ($l_ali_length < $min_end_length)){
        $l_stat=0;
      }
    }
  
    if($f_stat == 0){
      shift @this_paths;
    }
    if($l_stat == 0){
      pop @this_paths;
    }
    push @bandage,join(',',@this_paths);
  }
  close infile2;
  remove_nodes();
  my $del_num1=-1;
  while($del_num1 != 0){

    my @delete_nos2;
    foreach my $this_seq(@seq_cp_no){
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
	
      if($after * $before == 0){
        push @delete_nos2,$this_seq;
        print "Will be deleted (seperate node): ".$seq_name[$this_seq]."\n";
	      next;
      }
	    else{
        my @linked=grep {/\b${temp_target}\b/} @bandage;
        my $conn_with_long=0;
        foreach my $link_line(@linked){
          foreach my $long_name(@seq_long_mt_name){
            if($link_line =~ /\b${long_name}\b/){
              $conn_with_long ++;
              print "Detected long mt contigs for cp-segs: ".$seq_name[$this_seq]." -- ${long_name}\n";
            }
          }
        }
	  
	      if ($conn_with_long < $min_path_no){
          push @delete_nos2,$this_seq;
          print "Will be deleted (only MT node): ".$seq_name[$this_seq]."\n";
        }
	    }

    }

    remove_nodes(@delete_nos2);#注意这里接受的是数组的下标
    $del_num1 = @delete_nos2;
  }

  $outP=basename($infile1) if length($outP) == 0;

  open outfile1,">${outP}.filtered.gfa";

  print outfile1 join("\n",@seqs);
  print outfile1 "\n";
  print outfile1 join("\n",@links);
  close outfile1;

  print "All completed!\n";
}

sub zhongweishu {
  my @sorted=sort {$a <=> $b} @_;
  my $middle_id=int(@sorted/2);
  return $sorted[$middle_id];
}

#只接受数组下标的列表
sub remove_nodes{
  my @sorted_nos=sort {$b <=> $a} @_;
	
  foreach my $this_no(@sorted_nos){
	   my $temp_offset=$this_no;
	   my $temp_name=$seq_name[$this_no];
	 
     print "Deleting the node:".$seq_name[$this_no]." [length(".$seq_length[$this_no].") depth(".$seq_depth[$this_no].")]\n";
     #if($seq_length[$this_no] > 1000 && $seq_depth[$this_no] > 15 && $seq_depth[$this_no] < 100){#Attention! This may be changed when used in other studies.
     #  print "Warning: This node may should not be deleted: length(".$seq_length[$this_no].") depth(".$seq_depth[$this_no].") \n";
     #}

    splice(@seqs,$temp_offset,1);
	  splice(@seq_name,$temp_offset,1);
	  splice(@seq_length,$temp_offset,1);
	  splice(@seq_depth,$temp_offset,1);
	 
	  my @new_links=grep {$_ !~ /\t${temp_name}\t/} @links;
	  @links=@new_links;
  }

  #only length >= 350 will be used for the thread calculation
  my @threads_depths=map {$seq_depth[$_]} (grep {$seq_length[$_] >= 350} 0..$#seq_length);

  my $thread=zhongweishu(@threads_depths)*2.2;
  print "The threshold was changed into: $thread \n";
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

#rm_bubble_cp，理论上说此步骤是由于gmapper区分不开cp和mt同源ctg，而采取的临时性优化措施。但强行删除pt片段可能会导致真正来源于pt的reads被错误鉴定为mt来源，因此并不推荐使用此选项。
sub rm_bb_cp{
  my $path1=shift @_;
  my $path2=$path1;
#  foreach (@seq_cp_no){
#    my $cp_name=$seq_name[$_];

    while($path1 =~ /\{(\@*)([^\}]+)\}/g){
      my $pf=$`;
      my $pl=$';
	  my @t_nodes=split(';',$2);
	  my @non_pt;
	  foreach my $t_pt(@t_nodes){
        my $t_s=0;
	    foreach my $nn(@seq_cp_no){
          my $cp_name=$seq_name[$nn];
		  if($t_pt =~ /\b${cp_name}[\+\-]/){
            $t_s ++;
			last;
		  }		  
	    }
		push @non_pt,$t_pt if $t_s == 0;
	  }
      #my @non_pt=grep {$_ !~ /${cp_name}[\+\-]/} (split(';',$2));
      if(@non_pt >1){
        $path2=$pf.'{'.$1.join(";",@non_pt).'}'.$pl;
      }
      elsif(@non_pt == 1){
        $path2=$pf.$non_pt[0].$pl;
      }
      else{
        $path2=$path1;
      }
    }
#  }
  return $path2;
}

1;