#updated on Jul 8, 2022
#updated on Jul 29, 2022
#updated on Mar 29, 2023
#updated on Mar 19, 2024
#update to version 2.0 on Feb 13, 2025; PathCorr is avaiable now.
#updated at Feb 23, 2025: a small bug was fixed
#updated at Feb 27, 2025: add a params for applying PathCorr with all nodes; ；a bug detected: repeated link additons unsolved
#updated at Jul 24, 2025: add a params for removing NUMT-like nodes with assessement of depth
#updated at Aug 6, 2025: fix a small bug about the name of added nodes
#updated at Aug 7, 2025: add a params for simplifying the complex paths.
#updated at Sep 7, 2025: improved the simplification for bubbles. 
#updated at Sep 8, 2025: fix a bug for adding nodes.
#updated at Dec 6, 2025: optimized the strategy for simplifying the paths.
#updated at Dec 9, 2025: the re-connection of deleted paths was set to be optional.
#Recommendation params for gSim pipeline: the edge length for contig and read mapping (gMap params) should be longer than the overlap length between nodes in the graph
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
my $file3;
my $min_path_no;
my $min_ratio=2;
my $min_end_length;

my @seqs;
my @links;

my @seq_cp_no;
my @seq_mt_no;
my @seq_long_mt_name;

my @seq_name;
my @seq_length;
my @seq_depth;
my @seq_seq;
my %ori_rev=('+'=>'-','-'=>'+');
my @bandage;
my $out_pre;

my $path_num4add=20;#新增path的最小连接数
my $pathCorr;#对path进行矫正

my %seqn_uniq;
my $rm_fake=0;
my $dep_range='10:30';

sub remove_fake_node{
  my ($in,$c_all,$d_range)=@_;
  my ($d_min,$d_max)=split(/:/,$d_range);
  my %dep;
  open(IN1,$in) || die $!;
  foreach (<IN1>){
    chomp;
    my ($node,undef,$d)=split(/\t/);
    if(exists $dep{$node}){
      push @{$dep{$node}},$d;
    }
    else{
      $dep{$node}=[$d];
    }
  }
  close IN1;

  my @del_nodes;
  foreach my $nn(keys %dep){
    my $sum=0;
    my $size=@{$dep{$nn}};
    foreach my $i(@{$dep{$nn}}){
      $sum += $i;
    }
    my $avg=$sum/$size;
    if($avg < $d_min){
      push @del_nodes,$nn;
      warn "Warnning: Found a poorly supported node ($nn: depth < $d_min), which strongly indicates a mismatch of the gfa file and the long reads dataset!\n";
    }
    elsif($avg < $d_max){
      push @del_nodes,$nn;
      warn "Detect a putative contaminant node and it will be deleted: $nn \n";
    }
    else{
      ;
    }
  }
  #remove
  my $regex = join '\t|\t', map quotemeta, @del_nodes;
  $regex='\t'.$regex.'\t';
  my @filtered = grep { ! /$regex/ } @{$c_all};
  @{$c_all} = @filtered;
}

sub pathCount{#输入参数均为指针型变量
  my ($linkInfo,$pathInfo)=@_;
  my @linkcount;
  foreach (@{$linkInfo}){
    my @lk_info=split(/\t/);
    shift @lk_info;
    my ($ctg1,$ori1,$ctg2,$ori2,undef)=@lk_info;
    my $lk_code1='\b'.$ctg1.'\\'.$ori1.',\b'.$ctg2.'\\'.$ori2;
    $ori1 =~ tr/+-/-+/;
    $ori2 =~ tr/+-/-+/;
    my $lk_code2='\b'.$ctg2.'\\'.$ori2.',\b'.$ctg1.'\\'.$ori1;
    my @lk_paths=grep {/$lk_code1|$lk_code2/} @{$pathInfo};
    my $lk_num=@lk_paths;
    push @linkcount,[(@lk_info[0..3],$lk_num)];
  }
  return @linkcount;
}

sub rev_comp{
  my $seq=shift;
  my $out = reverse($seq);
  $out =~ tr/atcgATCG/tagcTAGC/;
  return $out;
}

sub find_neighbors{
  my ($in_lks,$in_sname)=@_;
  my @before_names;
  my @before_ori;
  my @after_names;
  my @after_ori;

  foreach my $this_link(@{$in_lks}){
    my (undef,$node1,$oritation1,$node2,$oritation2,undef)=split(/\t/,$this_link);
	  if ($node1 eq $in_sname && $oritation1 eq '+'){
      push @after_names,$node2;
      push @after_ori,$oritation2;
	  }
	  elsif($node2 eq $in_sname && $oritation2 eq '-'){
      push @after_names,$node1;
      push @after_ori,$ori_rev{$oritation1};
	  }
	  elsif($node1 eq $in_sname && $oritation1 eq '-'){
      push @before_names,$node2;
      push @before_ori,$ori_rev{$oritation2};
    }
    elsif($node2 eq $in_sname && $oritation2 eq '+'){
      push @before_names,$node1;
      push @before_ori,$oritation1;
	  }
    else{
      die "Fatal error detected again!\n";
	  }
  }
  return (\@before_names,\@before_ori,\@after_names,\@after_ori);
}

sub filterGraph{
  ($infile1,$infile2,$min_path_no,$min_end_length,$out_pre,$path_num4add,$pathCorr,$rm_fake,$file3,$dep_range)=@_;
   $pathCorr = 'off' if not defined  $pathCorr;

  open (infile1,$infile1) || die "Cannot open the file: $infile1\n";
  open (infile2,$infile2) || die "Cannot open the file: $infile2\n";

  my @content=<infile1>;
  chomp(@content);
  close infile1;
  
  remove_fake_node($file3,\@content,$dep_range) if $rm_fake > 0;

  @seqs=grep {/^S/} @content;
  @links=grep {/^L/} @content;

  my $olen;
  if($links[0] =~ /([0-9]+)M\s*$/){
    $olen = $1;
  }

  foreach my $each_info(@seqs){
    $each_info =~ /^S\t([^\t]+)\t([^\t]+)\t/;
    my $name=$1;
    push @seq_seq,$2;
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
    $seqn_uniq{$name}=$name;
    push @seq_length,$length;
    push @seq_depth,$temp_dep;
  }

  my $file_head=<infile2>;
  chomp(my @file2=<infile2>);
  @bandage=pathFilter::filterPath(\@file2,\@seq_name,\@seq_length,$min_end_length);

  #增加统计contig间连接数，同时删除低支持度连接
  open out_supported_links,">${out_pre}.supported_links";
  my @low_links;
  my @remain_links;
  my @now_link_count=pathCount(\@links,\@bandage);
  foreach my $lcn(0..$#now_link_count){
    print out_supported_links join("\t",@{$now_link_count[$lcn]})."\n";
    if($now_link_count[$lcn]->[-1] < $min_path_no){
      push @low_links,$now_link_count[$lcn];
    }
    else{
      #$links[$lcn] =~ /([0-9]+M)\s+$/;
      #my $this_olen=$1;
      push @remain_links,$links[$lcn];
    }  
  }
  close out_supported_links;
  @links=@remain_links;

  #增加对错误连接的处理能力
  #筛选低支持度连接，列入删除名单（上步已完成）；统计低支持度连接两侧contig端点与其他所有contig端点之间的连接数，对通过阈值的连接例如增加列表，暂时以10个N隔开。
  #1-生成待分析links,同时避免对lk-node两侧的links进行额外连接-这类短序列错误应该交给gCorr去解决
  #遍历所有可能的link组合，同时记录lk-node侧翼的links
  my @links_test;
  my @lk_links;
  foreach my $t_no(0..$#seq_name){
    my $t_node=$seq_name[$t_no];
    if($pathCorr eq 'all'){
      foreach my $i_no($t_no..$#seq_name){
        my $i_node=$seq_name[$i_no];
        push @links_test,join("\t",('L',$i_node,'+',$t_node,'+'));
        push @links_test,join("\t",('L',$i_node,'+',$t_node,'-'));
        push @links_test,join("\t",('L',$i_node,'-',$t_node,'+'));
        push @links_test,join("\t",('L',$i_node,'-',$t_node,'-'));
      }
    }
    elsif($pathCorr eq 'low'){
      foreach my $t_link(@low_links){
        push @links_test,join("\t",('L',@{$t_link}[0..1],$t_node,'+'));
        push @links_test,join("\t",('L',@{$t_link}[0..1],$t_node,'-'));
        push @links_test,join("\t",('L',$t_node,'+',@{$t_link}[2..3]));
        push @links_test,join("\t",('L',$t_node,'-',@{$t_link}[2..3]));
      }
    }
    elsif($pathCorr eq 'off'){
      ;
    }
    else{
      die "Error: Wrong params for corrPath option: ONLY all/low/off allowed.";
    }

    #记录lk-node侧翼的links
    next if ($t_node !~ /^add_node_/) || ($seq_length[$t_no] > 2 * $olen + 20);
    my @l_hits=grep {/\t${t_node}\t/} @links;
    my ($bn,$bo,$an,$ao)=find_neighbors(\@l_hits,$t_node);
    foreach my $il(0..$#{$bn}){
      foreach my $jl(0..$#{$an}){
        push @lk_links,join("\t",('L',$bn->[$il],$bo->[$il],$an->[$jl],$ao->[$jl]));
      }
    }
  }

  my @links_test_count = pathCount(\@links_test,\@bandage);
  my @links_add1 = grep {$_->[-1] >= $path_num4add} @links_test_count;
  #过滤掉已有的links
  my @links_add2;
  foreach my $test_link(@links_add1){
    my $lcode1=join("\t",('L',$test_link->[0],'\\'.$test_link->[1],$test_link->[2],'\\'.$test_link->[3]));
    my $lcode2=join("\t",('L',$test_link->[2],'\\'.$ori_rev{$test_link->[3]},$test_link->[0],'\\'.$ori_rev{$test_link->[1]}));

    my $ln1=grep {/$lcode1/} @links,@lk_links;
    my $ln2=grep {/$lcode2/} @links,@lk_links;

    push @links_add2,$test_link if $ln1 == 0 && $ln2 == 0;
  }

  #2-增加path
  my $add_no=0;
  foreach my $lk_add(@links_add2){
    #新增桥梁node
    $add_no ++;
    my $add_name='add_node_'.$add_no;
    while(exists $seqn_uniq{$add_name}){
      $add_no ++;
      $add_name='add_node_'.$add_no;
    }
    my $add_seq;
    #$links[0] =~ /([0-9]+)M\s*$/;
    #my $olen=$1;
    my $addn='N' x 20;
    my $add_len=$olen * 2 + 20;

    my @nhists1=grep {$seq_name[$_] eq $lk_add->[0]} 0..$#seq_name;
    my @nhists2=grep {$seq_name[$_] eq $lk_add->[2]} 0..$#seq_name;
    die "Error: Fatal error detected!\n" if @nhists1 != 1 || @nhists2 != 1;

    my $pnode1_no=shift @nhists1;
    my $pnode2_no=shift @nhists2;
    my $add_dep= ($seq_depth[$pnode1_no]+$seq_depth[$pnode2_no])/2;

    if($lk_add->[1] eq '+'){
      $add_seq = substr($seq_seq[$pnode1_no],0 - $olen);
    }
    elsif($lk_add->[1] eq '-'){
      $add_seq = rev_comp(substr($seq_seq[$pnode1_no],0,$olen));
    }
    else{
      die "Error: Fatal error!";
    }

    if($lk_add->[3] eq '+'){
      $add_seq .= $addn.substr($seq_seq[$pnode2_no],0,$olen);
    }
    elsif($lk_add->[3] eq '-'){
      $add_seq .= $addn.rev_comp(substr($seq_seq[$pnode2_no],0 - $olen));
    }
    else{
      die "Error: Fatal error!";
    }

    push @links,join("\t",('L',@{$lk_add}[0..1],$add_name,'+',$olen.'M'));
    push @links,join("\t",('L',$add_name,'+',@{$lk_add}[2..3],$olen.'M'));

    push @seq_name,$add_name;
    push @seq_seq,$add_seq;
    push @seq_length,$add_len;
    push @seq_depth,$add_dep;

    $seqn_uniq{$add_name}=$add_name;
    push @seqs,join("\t",('S',$add_name,$add_seq,'DP:f:'.$add_dep));
  }

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
  
      my $next_stat=0;
      if(@link_hits < 4){
        $next_stat ++;
      }
  
      my ($before_n,$before_o,$after_n,$after_o)=find_neighbors(\@link_hits,$temp_target);

      if(@{$before_n} < 2 || @{$after_n} <2){
        $next_stat ++;
      }
      
      #遍历连接对，方法为两个foreach
      my %true_link;
      my %minor_p_num;
      #my %true_link_num;
    
      foreach my $l_no(0..$#{$before_n}){
        my $left=$before_n->[$l_no];
        my $l_ori=$before_o->[$l_no];
        foreach my $r_no(0..$#{$after_n}){
          my $right=$after_n->[$r_no];
          my $r_ori=$after_o->[$r_no];
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

	        if($paths_no >= $min_path_no){
            #my $tag_null='';
            #$tag_null='disable' if $paths_no < $min_ratio * $min_path_no;
            my $tag_num=":$paths_no";
            if(exists $true_link{'left'}{$left.$l_ori}){
              push @{$true_link{'left'}{$left.$l_ori}},$right.$r_ori;
            }
            else{
              $true_link{'left'}{$left.$l_ori}=[$right.$r_ori];
            }
          
            if(exists $true_link{'right'}{$right.$r_ori}){
              push @{$true_link{'right'}{$right.$r_ori}},$left.$l_ori.$tag_num;
            }
            else{
              $true_link{'right'}{$right.$r_ori}=[$left.$l_ori.$tag_num];
            }
	        }
          else{
            $minor_p_num{'left'}{$left.$l_ori} = exists $minor_p_num{'left'}{$left.$l_ori} ? max($minor_p_num{'left'}{$left.$l_ori},$paths_no) : $paths_no;
            $minor_p_num{'right'}{$right.$r_ori} = exists $minor_p_num{'right'}{$right.$r_ori} ? max($minor_p_num{'right'}{$right.$r_ori},$paths_no) : $paths_no;
          }
	      }
      }

      next if $next_stat > 0;
      my $add_stat=0;
  
      foreach (0..$#{$before_n}){
        my $this_info=$before_n->[$_].$before_o->[$_];
        next if not exists $true_link{'left'}{$this_info};
        my $right_info=$true_link{'left'}{$this_info}->[0];
        my ($left_info,$p_num)=split(/:/,$true_link{'right'}{$right_info}->[0]);
        my $mp_num = exists $minor_p_num{'left'}{$this_info} ?  $minor_p_num{'left'}{$this_info} : -1;
        $mp_num = 0.1 if $mp_num == 0;
        $minor_p_num{'left'}{$this_info} = 0.1 if exists $minor_p_num{'left'}{$this_info} && $minor_p_num{'left'}{$this_info} == 0;
        if(($this_info eq $left_info) && $p_num / $mp_num >= $min_ratio && (@{$true_link{'left'}{$left_info}} == 1) && (@{$true_link{'right'}{$right_info}} == 1)){
          if($add_stat == 0){
            $add_stat =1;
          }
          else{
            add_copy($left_info,$this_seq,$right_info);
            print "Node has been copied: $temp_target \n";
            $add_num1=1;
          }
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