#updated in version 1.20: process the alternative nodes and remove paths that can not cover the whole read.
#updated in version 1.25: Added a new option for control the allowed max gap size in the blade region
#updated in version 1.26: fixed a small bug in the blastn command
#updated in version 1.29: try to resolve the bubble path
#updated in version 1.292: adjust the min alp to 1 when it is more than 1 in the removeFake function
#updated in version 1.295: added a new option min_iden for reads filtering;
#updated in version 1.296: added a new option min_cov_read for reads filtering;
#updated in version 1.30: added statistics about the reads located in the intra-region of contigs;
#updated in version 1.31: output alignment positions of passed reads;
#updated in version 1.40: add statistics about the 3-depth of contigs by overlaps of passed reads
#updated in version 1.50: add statistics about the align information of each passed read; fixed several bugs
#updated in version 1.51: provide stat to GraphCorrector;
#updated in version 1.511: Add more explanations to the usage discription
#updated in version 1.512: fixed the error in minIden caculation for paf format.
#updated in version 1.513: fixed the error in maxEdgeSize usage.
#updated in version 1.514: applying the filtering options for -cd function.
#updated in version 1.515: fixed the error in gap caculation for single path; replaced the maxOffset option with maxOffset1 and maxOffset2
#updated in version 1.516: support minimap2. fixed the error caused by maxOffset1 option.
#updated in version 1.517: fixed a important bug that the paf format codes misplaced the query and subject colums.
#updated in version 1.518: fixed some bugs in gap caculation.
#updated in version 1.519: fixed a bug which could cause mistake in the obtained paths.
#updated in version 1.520: fixed a bug which could cause wrong paths from the edge offset params;fixed a mistake in the values of maxoffset;update the mapping strategy of minimap2 pipeline.
#updated in version 1.521: fixed a bug when there is no contigs < 500bp
#updated in version 1.522: changed the parms of blastn
#updated in version 1.523: added a new option strictBubble to avoid worng bubbles
#updated in version 1.524: added a method to remove a ctg-read alignment that was mapped in the inter-region of another alignment
#updated in verison 1.525: fixed a big bug in the readGfa module
package graphMapper;
use strict;
#use warnings;
use List::MoreUtils qw(any uniq);
use List::Util 'max';
use List::Util 'min';
use List::Util 'sum';
use File::Basename;
use Getopt::Long;
use Getopt::Long qw(GetOptionsFromArray);
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";
use GSAT::graphIO;
#use Getopt::Long::Configure qw(bundling no_ignore_case);
#use Smart::Comments;

#更改mapping策略，minimap2仅用于>500bp的ctg
#因此重整mapping代码架构
my @align_info;
my @sub_ori;
my ($no_Qid,$no_Qlen,$no_Qs,$no_Qe,$no_Sid,$no_Slen,$no_Ss,$no_Se,$no_ori,$no_alp,$no_ai)=(0..10);
my %qlength;
my %qcoverage;
my %slength;
my $b7fmt='"7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"';
my $maxOffsetDefault=10;
my $b7_file;
my $paf_file;
my $read_file;
my $gfa_file;
my $out_prefix;
my $blastnThreads;
my $min_read_length=1000;
my $max_offset=10;
my $max_offset1;
my $max_offset2;
my $comb_max_dis=15;
my $max_edge_size=60;
my $align_stat=0;
my $max_iden_gap=1;
my $min_iden=0.85;
my $min_cov_read=0.9;
my $path_cov_min=0.9;
my $depth_stat=0;
my $only_depth=0;
my $cd_filter=0;
my $minimap2;
my $maxBounderRatio=0.1;
my $strict_bubble='yes';

sub readmap{#格式为paf/b7,\@alignfile#注意：<50bp的比对结果会被忽略。
  my ($fmt,$data)=@_;
  my $min_ali_len=50;

  my $query_id_no=0;
  my $subject_id_no=1;
  my $query_start_no=6;
  my $query_end_no=7;
  my $subject_start_no=8;
  my $subject_end_no=9;
  my $query_len_no=12;
  my $subject_len_no=13;
  my $subject_ori_no;

  my $ali_length_no=3;
  my $ali_matchs_no;
  my $ali_iden_no=2;

  if($fmt eq 'paf'){
    $query_id_no=0;
    $query_start_no=2;
    $query_end_no=3;
    $subject_id_no=5;
    $subject_start_no=7;
    $subject_end_no=8;
    $query_len_no=1;
    $subject_len_no=6;
    $subject_ori_no=4;
    $ali_length_no=10;
    $ali_matchs_no=9;
  }
  else{
    ;
  }

  foreach my $line(@{$data}){
    next if $line =~ /^#/;
    my @line_info=split(/\t/,$line);
    die"Error: The standard format 7 of the blast results is not supported!\nAlternatively the -outfmt code must be like this: -outfmt '7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'\n" if $fmt eq 'b7' && @line_info < 14;
  
    $qlength{$line_info[$query_id_no]}=$line_info[$query_len_no];
    $slength{$line_info[$subject_id_no]}=$line_info[$subject_len_no];
  
    my @selected_info=@line_info[$query_id_no,$query_len_no,$query_start_no,$query_end_no,$subject_id_no,$subject_len_no,$subject_start_no,$subject_end_no];
  
    my $temp_ori;
    my $temp_ali_len=$line_info[$ali_length_no];
    my $temp_ali_iden;
    if($fmt eq 'b7'){
      $temp_ali_iden=$line_info[$ali_iden_no];

      if($line_info[$subject_start_no] > $line_info[$subject_end_no]){
        ($selected_info[$no_Ss],$selected_info[$no_Se]) = ($selected_info[$no_Se],$selected_info[$no_Ss]);
        $temp_ori='-';
      }
      else{
        $temp_ori='+';
      }
    }
    else{
      $temp_ali_iden=$line_info[$ali_matchs_no]/$temp_ali_len*100;
      $temp_ori=$line_info[$subject_ori_no];
    }
  
    next if ($temp_ali_iden < $min_iden) || ($temp_ali_len < $min_ali_len);
    my $temp_ali_lengthprop=($line_info[$subject_end_no]-$line_info[$subject_start_no]+1)/$line_info[$subject_len_no];

    my @old_qc;
    @old_qc=@{$qcoverage{$line_info[$query_id_no]}{'dep'}} if exists $qcoverage{$line_info[$query_id_no]};
    $qcoverage{$line_info[$query_id_no]}={calc_cov(@line_info[$query_len_no,$query_start_no,$query_end_no],\@old_qc)};
  
    push @sub_ori,$temp_ori;
    push @align_info,[@selected_info,$temp_ori,$temp_ali_lengthprop,$temp_ali_iden];
  }
}

sub bounderCheck{#尚待在代码中使用,输入参数为ctg长度
  my $len=shift @_;
  my $bound=min($max_offset,$maxBounderRatio * $len);
  return $bound;
}

sub GetEdgeStat{
  my $the_dis=shift @_;
  #my $max_ln_length=shift @_;
  my $mid_point=shift @_;
  my $the_offset=shift @_;
  my $edge_code;
  my $code_prefix='';

  if($the_dis >= $mid_point + $the_offset){
    if($the_dis > 0){
        $edge_code = 9;#too far;Gap detected
    }
    else{
        $edge_code = 8;#too far; But no gap was observed.
    } 
    $code_prefix="($the_dis)";#right far than the mid-point
  }
  elsif($the_dis < $mid_point - $the_offset){
    $edge_code = -9;#too close
    $code_prefix="[$the_dis]";#left far than the mid-point
  }
  else{
    $edge_code = 1;#good
  }
  return ($edge_code,$code_prefix);
}

sub calc_cov{
  my ($seq_len,$ali_s,$ali_e,$old)=@_;
  my @seq_dep;
  if(@{$old} != $seq_len){
    if(@{$old} == 0){
      @seq_dep=(0) x $seq_len;
    }
    else{
      die "Error: Different lengths found for the same sequence: ".$#{$old}." / $seq_len !\n";
    }
  }
  else{
    @seq_dep=@{$old};
  }

  foreach ($ali_s-1 .. $ali_e-1){
    $seq_dep[$_]++;
  }

  my $ali_num=grep {$_ > 0} @seq_dep;
  my $ali_prop=$ali_num/$seq_len;
  my %seq_stat=('cov'=>$ali_prop,'dep'=>[@seq_dep]);

  return %seq_stat;
}

sub phaseOpt{
  GetOptionsFromArray($_[0],
        'blast7File=s'       => \$b7_file,
        'pafFile=s'          => \$paf_file,
        'readFile=s'         => \$read_file,
        'gfaFile=s'          => \$gfa_file,
        'out=s'              => \$out_prefix,
        'minRead=i'          => \$min_read_length,
        'maxOffset1=i'       => \$max_offset1,
        'maxOffset2=i'       => \$max_offset2,
        'maxCombDis=i'       => \$comb_max_dis,
        'maxEdgeSize1=i'     => \$max_edge_size,
        'maxEdgeSize2=i'     => \$max_offset,
        'maxIdenGap=f'       => \$max_iden_gap,
        'minIden=f'          => \$min_iden,
        'minCovofRead=f'     => \$min_cov_read,
        'minCovbyPath=f'     => \$path_cov_min,
        'strictBub=s'        => \$strict_bubble,
        'align!'             => \$align_stat,
        'depth!'             => \$depth_stat,
        'calDepth|cd!'       => \$only_depth,
        'filterPaths!'       => \$cd_filter,
        'minimap2=s'         => \$minimap2,
        'maxBounderRatio=f'  => \$maxBounderRatio,
        );
  my $option_lows=1;
  my $err_code;

  if(! defined($gfa_file)){
    $option_lows = 0;
    $err_code='Fatal: No gfa file!';
  }
  elsif($align_stat > 0 && (! defined($read_file))){
    $option_lows = 0;
    $err_code='Fatal: No read file specified!';
  }
  elsif($align_stat > 0 && (defined($b7_file) || defined($paf_file))){
    $option_lows = 0;
    $err_code='Fatal: More than 1 alignment file specified!';
  }
  elsif($align_stat == 0 && $only_depth == 0 && (! defined($b7_file)) && (! defined($paf_file))){
    $option_lows = 0;
    $err_code='Fatal: No alignment file specified!';
  }
  elsif(defined($b7_file) && defined($paf_file)){
    $option_lows = 0;
    $err_code='Fatal: More than 1 alignment file specified!';
  }
  elsif($only_depth > 0 && (! defined($out_prefix))){
    $option_lows = 0;
    $err_code='Fatal: -out option is required!';
  }
  elsif($cd_filter > 0 && (defined($min_iden) || defined($min_cov_read))){
    $err_code='';
    warn 'Warning: --minIden ['.$min_iden.'] and --minCovofRead ['.$min_cov_read.'] options will be ignored when -f option is applied! Only --minRead and --minCovbyPath options can be applied in this function now.\n';
  }
  elsif(defined($max_offset1) && defined($max_offset2)){
    $option_lows = 0;
    $err_code='Fatal: --maxOffset2 is not allowed when --maxOffset1 is specified!';
  }
  else{
    ;
  }
  if (! $option_lows){
    print "\n  $err_code \n";
    exit(0);
  }
}

my $gfa_S;
my $gfa_L;
sub gmap{
  phaseOpt(@_);
  $max_offset1=$maxOffsetDefault if !(defined($max_offset1) || defined($max_offset2));
  $out_prefix = basename($gfa_file) if ! defined($out_prefix);

    open(gfaFile,"$gfa_file") || die "Error: Cannot open the file: $gfa_file \n";
  chomp(my @gfa_content=<gfaFile>);
  
  my ($gfa_S,undef,$gfa_L)=graphIO::readGfa(\@gfa_content,'SL');
  my $gfa_lnk_k=${$gfa_L}[0]->[-1];
  
  my $offset_center;
  my $offset_length;
  my $ctg_start;
  if(defined($max_offset1)){
    $offset_center=1 - $gfa_lnk_k;
    $offset_length=$max_offset1;
    $ctg_start=$gfa_lnk_k;
  }
  else{
    $offset_center=0;
    $offset_length=$max_offset2;
    $ctg_start=1;
  }
  
  goto CalDepth if ($only_depth > 0);
  
  my $long_num=0;
  my $short_num=0;
  if($align_stat > 0){
    #open(read_file,"$read_file") || die "Error: Cannot open the file: $read_file \n";
    open(gfaFas1,">${out_prefix}.part1.fas");
    open(gfaFas2,">${out_prefix}.part2.fas");
  
    foreach my $each(sort keys %{$gfa_S}){
      if(${$gfa_S}{$each}{'len'} >= 500){
        print gfaFas1 '>'.$each."\n".${$gfa_S}{$each}{'seq'}."\n";
        $long_num ++;
      }
      else{
        print gfaFas2 '>'.$each."\n".${$gfa_S}{$each}{'seq'}."\n";
        $short_num ++;
      }  
    }
    close gfaFas1;
    close gfaFas2;
  
    if(defined($minimap2) && $long_num > 0){
      my $r_fmt;
      if($minimap2 eq 'hifi'){
        $r_fmt='map-hifi';
      }
      elsif($minimap2 eq 'clr'){
        $r_fmt='map-pb';
      }
      elsif($minimap2 eq 'ont'){
        $r_fmt='map-ont';
      }
      else{
        die"Error: invalid values for --minimap2 option!\n";
      }
      
      print "[info] align long ctgs to the long reads with minimap2.\n";
      `minimap2 -x $r_fmt --cs -t 8 -o ${out_prefix}.part1.reads.paf ${out_prefix}.part1.fas $read_file`;
  
      $paf_file="${out_prefix}.part1.reads.paf";
      open (align_file,$paf_file) || die "Error: Cannot open the file: $paf_file ! \n";
      chomp(my @mapped=<align_file>);
      readmap('paf',\@mapped);
      splice @mapped;
      close align_file;
  
      if($short_num > 0){
        print "[info] align short ctgs to the long reads with blastn.\n";
        `blastn -query $read_file -subject ${out_prefix}.part2.fas -task blastn -dust no -outfmt $b7fmt -out ${out_prefix}.part2.reads.b7`;#注意这里默认返回的最大hits数量是500条
        $b7_file="${out_prefix}.part2.reads.b7";
        open (align_file,$b7_file) || die "Error: Cannot open the file: $b7_file ! \n";
        chomp(my @mapped=<align_file>);
        readmap('b7',\@mapped);
        splice @mapped;
        close align_file;
      }
      else{
        warn "Warning: No contigs < 500 bp were found!\n";
      }
    }
    elsif($short_num + $long_num > 0){
      print "[info] align all ctgs to the long reads with blastn.\n";
      `cat ${out_prefix}.part1.fas ${out_prefix}.part2.fas > ${out_prefix}.fas`;
      `blastn -query $read_file -subject ${out_prefix}.fas -dust no -outfmt $b7fmt -out ${out_prefix}.reads.b7`;#注意这里默认返回的最大hits数量是500条
  
      $b7_file="${out_prefix}.reads.b7";
      open (align_file,$b7_file) || die "Error: Cannot open the file: $b7_file ! \n";
      chomp(my @mapped=<align_file>);
      readmap('b7',\@mapped);
      splice @mapped;
      close align_file;
    }
    else{
      die "Error: No contigs found in the gfa file!\n";
    }
  }
  else{
    if(defined($paf_file)){
      open (align_file,$paf_file) || die "Error: Cannot open the file: $paf_file ! \n";
    }
    elsif(defined($b7_file)){
      open (align_file,$b7_file) || die "Error: Cannot open the file: $b7_file ! \n";
    }
    else{
      die"Error: No valid .paf or .b7 file was found!\n";
    }
    chomp(my @mapped=<align_file>);
    readmap('b7',\@mapped);
    splice @mapped;
    close align_file;
  }
  
  open out_path,">${out_prefix}.mapping.paths";
  print out_path "Read_id\tRead_length\tGapped_length\tMapped_ratio\tMapped_path\tFirstCtgAliPos\tLastCtgAliPos\tMappedCtgs\tAdjustedReadAliPos\tAdjustedCtgAlipos\n";
  open out_reads_stat,">${out_prefix}.mt.reads";
  print out_reads_stat "Read_id\tRead_length\tRead_mt_coverage\tMappedEdgesOfContigs\tMappedContigs\n";
  
  foreach my $this_read(keys %qlength){
    next if $qlength{$this_read} < $min_read_length;
    next if $qcoverage{$this_read}{'cov'} < $min_cov_read;
    #print out_reads_stat $this_read."\t".$qlength{$this_read}."\t".$qcoverage{$this_read}{'cov'}."\n";
  
    my @sel_rows=grep {$align_info[$_]->[$no_Qid] eq $this_read} 0..$#align_info;
    my @these_ctgs=uniq (map {$align_info[$_]->[$no_Sid]} @sel_rows);
  
    my @this_merged;
    #my @this_merged_ori;
    foreach my $this_ctg(@these_ctgs){
      my @temp_m_f;
      my @temp_m_r;
      #正向
      #my @temp_no_stat;
      my @temp_no_f=sort {$align_info[$a]->[$no_Ss] <=> $align_info[$b]->[$no_Ss]} (grep {$align_info[$_]->[$no_Sid] eq $this_ctg && $sub_ori[$_] eq '+'} @sel_rows);
      if(@temp_no_f >=2){
        foreach my $i(0..$#temp_no_f){
          my @temp_merge=@{$align_info[$temp_no_f[$i]]};
          foreach my $j($i+1..$#temp_no_f){
            my $dis_len1 = $align_info[$temp_no_f[$j]]->[$no_Ss] - $temp_merge[$no_Se];#ctg比对，可正可负
            my $dis_len2 = $align_info[$temp_no_f[$j]]->[$no_Qs] - $temp_merge[$no_Qe];#read比对，可正可负
            my $hit_len1 = $align_info[$temp_no_f[$j]]->[$no_Se] - $align_info[$temp_no_f[$j]]->[$no_Ss]+1;
            my $hit_len2 = $temp_merge[$no_Se] - $temp_merge[$no_Ss]+1;
            my $min_ctg_len=$slength{$this_ctg}-bounderCheck($slength{$this_ctg});
  
            if((!($hit_len1 >= $min_ctg_len && $hit_len2 >= $min_ctg_len)) && abs($dis_len1) <= $comb_max_dis && abs($dis_len2) <= $comb_max_dis){
              $temp_merge[$no_Qe]=max($align_info[$temp_no_f[$j]]->[$no_Qe],$temp_merge[$no_Qe]);
              $temp_merge[$no_Se]=max($align_info[$temp_no_f[$j]]->[$no_Se],$temp_merge[$no_Se]);
              #合并会导致align_info的-2, -1两个元素值，即比对长度比例和比对准确度，不再准确，因此也需要进行重新计算#$no_alp,$no_ai
              $temp_merge[$no_alp]=($temp_merge[$no_Se]-$temp_merge[$no_Ss]+1)/$slength{$this_ctg};
              #合并后的iden难以计算，暂时取近似值，谨慎使用合并相关的参数
              $temp_merge[$no_ai]=$hit_len1/($hit_len1+$hit_len2)*$align_info[$temp_no_f[$j]]->[$no_ai] + $hit_len2/($hit_len1+$hit_len2)*$temp_merge[$no_ai];
            }
            else{
              ;
            }
          }
          push @temp_m_f,[@temp_merge];
        }
      }
      elsif(@temp_no_f ==1){
        push @temp_m_f,[@{$align_info[$temp_no_f[0]]}];
      }
      else{
        ;
      }
  
      #去重 / 去碎片，其他
      my @the_m_f;
      foreach my $each (@temp_m_f){
        my $retain=0;
        my $ctg_edge_max=bounderCheck($slength{$each->[$no_Sid]});
        if($each->[$no_Qs] <= $max_edge_size && $each->[$no_Se] >= ($slength{$each->[$no_Sid]} - $ctg_edge_max)){
          $retain=1;
        }
        elsif($each->[$no_Ss] <= $ctg_edge_max && $each->[$no_Qe] >= ($qlength{$this_read} - $max_edge_size)){
          $retain=1;
        }
        elsif(($each->[$no_Ss] <= $ctg_edge_max && $each->[$no_Se] >= ($slength{$each->[$no_Sid]} - $ctg_edge_max)) || ($each->[$no_Qs] <= $max_edge_size && $each->[$no_Qe] >= ($qlength{$each->[$no_Qid]} - $max_edge_size))){#无论ctg或者read全长比上，均保留
          $retain=2;
        }
        else{
          next;
        }
  
        $retain=2 if $retain == 1 && ($each->[$no_Qe] - $each->[$no_Qs]) > ($gfa_lnk_k + $ctg_edge_max);
  
        push @the_m_f,[@{$each}] if $retain == 2;
      }
  
      #反向
      my @temp_no_r=sort {$align_info[$a]->[$no_Ss] <=> $align_info[$b]->[$no_Ss]} (grep {$align_info[$_]->[$no_Sid] eq $this_ctg && $sub_ori[$_] eq '-'} @sel_rows);
      if(@temp_no_r >=2){
        foreach my $i(0..$#temp_no_r){
          my @temp_merge=@{$align_info[$temp_no_r[$i]]};
          foreach my $j($i+1..$#temp_no_f){
            my $dis_len1 = $align_info[$temp_no_r[$j]]->[$no_Ss] - $temp_merge[$no_Se];
            my $dis_len2 = $temp_merge[$no_Qs] - $align_info[$temp_no_r[$j]]->[$no_Qe];#q和s的se是反向对应的#updated in v1.520
            my $hit_len1 = $align_info[$temp_no_r[$j]]->[$no_Se] - $align_info[$temp_no_r[$j]]->[$no_Ss];
            my $hit_len2 = $temp_merge[$no_Se] - $temp_merge[$no_Ss];
            my $min_ctg_len=$slength{$this_ctg}-bounderCheck($slength{$this_ctg});
  
            if((!($hit_len1 >= $min_ctg_len && $hit_len2 >= $min_ctg_len))&& abs($dis_len1) <= $comb_max_dis && abs($dis_len2) <= $comb_max_dis){
              $temp_merge[$no_Qs]=min($align_info[$temp_no_r[$j]]->[$no_Qs],$temp_merge[$no_Qs]);#q和s的se是反向对应的
              $temp_merge[$no_Se]=max($align_info[$temp_no_r[$j]]->[$no_Se],$temp_merge[$no_Se]);
              #合并会导致align_info的-2, -1两个元素值，即比对长度比例和比对准确度，不再准确，因此也需要进行重新计算#$no_alp,$no_ai
              $temp_merge[$no_alp]=($temp_merge[$no_Se]-$temp_merge[$no_Ss]+1)/$slength{$this_ctg};
              #合并后的iden难以计算，暂时取近似值，谨慎使用合并相关的参数
              $temp_merge[$no_ai]=$hit_len1/($hit_len1+$hit_len2)*$align_info[$temp_no_f[$j]]->[$no_ai] + $hit_len2/($hit_len1+$hit_len2)*$temp_merge[$no_ai];
            }
          }
          push @temp_m_r,[@temp_merge];
        }
      }
      elsif(@temp_no_r ==1){
        push @temp_m_r,[@{$align_info[$temp_no_r[0]]}];
      }
      else{
        ;
      }
  
      #去重 / 去碎片，其他
      my @the_m_r;
      foreach my $each (@temp_m_r){
        my $retain=0;
        my $ctg_edge_max=bounderCheck($slength{$each->[$no_Sid]});
        if($each->[$no_Qs] <= $max_edge_size && $each->[$no_Ss] <= $ctg_edge_max){
          $retain=1;
        }
        elsif($each->[$no_Se] >= ($slength{$each->[$no_Sid]} - $ctg_edge_max) && $each->[$no_Qe] >= ($qlength{$this_read} - $max_edge_size)){
          $retain=1;
        }
        elsif(($each->[$no_Ss] <= $ctg_edge_max && $each->[$no_Se] >= ($slength{$each->[$no_Sid]} - $ctg_edge_max)) || ($each->[$no_Qs] <= $max_edge_size && $each->[$no_Qe] >= ($qlength{$each->[$no_Qid]} - $max_edge_size))){#无论ctg或者read全长比上，均保留
          $retain=2;
        }
        else{
          next;
        }
  
        $retain=2 if $retain == 1 && ($each->[$no_Qe] - $each->[$no_Qs]) > ($gfa_lnk_k + $ctg_edge_max);
  
        push @the_m_r,[@{$each}] if $retain == 2;
      }
  
      #整理结果赋值，并去除可能存在的空元素
      foreach (@the_m_f){
        if(@{$_} > 0){
          push @this_merged,[@{$_}];
          #push @this_merged,[@{$_},'+'];
          #push @this_merged_ori,'+';
        }
      }
      foreach (@the_m_r){
        if(@{$_} > 0){
          push @this_merged,[@{$_}];
          #push @this_merged,[@{$_},'-'];
          #push @this_merged_ori,'-';
        }
      }
    }
  
    #排列，编码，输出path
    my @sorted_merged=sort {$a->[$no_Qs] <=> $b->[$no_Qs]} @this_merged;
    my @ctg_with_ends=uniq (map {$_->[$no_Sid]} @sorted_merged);
    my $edge_ali=@ctg_with_ends;#是否比对到末端
     
    print out_reads_stat $this_read."\t".$qlength{$this_read}."\t".$qcoverage{$this_read}{'cov'}."\t".$edge_ali."\t".join(';',@ctg_with_ends)."\n";
    next if @ctg_with_ends <1;
  
    #包含性鉴定与过滤
    my @sorted_merge_final;
    my ($test_s,$test_e)=($sorted_merged[0]->[$no_Qs],$sorted_merged[0]->[$no_Qe]);
    if(@sorted_merged >= 2){
      #my @rm_nos;
      foreach my $test_info(@sorted_merged){
        if(($test_info->[$no_Qs] > $test_s) && ($test_info->[$no_Qe] < $test_e)){
          #push @rm_nos,$test_no;
        }
        else{
          push @sorted_merge_final,$test_info;
  	      ($test_s,$test_e)=($test_info->[$no_Qs],$test_info->[$no_Qe]);
        }
      }
      @sorted_merged=@sorted_merge_final;
    }  
  
    my $ali_head_ctg_pos;
    my $ali_tail_ctg_pos;
    #处理多个contigs匹配于同一位置的问题
    #通过grep，按Qs对contigs进行统一处理；如果存在多个同样Qs，对其与后续的contig进行连接分析-两个都通过，则均保留且输出其比对信息-仅有1个通过，则去掉其他-都不通过，则（！）
    #设计原则为：同一个contig放同一行；path通过符号体系展示bubble结构；最好能输出bubble的具体标识，{}作为标志，便于后续与专门处理的代码对接
    #对ctg排序规则进行重新调整：末端需符合两个标准之一（无须额外判断），内部需符合内含标准（判断，如果不符合，则直接排除掉）。
    my $now_no=0;
    my @now_the_nos;
    my @path_codes; #($sorted_merged[0]->[$no_Sid].$sorted_merged[0]->[$no_ori]);
    my $gapped_length=0;
    my $read_len=$qlength{$this_read};
    #my $max_lnk_length=$gfa_lnk_k + $max_offset;
    my $target_last_Qe = $ctg_start;
    my @passed_merge_no;
    my $align_r_mins=min(map {$_->[$no_Qs]} @sorted_merged);
    my $align_r_maxs=max(map {$_->[$no_Qs]} @sorted_merged);
    #my $align_r_maxe=max(map {$_->[$no_Qe]} @sorted_merged);
  
    while($now_no <= $#sorted_merged){
      my $this_no=$now_no;
      $now_no ++;
      my $f_dis;
      my $f_stat;
      my $f_code_p;
      #匹配在read内部或者path内部，都必须满足内含标准
      if((($sorted_merged[$this_no]->[$no_Qs] > $max_edge_size) || ($sorted_merged[$this_no]->[$no_Qs] > $align_r_mins)) && (($sorted_merged[$this_no]->[$no_Qe] < ($read_len - $max_edge_size)) || ($sorted_merged[$this_no]->[$no_Qs] < $align_r_maxs))){
        if(($sorted_merged[$this_no]->[$no_Ss] > bounderCheck($slength{$sorted_merged[$this_no]->[$no_Sid]})) || ($sorted_merged[$this_no]->[$no_Se] < $slength{$sorted_merged[$this_no]->[$no_Sid]} - bounderCheck($slength{$sorted_merged[$this_no]->[$no_Sid]}))){
        	print "warning: this ctg [ ".$sorted_merged[$this_no]->[$no_Sid]." ] with be skipped for the read: $this_read \n";
  	      next;
        }
      }
  
      push @now_the_nos,$this_no;
  
      #控制索引值不能超过边界，否则每次引用都会自动创建新元素
      if($now_no <= $#sorted_merged){
        next if $sorted_merged[$now_the_nos[0]]->[$no_Qs] == $sorted_merged[$now_no]->[$no_Qs];
      }
  
      print "Now is the read: $this_read ; and ctg_nos: @now_the_nos \n";#only for test
  
      if(@now_the_nos == 1){#如果只有一个，开头结尾也要判断gap大小。此处仅判断了开头。【记上末尾的gap】
        my $the_no=$now_the_nos[0];
        $ali_head_ctg_pos=$sorted_merged[$the_no]->[$no_Ss].'-'.$sorted_merged[$the_no]->[$no_Se] if ! defined($ali_head_ctg_pos);
        $ali_tail_ctg_pos=$sorted_merged[$the_no]->[$no_Ss].'-'.$sorted_merged[$the_no]->[$no_Se];
         
        my $the_dis=$sorted_merged[$the_no]->[$no_Qs] - $target_last_Qe;        
        my ($ln_stat,$code_pre)=GetEdgeStat($the_dis,$offset_center,$offset_length);
        
        #my $end_dis;
        #if($now_no > $#sorted_merged){
        #  $end_dis=$read_len - $sorted_merged[$the_no]->[$no_Qe];#单路径-处理末端gap
        #}
        #my ($end_ln_stat,$end_code)=GetEdgeStat($end_dis,$offset_center,$offset_length);
  
        push @path_codes,$code_pre.$sorted_merged[$the_no]->[$no_Sid].$sorted_merged[$the_no]->[$no_ori];
        push @passed_merge_no,$the_no;
        $gapped_length += $the_dis if $ln_stat == 9;
        #$gapped_length += $end_dis if $end_ln_stat == 9;
  
        $target_last_Qe=$sorted_merged[$the_no]->[$no_Qe];
      }
      else{
        if($target_last_Qe == 0){#仅保留作为备用，实际上不起作用。
          $f_stat=1;
          $f_code_p='';
        }
        else{
          #前面
          $f_dis=$sorted_merged[$now_the_nos[0]]->[$no_Qs] - $target_last_Qe;
          ($f_stat,$f_code_p)=GetEdgeStat($f_dis,$offset_center,$offset_length);
          $gapped_length += $f_dis if $f_stat == 9;
        }
  
        #考虑通过比较bubble结构的相似性，block百分比和iden都要最大，才能排除
        if($max_iden_gap < 1){
          @now_the_nos=removeFake(\@sorted_merged,@now_the_nos);#暂不确定这一步该放在哪个环节。
          die"all nos were removed!\n" if @now_the_nos == 0;#test
        }
  
        my @record_codes=map {$sorted_merged[$_]->[$no_Sid].$sorted_merged[$_]->[$no_ori]} @now_the_nos;
        my $record_code='{'.$f_code_p.join(';',@record_codes).'}';
        my $record_pos=join(';',map {$sorted_merged[$_]->[$no_Ss].'-'.$sorted_merged[$_]->[$no_Se]} @now_the_nos);
        my $rec_last_Qe=max(map {$sorted_merged[$_]->[$no_Qe]} @now_the_nos);
        my $last_Qe=$rec_last_Qe;
  
        #后面
  #      if($now_no <= $#sorted_merged){
          my @the_nos_stat;
          my $temp_dis;
          my $temp_stat;
          foreach my $the_n(@now_the_nos){
            if($now_no <= $#sorted_merged){
              $temp_dis=$sorted_merged[$now_no]->[$no_Qs] - $sorted_merged[$the_n]->[$no_Qe];
              ($temp_stat,undef)=GetEdgeStat($temp_dis,$offset_center,$offset_length);
            }
            else{
              $temp_dis=$read_len - $sorted_merged[$the_n]->[$no_Qe];#多路径-处理末端gap
              ($temp_stat,undef)=GetEdgeStat($temp_dis,0,$offset_length);
            }
            #($temp_stat,undef)=GetEdgeStat($temp_dis,$offset_center,$offset_length);
            push @the_nos_stat,$temp_stat;
          }
  
          my @passed_nos=map {$now_the_nos[$_]} (grep {$the_nos_stat[$_] == 1} 0..$#now_the_nos);
  
          if($strict_bubble =~ /^yes$|^Yes$|^Y$|^y$/){
            my $max_end=max(map {$sorted_merged[$_]->[$no_Qe]} @passed_nos);
            @passed_nos=grep {$sorted_merged[$_]->[$no_Qe] == $max_end} @passed_nos;
          }
  
          if(@passed_nos >=1){
            my @the_codes=map {$sorted_merged[$_]->[$no_Sid].$sorted_merged[$_]->[$no_ori]} @passed_nos;
            my $the_code=$f_code_p.join(';',@the_codes);
            $the_code = "{$the_code}" if @passed_nos >=2;
            push @path_codes,$the_code;
            push @passed_merge_no,$passed_nos[0] if @passed_nos == 1; #存在备选形式的ctg将不会被记录为可信结果。
            #pos_end
            $last_Qe=max(map {$sorted_merged[$_]->[$no_Qe]} @passed_nos);
            $ali_head_ctg_pos=join(';',map {$sorted_merged[$_]->[$no_Ss].'-'.$sorted_merged[$_]->[$no_Se]} @passed_nos) if ! defined($ali_head_ctg_pos);
            $ali_tail_ctg_pos=join(';',map {$sorted_merged[$_]->[$no_Ss].'-'.$sorted_merged[$_]->[$no_Se]} @passed_nos);
          }
          else{          
            push @path_codes,'{@'.$f_code_p.join(';',@record_codes).'}';#contig连接关系无法确认最佳路径
          }
  #      }
  #      else{
  #        push @path_codes,$record_code;
  #        push @passed_nos,@now_the_nos;#注意：对于buble结构的矫正存在风险！因此，建议只采信无bubble结构的矫正结果。
  #        $ali_head_ctg_pos=$record_pos if ! defined($ali_head_ctg_pos);
  #        $ali_tail_ctg_pos=$record_pos;
  #        $last_Qe=
  #      }
  
        $target_last_Qe = $last_Qe;
      }
  
      splice @now_the_nos;
    }
  
    my $end_dis = $read_len - $target_last_Qe;
    my ($end_stat,$end_code)=GetEdgeStat($end_dis,0,$offset_length);
    if($end_stat == 9){
      $gapped_length += $end_dis;
      $path_codes[-1] .= $end_code;
    }
  
    my $cov_ratio=($read_len - $gapped_length)/$read_len;
    my @pass_Sid;
    my @pass_ori;
    my @pass_Qse;
    my @pass_Sse;
    foreach my $pass_no(@passed_merge_no){
      push @pass_Sid,$sorted_merged[$pass_no]->[$no_Sid];
      push @pass_ori,$sorted_merged[$pass_no]->[$no_ori];
      push @pass_Qse,$sorted_merged[$pass_no]->[$no_Qs].','.$sorted_merged[$pass_no]->[$no_Qe];
      push @pass_Sse,$sorted_merged[$pass_no]->[$no_Ss].','.$sorted_merged[$pass_no]->[$no_Se];
    }
    my $out_ali=join(';',@pass_Sid)."\t".join(';',@pass_ori)."\t".join(';',@pass_Qse)."\t".join(';',@pass_Sse);
    #@passed_align_info = (@passed_align_info,(map {$sorted_merged[$_]} @passed_merge_no));
    print out_path $this_read."\t".$read_len."\t".$gapped_length."\t".$cov_ratio."\t".join(',',@path_codes)."\t${ali_head_ctg_pos}\t${ali_tail_ctg_pos}\t${out_ali}\n" if @path_codes > 0 && $cov_ratio >= $path_cov_min;
    
  }
  close out_path;
  close out_reads_stat;
  
  CalDepth:
  if($depth_stat >0 || $only_depth > 0){
    #my $intra_c=4;
    #my $intra_p=5;
    my $r_len=1;
    my $p_co=3;
    my $path_p=4;
    my $path_hs=5;
    my $path_te=6;
    #open (mt_read,"${out_prefix}.mt.reads") || die "Error: Cannot open the ${out_prefix}.mt.reads \n";
    open (mt_paths,"${out_prefix}.mapping.paths") || die "Error: Cannot open the ${out_prefix}.mapping.paths \n";
  
    #my @intra_cp=map {chomp($_);[(split(/\t/))[$intra_c,$intra_p]]} <mt_read>;
    #@intra_cp=grep {$_->[0] =~ /\(.+\)/} @intra_cp;
  
    #foreach (@intra_cp){
    #  $_->[0] =~ s/[0-9]+\(|\)//g;
    #}
    my @path_cp;
    my $useless=<mt_paths>;
    foreach my $line(<mt_paths>){
      chomp($line);
      my ($r_length,$p_cov,$this_path,$this_s,$this_e)=(split(/\t/,$line))[$r_len,$p_co,$path_p,$path_hs,$path_te];
      
      #应用过滤参数
      if($cd_filter > 0 && ($r_length < $min_read_length || $p_cov < $path_cov_min)){
        next;
      }
      
      my @path_allc=grep {$_ !~ /^\{@/} (split(/,/,$this_path));
      my @this_ee=split(/;/,$this_e);
      my @this_ss=split(/;/,$this_s);
      my $ctg_no=-1;
  
      foreach my $c_no(0..$#path_allc){
        my $temp_c=$path_allc[$c_no];
        $temp_c =~ s/\([^\)]+\)//g;
        $temp_c =~ s/\[[^\]]+\]//g;
        $temp_c =~ s/[\{\}\+\-]//g;
        next if length($temp_c) <= 0;
        #print $temp_c;die;
        $ctg_no ++;
        my @this_cc=split(/;/,$temp_c);
  
        foreach my $cc(0..$#this_cc){
          if($ctg_no == 0){
            push @path_cp,[$this_cc[$cc],$this_ss[$cc]];
          }
          elsif($c_no == $#path_allc){
            push @path_cp,[$this_cc[$cc],$this_ee[$cc]];
          }
          else{
            push @path_cp,[$this_cc[$cc],'1-'.${$gfa_S}{$this_cc[$cc]}{'len'}];
          }
        }     
      }
    }
    #close mt_read;
    close mt_paths;
    #print join(";",map {$_->[0]} @intra_cp[0..4])."\n";#test
    #my %mt_depth=calcu_depth(@intra_cp,@path_cp);
    my %mt_depth=calcu_depth(@path_cp);
  
    open out_depth,">${out_prefix}.passed.depth";
    foreach my $each_ctg(sort keys %mt_depth){
      foreach my $each_pos(0..$#{$mt_depth{$each_ctg}}){
        my $each_base=$each_pos+1;
        print out_depth "${each_ctg}\t${each_base}\t".$mt_depth{$each_ctg}->[$each_pos]."\n";
      }
    }
    close out_depth;
  }
}

sub removeFake{#移除条件：匹配长度百分比必须最低；可选项，iden小于最优值的差异。原则：只要与集合中最高iden的差异超过阈值即可。注意：该功能会导致bubble结构与uni-contig区域reads的筛选标准不一致，影响深度统计结果，因此不建议使用。
  my $the_merged=shift;
  my @the_nos=@_;
  my $max_iden=max(map {${$the_merged}[$_]->[$no_ai]} @the_nos);
  my $min_alp=min(map {${$the_merged}[$_]->[$no_alp]} @the_nos);
  my @sel_nos=grep {my $adj_this_alp=${$the_merged}[$_]->[$no_alp];$adj_this_alp = 1 if $adj_this_alp > 1;!(($max_iden - ${$the_merged}[$_]->[$no_ai] >= $max_iden_gap) && ($adj_this_alp <= $min_alp))} @the_nos;

  return @sel_nos;
}

sub calcu_depth{
  my %result_dep;
  foreach my $this_rec(@_){
    my ($ctg_n,$ctg_ali)=@{$this_rec};######why empty  
    #print join("-",@{$this_rec})."\n";#test
    $result_dep{$ctg_n}=[map {'0'} (0 .. (${$gfa_S}{$ctg_n}{'len'} - 1))] if not exists $result_dep{$ctg_n};
    #print $ctg_ali."\n";#test
    my ($ctg_a_s,$ctg_a_e)=split(/\-/,$ctg_ali);
    foreach ($ctg_a_s-1 .. $ctg_a_e-1){
      $result_dep{$ctg_n}->[$_] ++;
    }
  }
  return %result_dep;  
}

1;