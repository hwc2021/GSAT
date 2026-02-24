#generated at Feb 25, 2025
#updated at Feb 28, 2025
#updated at Mar 5, 2025
#updated at Aug 18, 2025
#updated at Aug 23, 2025
#updated at Nov 30, 2025
#updated at Dec 3, 2025

package graphClip;

use strict;
use warnings;
use List::MoreUtils qw(any uniq);
use List::Util 'max';
use List::Util 'min';
use FindBin;
use File::Basename;
use Getopt::Long qw(GetOptionsFromArray);
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";
use GSAT::graphIO;

my @sub_ori;
my @align_info;
my ($no_Qid,$no_Qlen,$no_Qs,$no_Qe,$no_Sid,$no_Slen,$no_Ss,$no_Se,$no_ori,$no_alp,$no_ai)=(0..10);
my %qlength;
my %slength;
my $min_iden = 0.8;
my %qcoverage;

sub readmap2{#格式为paf/b7,\@alignfile#注意：<50bp的比对结果会被忽略。
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
    elsif($fmt eq 'paf'){
      $temp_ali_iden=$line_info[$ali_matchs_no]/$temp_ali_len*100;
      $temp_ori=$line_info[$subject_ori_no];
      $selected_info[$no_Ss] += 1;#covert to 1-based position
      $selected_info[$no_Qs] += 1;#covert to 1-based position
    }
    else{
      die "Error: Wrong alignment file type!\n";
    }
  
    next if ($temp_ali_iden < 100*$min_iden) || ($temp_ali_len < $min_ali_len);
    my $temp_ali_lengthprop=($line_info[$subject_end_no]-$line_info[$subject_start_no]+1)/$line_info[$subject_len_no];

    my @old_qc;
    @old_qc=@{$qcoverage{$line_info[$query_id_no]}{'dep'}} if exists $qcoverage{$line_info[$query_id_no]};
    $qcoverage{$line_info[$query_id_no]}={calc_cov(@line_info[$query_len_no,$query_start_no,$query_end_no],\@old_qc)};
  
    push @sub_ori,$temp_ori;
    push @align_info,[@selected_info,$temp_ori,$temp_ali_lengthprop,$temp_ali_iden];
  }
  #return(\@sub_ori,@align_info);
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

my $mt_read_info;
my $gfa_file;
my $b7_file;
my $paf_file;
my $min_path_num;
my $min_len=100;
my $out_prefix;
#my $offset_len = 1;#step_size
#my $obs_win=200;
my $top_hits=100;
my $top_min=20;

my $gfa_S;
my $gfa_L;

sub phaseOpt{
  GetOptionsFromArray($_[0],
        'mtReadList=s'     => \$mt_read_info,
        'gfaFile=s'      => \$gfa_file,
        'b7File=s'      => \$b7_file,
        'pafFile=s'      => \$paf_file,
        'out=s'          => \$out_prefix,
        'minPathNum=f'      => \$top_min,
        'minNodeLen=f'      => \$top_hits,
  );

  #$outp_name=basename($out_prefix);
  
  my $option_lows=1;
  my $err_code;
  
  if(! defined($gfa_file)){
    $option_lows = 0;
    $err_code='Fatal: No gfa file!';
  }
  elsif(! defined($mt_read_info)){
    $option_lows = 0;
    $err_code='Fatal: No file for mt-reads information specified!';
  }
  elsif((! defined($b7_file)) && (!defined($paf_file))){
    $option_lows = 0;
    $err_code='Fatal: No alignment file specified!';
  }
  else{
    ;
  }
  
  if (! $option_lows){
    print "\n  $err_code \n";
    exit(0);
  }
}

sub clip_graph{
  phaseOpt(@_);
  $out_prefix = 'gsat' if ! defined($out_prefix);
  warn "Reading the alignment file ... \n";

  if(defined($b7_file)){
    open(align_file,$b7_file) || die "Cannot open the file: $b7_file\n";
    chomp(my @mapped=<align_file>);
    readmap2('b7',\@mapped);
    splice @mapped;
    close align_file;
  }

  if(defined($paf_file)){
    open(align_file,$paf_file) || die "Cannot open the file: $paf_file\n";
    chomp(my @mapped=<align_file>);
    readmap2('paf',\@mapped);
    splice @mapped;
    close align_file;
  }

  open(infile1,$mt_read_info) || die "Cannot open the file: $mt_read_info\n";#mt_reads
  open(infile2,$gfa_file) || die "Cannot open the file: $gfa_file\n";#gfa_file

  my @content=<infile2>;
  #chomp(@content);
  @content = map {s/\s+$//;$_} @content;
  close infile2;
  ($gfa_S,undef,$gfa_L)=graphIO::readGfa(\@content,'SL');
  my %sel_rid;
  foreach my $line(<infile1>){
    chomp($line);
    my $rid=(split(/\t/,$line,2))[0];
    $sel_rid{$rid}=1;
  }
  close infile1;
  my $link_ovl=${$gfa_L}[0]->[-1];

  warn "Extracting the mapping information for mt-reads ... \n";
  my @sel_rno=grep {exists $sel_rid{$align_info[$_]->[$no_Qid]}} 0..$#align_info;

  #收集处理后比对结果，分contig进行断点评估
  my %ctg_break;
  foreach my $rno(@sel_rno){
    my $ctg_ali=$align_info[$rno]->[$no_Sid];
    next if $align_info[$rno]->[$no_Se] - $align_info[$rno]->[$no_Ss] < $min_len;
    $ctg_break{$ctg_ali}=[] if not exists $ctg_break{$ctg_ali};
    push @{$ctg_break{$ctg_ali}},[$align_info[$rno]->[$no_Ss],'l',$rno];
    push @{$ctg_break{$ctg_ali}},[$align_info[$rno]->[$no_Se],'r',$rno];
  }

  foreach my $ctg(keys %ctg_break){
    my %cn_freq;
    foreach my $cn(@{$ctg_break{$ctg}}){
      if(exists $cn_freq{$cn->[0]}){
        $cn_freq{$cn->[0]} ++;
      }
      else{
        $cn_freq{$cn->[0]} = 1;
      }
    }
    my @break_sort = map {[$_,$cn_freq{$_}]} (sort {$a <=> $b} keys %cn_freq);
    warn "Find break points for $ctg : ".@break_sort."\n";

    my @break_tops;
    my $mid_lk=$link_ovl * 0.5;
    my $mid_max= $mid_lk > 50 ? $mid_lk:50;
    my $start=1;
    my $end=$top_hits;
    my $step=1;
    my $top_count=0;
    my $topp_pre=0;
    my $topf_pre=0;
    my $run_num= ${$gfa_S}{$ctg}{'len'} - $mid_max + 1;
    #my $run_num=$break_sort[-1]->[0] - $top_hits + 1;
    #my @pos_tops;
    foreach my $rn (1..$run_num){
      #$start += $step * ($rn -1);
      #$end += $step * ($rn -1);
      my $max_p=0;
      my $max_f=0;
      foreach my $pos($start..$end){
        if(exists $cn_freq{$pos} && $cn_freq{$pos} > $max_f){
          $max_f = $cn_freq{$pos};
          $max_p = $pos;
        }
        else{
          ;
        }
      }
      #warn "Pos: $rn ; max_p: $max_p ; max_f: $max_f \n";
      if($max_p > 0 && $max_p == $topp_pre){
        $top_count++;
        push @break_tops,$max_p if $top_count == $top_hits && $max_f >= $top_min;
      }
      else{
        $top_count=1;
        $topp_pre=$max_p;
        $topf_pre=$max_f;
      }
      $start += $step;
      $end += $step;
    }
    warn "Find top-breaks for $ctg: @break_tops \n";
    
    #clip_ctg
    my $fragment_no=0;
    my $fs;
    my $fl=0;
    #my $fe=${$gfa_S}{$ctg}{'len'};
    #my $break_num=0;
    foreach my $fn(0..$#break_tops){
      my $fb_p=$break_tops[$fn];
      $fs = $fl;
      #my $mid_lk=$link_ovl * 0.5;
      #my $mid_max= $mid_lk > 50 ? $mid_lk:50;
      my $mid_int=$link_ovl+$mid_max;
      my $fnext;
      if($fn >= $#break_tops){
        $fnext=${$gfa_S}{$ctg}{'len'};
      }
      else{
        $fnext=$break_tops[$fn+1];
      }
      my $leftInt=$fb_p - $fs;
      my $rightInt=$fnext - $fb_p;
      my $skip=0;
#      $fl=$fb_p;

      if($fb_p < $mid_max ||  $leftInt < $mid_max || $rightInt < $mid_max || ${$gfa_S}{$ctg}{'len'} - $fb_p < $mid_max){
        ;
      }
      elsif($leftInt > $mid_int &&  $rightInt > $mid_int){
        $fragment_no ++;
        unshift @content,"S\t${ctg}_fragment${fragment_no}\t".substr(${$gfa_S}{$ctg}{'seq'},$fs,$leftInt)."\tDP:f:".${$gfa_S}{$ctg}{'dep'};
        $fragment_no ++;
        unshift @content,"S\t${ctg}_fragment${fragment_no}\t".substr(${$gfa_S}{$ctg}{'seq'},$fb_p-$link_ovl,$link_ovl * 2)."\tDP:f:".${$gfa_S}{$ctg}{'dep'};
        $fl=$fb_p;
      }
      elsif($rightInt > $mid_int){
        $fragment_no ++;
        unshift @content,"S\t${ctg}_fragment${fragment_no}\t".substr(${$gfa_S}{$ctg}{'seq'},$fs,$leftInt + $link_ovl)."\tDP:f:".${$gfa_S}{$ctg}{'dep'};
        $fl=$fb_p;
      }
      elsif($leftInt > $mid_int){
        $fragment_no ++;
        unshift @content,"S\t${ctg}_fragment${fragment_no}\t".substr(${$gfa_S}{$ctg}{'seq'},$fs,$leftInt)."\tDP:f:".${$gfa_S}{$ctg}{'dep'};
        $fl=$fb_p-$link_ovl;
      }
      else{
        ;
      }

      if($fn == $#break_tops && $fb_p < ${$gfa_S}{$ctg}{'len'}){
        $fragment_no ++;
        unshift @content,"S\t${ctg}_fragment${fragment_no}\t".substr(${$gfa_S}{$ctg}{'seq'},$fl)."\tDP:f:".${$gfa_S}{$ctg}{'dep'};
      }
    }

    if($fragment_no > 0){
      foreach (@content){
        s/(^S\t$ctg\t)/\#$1/;
      }

      foreach my $fno(1..$fragment_no-1){
        push @content,"L\t${ctg}_fragment${fno}\t+\t${ctg}_fragment".($fno+1)."\t+\t${link_ovl}M";
      }

      foreach my $c_line(@content){
        next if $c_line !~ /^L/;
        my @cs=split(/\t/,$c_line);
        if($cs[1] eq $ctg && $cs[2] eq '+'){
          $cs[1]=$ctg.'_fragment'.$fragment_no;
        }
        elsif($cs[1] eq $ctg && $cs[2] eq '-'){
          $cs[1]=$ctg.'_fragment1';
        }
        else{
          ;
        }
        if($cs[3] eq $ctg && $cs[4] eq '+'){
          $cs[3]=$ctg.'_fragment1';
        }
        elsif($cs[3] eq $ctg && $cs[4] eq '-'){
          $cs[3]=$ctg.'_fragment'.$fragment_no;
        }
        else{
          ;
        }
        $c_line=join("\t",@cs);
      }
    }
    warn 'Fragmental nodes after applying top breaks: '.$fragment_no."\n";
  }

  open(OUTG,">${out_prefix}.cleaved.gfa");
  foreach (@content){
    next if $_ !~ /^[SPL]/;
    print OUTG $_."\n";
  }
  close OUTG;
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
