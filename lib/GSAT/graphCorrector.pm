#updated in version 1.03: removed the gaps remained in the corrected gfa file.
#updated in version 1.031: fixed a big bug in the readGfa module.
#updated at Aug 19,2022: fixed a bug related to GraphIO module.
#updated at 20240430: fixed some bugs when record the variations.
#updated at 20250217: enhanced the correction effeciency on N-nodes.
#updated at 20250217: enable the multiprocess support.
#updated at 20250221: Adjust the criteria for determining the best corrected genotype in nstat-contig.
#updated at 20250305: fix a small bug.
#updated at 20250407: fix a small bug.
#updated at 20250731: try to reduce the memory usage.
#updated at 20250905: improve the alignment of nctgs.
#updated at 20260115: Adjust the criteria for correcting the old ONT reads
#注意skipHPM选项尚未在菜单中启用

package graphCorrector;

use strict;
#use warnings;
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";
use GSAT::graphIO;
use List::MoreUtils qw(any uniq);
use List::Util 'max';
use List::Util 'min';
use List::Util 'sum';
use File::Basename;
use Getopt::Long;
use Getopt::Long qw(GetOptionsFromArray);
use Parallel::ForkManager;
use Storable;
use File::Temp qw(tempdir);

my $read_file;
my $gfa_file;
my $path_file;
my $out_prefix;
my $outp_name;
my $true_ratio=0.9;
my $threads=8;
my $minPath=10;
my $hp_stat=0;

sub phaseOpt{
  GetOptionsFromArray($_[0],
        'readFile=s'     => \$read_file,
        'gfaFile=s'      => \$gfa_file,
        'pathFile=s'     => \$path_file,
        'out=s'          => \$out_prefix,
        'minReadProp=f'  => \$true_ratio,
        'skipHPM=s'      => \$hp_stat,  #for ONT reads   
  );

  $outp_name=basename($out_prefix);
  
  my $option_lows=1;
  my $err_code;
  
  if(! defined($gfa_file)){
    $option_lows = 0;
    $err_code='Fatal: No gfa file!';
  }
  elsif(! defined($read_file)){
    $option_lows = 0;
    $err_code='Fatal: No read file specified!';
  }
  elsif(! defined($path_file)){
    $option_lows = 0;
    $err_code='Fatal: No path file specified!';
  }
  else{
    ;
  }
  
  if (! $option_lows){
    print "\n  $err_code \n";
    exit(0);
  }
}

sub comp_rev{
  my $dna=shift @_;
  $dna =~ tr/atcgATCG/tagcTAGC/;
  my $seq = reverse($dna);

  return $seq;
}

sub call_var{
  my ($prid,$prseq,$mcn,$mrali,$mcali,$mori,$gs,$ctg_v,$ctg_d,$nctgs,$ncrseq,$t_dir) = @_;

  foreach (keys %{$gs}){
    ${$ctg_v}{$_}=[({})];
    ${$ctg_d}{$_}=[(0) x ${$gs}{$_}{'len'}];
  }
  foreach my $path_no(0..$#{$prid}){
    foreach my $this_ctg_no(0..$#{${$mcn}[$path_no]}){
      my $this_ctg=${$mcn}[$path_no]->[$this_ctg_no];
      my ($r_s,$r_e)=split(/,/,${$mrali}[$path_no]->[$this_ctg_no]);
      my ($c_s,$c_e)=split(/,/,${$mcali}[$path_no]->[$this_ctg_no]);
      foreach ($c_s-1 .. $c_e-1){
        ${$ctg_d}{$this_ctg}->[$_] ++;
      }
  
      my $c_seq=substr(${$gs}{$this_ctg}{'seq'},$c_s-1,$c_e-$c_s+1);
      my $r_seq;
      
      if(${$mori}[$path_no]->[$this_ctg_no] eq '-'){
        $r_seq=comp_rev(substr(${$prseq}[$path_no],$r_s-1,$r_e-$r_s+1));
      }
      else{
        $r_seq=substr(${$prseq}[$path_no],$r_s-1,$r_e-$r_s+1);
      }

      if(${$gs}{$this_ctg}{'nstat'} > 0){
        ${$nctgs}{$this_ctg}=${$gs}{$this_ctg}{'seq'} if not exists ${$nctgs}{$this_ctg};
        if(not exists ${$ncrseq}{$this_ctg}){
          ${$ncrseq}{$this_ctg}=[$r_seq];
        }
        else{
          push @{${$ncrseq}{$this_ctg}},$r_seq;
        }
        next;
      }

      # MAFFT alignment
      my ($temp_fh, $temp_fn) = File::Temp::tempfile(
        "${outp_name}_XXXXX", 
        DIR => $t_dir,
        SUFFIX => '.fas',
        UNLINK => 1
      );
      print $temp_fh ">ctg\n$c_seq\n>read\n$r_seq\n";
      close $temp_fh;
  
      my $mafft_command="mafft --maxiterate 1000 --globalpair --quiet $temp_fn";
      my @mafft_seq=`$mafft_command`;
      die "Fatal: the mafft software is not well installed or the input files are not matched each other. \n" if @mafft_seq < 4;
      my (undef,$ali_seq)=graphIO::readFas(\@mafft_seq);
      my $ali_c_seq=uc($ali_seq->[0]);
      my $ali_r_seq=uc($ali_seq->[1]);
      $c_seq=uc($c_seq);

      record_var($this_ctg,$c_seq,$c_s,$ali_c_seq,[$ali_r_seq],$ctg_v);
    }
  }
}

sub record_var{
  my ($thisctg,$cseq,$cs,$alicseq,$alirseq_in,$ctgvar)=@_;
#  my $nlarge = 20;#后期可以设为调整参数
#  my @var_large;

  #尝试变异检测和保存
  my @ref=split(//,uc($cseq));
  $alicseq = uc($alicseq);
  foreach my $alirseq(@{$alirseq_in}){
    my @alirseqs=split(//,uc($alirseq));
    my $ali_s=-1;
    my $ali_e=-1;
    my $ref_no=0;
    my $loc0=$ref_no+$cs-1;
    
    while($alicseq =~ /(\-*[atcgnATCGN])/g){
      my $now_cseq=$1;
      my $now_clen=length($now_cseq);
      $ali_s ++;
      $ali_e += $now_clen;    
      my $ref_nuc=$ref[$ref_no];
      die"Fatal: The reference and the mafft results is not matched [$ref_nuc : $1].\n" if $1 !~ /$ref_nuc/;
  
      my $alt_nuc=join('',@alirseqs[$ali_s .. $ali_e]);
      #去除gap符号
      if(length($alt_nuc) > 1 && $alt_nuc =~ s/\-//g){
        $alt_nuc = '-' if length($alt_nuc) < 1;
      }

#      if($now_clen >= $nlarge){
#        push @var_large,{'loc0' => $loc0};
#        my $varl_cseq=$now_cseq;
#        my $varl_cseq =~ s/\-/N/g;
#        record_var($thisctg,$varl_cseq,$loc0+1,$varl_cseq,[$alt_nuc],$var_large[-1]);
#        #next;
#      }
#      else{
        if($alt_nuc ne $ref_nuc){
          if(not exists ${${$ctgvar}{$thisctg}->[$loc0]}{$alt_nuc}){
            ${${$ctgvar}{$thisctg}->[$loc0]}{$alt_nuc} = 1;
          }
          else{
            ${${$ctgvar}{$thisctg}->[$loc0]}{$alt_nuc} ++;
          }
        }
#      }
        
      $ref_no ++;
      $loc0 ++;
      $ali_s=$ali_e;
    }
  }

#  foreach my $var_l(@var_large){
#    my $varl_total=@{$alirseq_in};
#    my $alt_seq=map {my $x=$_;my @corr=grep {${$x}{$_}/$varl_total > $true_ratio} (keys %{$x});if(@corr == 1){return $corr[0];}else{return 'N';}} @{${$var_l}{$thisctg}};
#    #$alt_seq =~ s/N//g;#默认取消该设置，后续可以作为输入参数
#    if(length($alt_seq) > 0){
#      my $true_loc0=${$var_l}{'loc0'};
#      ${${$ctgvar}{$thisctg}->[$true_loc0]}{$alt_seq}=$varl_total;
#    }
#  }
}

sub grapCorr{
  phaseOpt(@_);
  $out_prefix = basename($gfa_file) if ! defined($out_prefix);

  warn "Reading gfa file ...\n";
  open(gfaFile,"$gfa_file") || die "Error: Cannot open the file: $gfa_file \n";
  chomp(my @gfa_content=<gfaFile>);
  
  my ($gfa_S,undef,$gfa_L)=graphIO::readGfa(\@gfa_content,'SL');#readFas
  my $gfa_lnk_k=${$gfa_L}[0]->[-1];
  close gfaFile;
  warn "Reading path file ...\n";
  open(pathFile,"$path_file") || die "Error: Cannot open the file: $path_file \n";
  
  my @path_Rid;
  my @path_Rlen;
  my @mapped_ctg;#二维数组
  my @mapped_ori;#二维数组
  my @mapped_Rali;#二维数组
  my @mapped_Cali;#二维数组
  my %path_Rseq;
  my $not_use=<pathFile>;
  foreach (<pathFile>){
    chomp;
    next if length($_) == 0;
    my ($Rid,$Rlen,undef,undef,undef,undef,undef,$ctg,$ori,$Rali,$Cali)=split(/\t/);
    push @path_Rid,$Rid;
    push @path_Rlen,$Rlen;
    push @mapped_ctg,[(split(/;/,$ctg))];
    push @mapped_ori,[(split(/;/,$ori))];
    push @mapped_Rali,[(split(/;/,$Rali))];
    push @mapped_Cali,[(split(/;/,$Cali))];
  }
  close pathFile;
  warn "Reading reads file ...\n";
  open(readFile,"$read_file") || die "Error: Cannot open the file: $read_file \n";
  my %search_list = map { $_ => 1 } @path_Rid;  # 使用哈希快速查找
  print "Searching reads for paths: ".@path_Rid."\n";
  #my %path_Rseq;
  my $current_id;
  while (<readFile>) {
    chomp;
    if (/^>(\S+)/) {
        $current_id = $1;
        $path_Rseq{$current_id} = '' if exists $search_list{$current_id};
    } 
    elsif ($current_id && exists $search_list{$current_id}) {
        $path_Rseq{$current_id} .= $_;
    }
  }
  close readFile;
  my $r_count=keys %path_Rseq;
  print "Obtained reads for paths: ".$r_count."\n";
  die "Error: Sequnce labels in the read file and gfa file are not matched!\n" if @path_Rid != $r_count;
  
  #如果只矫正ctg的序列
  #my @ctg_var;#hash-'alt'=>'count'
  my %ctg_var;
  my %ctg_dep;
  foreach (keys %{$gfa_S}){
    $ctg_var{$_}=[({})];
    $ctg_dep{$_}=[(0) x ${$gfa_S}{$_}{'len'}];
  }
  my %nstat_ctgs;
  my %nstat_ctgrseq;

  # Create temporary directory
  my $tempdir = tempdir("${out_prefix}_XXXXX", CLEANUP => 1);

  # Split tasks into chunks
  my @all_paths = (0..$#path_Rid);
  my $chunk_size = int(@all_paths/$threads) + 1;
  my @index_chunks;
  while (@all_paths) {
    push @index_chunks, [splice @all_paths, 0, $chunk_size];
  }

  # Parallel processing
  my $pm = Parallel::ForkManager->new($threads);
  #my @result_files;
  my @result_ctgv;  # 存储子进程结果
  my @result_ctgd;
  my @result_nctgs;
  my @result_nrseqs;
  # 设置子进程结束时的回调
  $pm->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $signal, $core_dump, $data) = @_;
    return unless defined $data;
    my ($chunk_id, $ctgv, $ctgd, $nctgs, $nrseqs) = @$data;
    # 存储结果
    $result_ctgv[$chunk_id] = $ctgv;
    $result_ctgd[$chunk_id] = $ctgd;
    $result_nctgs[$chunk_id] = $nctgs;
    $result_nrseqs[$chunk_id] = $nrseqs;
  });

  foreach my $chunk_id (0..$#index_chunks) {
    $pm->start and next;
    
    # OPTIMIZED: 动态获取分块数据
    my $indices_ref = $index_chunks[$chunk_id];
    my @chunk_prid = @path_Rid[@$indices_ref];
    my @chunk_prseq = map { $path_Rseq{$_} } @chunk_prid;
    my @chunk_mcn = @mapped_ctg[@$indices_ref];
    my @chunk_mrali = @mapped_Rali[@$indices_ref];
    my @chunk_mcali = @mapped_Cali[@$indices_ref];
    my @chunk_mori = @mapped_ori[@$indices_ref];
    
    # 子进程局部变量
    my %chunk_ctgv;
    my %chunk_ctgd;
    my %chunk_nctgs;
    my %chunk_nrseqs;

    call_var(
        \@chunk_prid,
        \@chunk_prseq,
        \@chunk_mcn,
        \@chunk_mrali,
        \@chunk_mcali,
        \@chunk_mori,
        $gfa_S,
        \%chunk_ctgv,
        \%chunk_ctgd,
        \%chunk_nctgs,
        \%chunk_nrseqs,
        $tempdir
    );
    
    # OPTIMIZED: 使用临时变量传递结果
    $pm->finish(0, [
        $chunk_id, 
        \%chunk_ctgv, 
        \%chunk_ctgd,
        \%chunk_nctgs,
        \%chunk_nrseqs
    ]);
  }

  $pm->wait_all_children;
  
  undef @path_Rid;
  undef @mapped_ctg;
  undef @mapped_ori;
  undef @mapped_Rali;
  undef @mapped_Cali;

  warn "Merging the resulting data ... \n";

  foreach my $t_ctg (keys %{$gfa_S}){
    foreach my $t_pos(0..$#{$ctg_dep{$t_ctg}}){
      foreach my $t_cid (0..$#result_ctgv){
        my $chunk_data = $result_ctgv[$t_cid];
        foreach my $t_alt(keys %{$chunk_data->{$t_ctg}->[$t_pos]}){
          if(exists ${$ctg_var{$t_ctg}->[$t_pos]}{$t_alt}){
            ${$ctg_var{$t_ctg}->[$t_pos]}{$t_alt} += ${$chunk_data->{$t_ctg}->[$t_pos]}{$t_alt};
          }
          else{
            ${$ctg_var{$t_ctg}->[$t_pos]}{$t_alt} = ${$chunk_data->{$t_ctg}->[$t_pos]}{$t_alt};
          }
        }
      }
    }
    $ctg_dep{$t_ctg}= [(map {my $t_p=$_;my @ctgds=map {${$result_ctgd[$_]}{$t_ctg}->[$t_p]} 0..$#result_ctgd;sum(@ctgds)} 0..$#{$ctg_dep{$t_ctg}})];
  }

  %nstat_ctgs = %{shift @result_nctgs};
  %nstat_ctgrseq = %{shift @result_nrseqs};
  foreach my $c_no(0..$#result_nctgs){
    my %merge_ncs=(%nstat_ctgs,%{$result_nctgs[$c_no]});
    %nstat_ctgs = %merge_ncs;

    foreach my $t_nctg(keys %{$result_nrseqs[$c_no]}){
      if(exists $nstat_ctgrseq{$t_nctg}){
        my @merge_ncr=(@{$nstat_ctgrseq{$t_nctg}},@{${$result_nrseqs[$c_no]}{$t_nctg}});
        $nstat_ctgrseq{$t_nctg}=[@merge_ncr];
      }
      else{
        $nstat_ctgrseq{$t_nctg} = ${$result_nrseqs[$c_no]}{$t_nctg};
      }
    }
  }

  foreach my $nctg(keys %nstat_ctgs){
    my $nc_seq=$nstat_ctgs{$nctg};
    open temp_fas,">${out_prefix}.temp.fas";
    print temp_fas ">ctg\n".$nc_seq."\n";
    foreach my $t_no(0..$#{$nstat_ctgrseq{$nctg}}){
      print temp_fas ">read$t_no\n".$nstat_ctgrseq{$nctg}->[$t_no]."\n";
    }
    close temp_fas;

    my $mafft_command="mafft --auto --quiet ${out_prefix}.temp.fas";
    my @mafft_seq=`$mafft_command`;
    die "Fatal: the mafft software is not well installed or the input files are not matched each other. \n" if @mafft_seq < 4;
    my (undef,$ali_seq)=graphIO::readFas(\@mafft_seq);
    my $alic_seq=shift @{$ali_seq};

    $nc_seq=~/N/;
    my $n1st_no=$-[0];
    my $n1st_dep=$ctg_dep{$nctg}->[$n1st_no];

    #将最近的1个gap挪进N-array里去
    $alic_seq =~ s/([\-]+)([atcgATCG]{1,3}[nN\-]*[nN][nN\-]*)/$2$1/;
    $alic_seq =~ s/([nN\-]*[nN][nN\-]*[atcgATCG]{1,3})([\-]+)/$2$1/;

    $alic_seq =~ s/([nN\-]*[nN][nN\-]*)/'N' x length($1)/ge;
    $nc_seq=$alic_seq;
    $nc_seq=~s/\-//g;
    ${$gfa_S}{$nctg}{'seq'}=uc($nc_seq);
    my $now_len=length($nc_seq);
    my $add_len=$now_len-${$gfa_S}{$nctg}{'len'};
    ${$gfa_S}{$nctg}{'len'}=$now_len;

    my $c_dep=$ctg_dep{$nctg};
    $ctg_dep{$nctg}=[@{$c_dep}[0..$n1st_no],($n1st_dep) x $add_len, @{$c_dep}[$n1st_no+1..$#{$c_dep}]];

    #由于末端的位置难以简单界定统一，因此不对末端进行矫正，设置忽略末端比对的长度，后续可以作为参数调用
    my $non_corr_end_len=20;
    #my $non_corr_start0=$non_corr_end_len-1;
    my $c_start=1;#如果不设置末端长度，这里的初始值应当为1
    my $del_seq=substr($alic_seq,0,$non_corr_end_len);
    my $new_ali_len=length($alic_seq)-$non_corr_end_len*2;
    $alic_seq=substr($alic_seq,$non_corr_end_len,$new_ali_len);
    my @alir_seqs=map {substr($_,$non_corr_end_len,$new_ali_len)} @{$ali_seq};

    $del_seq =~ s/\-//g;
    my $del_len=length($del_seq);
    my $nnc_ctg=substr($nc_seq,$del_len,$now_len-$del_len*2);
    $c_start=$del_len+1;

    record_var($nctg,$nnc_ctg,$c_start,$alic_seq,\@alir_seqs,\%ctg_var);
  }
  
  #开始矫正；0，记录无变异；>true_count_min,直接矫正；>count_min,报警组装错误；否则仅记录稀有变异。输出变异记录文件，可以包括深度信息
  #记录文件：ctg\tpos\tref\tread_count\tvar_count\tstat[corrected|detected|identical]\talt_nuc[;]\talt_count[;]
  open out_gfa,">${out_prefix}.corrCtg.gfa";
  open out_record,">${out_prefix}.CtgVar.info";
  print out_record "ctg\tpos\tref\tread_count\tvar_count\tstat\talt_nuc\talt_count\n";
  my %corr_c_seq;
  foreach my $ctg(keys %ctg_dep){
    my @corr_seq;
    my $t_cseq=${$gfa_S}{$ctg}{'seq'};
    my @ctg_ref=split(//,$t_cseq);
    my $corr_seq;
    my %hp;
    if($hp_stat > 0){
      #鉴定同聚物区域
      my $minLen=5;
      while($t_cseq =~ /(A{$minLen,}|T{$minLen,}|C{$minLen,}|G{$minLen,})/gi){
        foreach my $ii($-[0] .. $+[0]){
          $hp{$ii}=1;
        }
      }
    }

    foreach my $loc0(0..$#{$ctg_dep{$ctg}}){
      my @alt_type=keys %{$ctg_var{$ctg}->[$loc0]};
      my @alt_count=map {${$ctg_var{$ctg}->[$loc0]}{$_}} @alt_type;
      my $dep_count=$ctg_dep{$ctg}->[$loc0];
      my @true_no=grep {$dep_count > $minPath && $alt_count[$_]/$dep_count >= $true_ratio} 0..$#alt_count;
      my $var_sum=int(sum(@alt_count));
      my $stat;
      my $ref0=$ctg_ref[$loc0];
  
      if(@true_no >1){
        warn "$ctg\t$loc0\t".$ref0."\t$dep_count\t$var_sum\tundefined\t".join(';',@alt_type)."\t".join(';',@alt_count)."\n";
        die "Error: Found two or more genotypes which should be corrected in the same locus. Please check your min_ratio_corr option!\n";
      }
      elsif(@true_no == 1){
        #判断是否位于同聚物区域
        my $alt_ss = $alt_type[$true_no[0]];
        if($hp_stat > 0 && exists $hp{$loc0} && ($alt_ss eq '-' || $alt_ss =~ /^[$ref0]+$/)){
          $stat = 'skipped';
          $corr_seq .= $alt_ss;
        }
        else{
          $stat = 'corrected';
          $corr_seq .= $alt_type[$true_no[0]];
        }
      }
      elsif($var_sum == 0){
        $stat = 'identical';
        $corr_seq .= $ref0;
      }
      elsif(${$gfa_S}{$ctg}{'nstat'} > 0 && $ref0 =~ /N|n/){
        my $max_c=max(@alt_count);
        my @max_alt_no=grep {$alt_count[$_] == $max_c} 0..$#alt_count;
        $corr_seq .= $alt_type[$max_alt_no[0]];
        $stat = 'corrected_low_qual';
      }
      else{
        $stat = 'detected';
        $corr_seq .= $ref0;
      }
      print out_record "$ctg\t$loc0\t".$ref0."\t$dep_count\t$var_sum\t$stat\t".join(';',@alt_type)."\t".join(';',@alt_count)."\n";
    }
    $corr_seq =~ s/\-//g;
    $corr_c_seq{$ctg}=$corr_seq;
  }
  close out_record;
  
  my @temp_content=grep {$_ !~ /^S/} @gfa_content;
  my @S_content=map {"S\t".$_."\t".$corr_c_seq{$_}."\tDP:f:".${$gfa_S}{$_}{'dep'}} (sort keys %ctg_dep);
  print out_gfa join("\n",(@S_content,@temp_content))."\n";
  close out_gfa;
}

1;