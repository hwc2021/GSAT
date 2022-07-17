#updated in version 1.03: removed the gaps remained in the corrected gfa file.
#updated in version 1.031: fixed a big bug in the readGfa module
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

my $read_file;
my $gfa_file;
my $path_file;
my $out_prefix;
my $true_ratio=0.9;
sub phaseOpt{
  GetOptionsFromArray($_[0],
        'readFile=s'     => \$read_file,
        'gfaFile=s'      => \$gfa_file,
        'pathFile=s'      => \$path_file,
        'out=s'          => \$out_prefix,
        'minReadProp=f'      => \$true_ratio,
  );
  
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
  $/='';
  my $r_info=<readFile>;
  $/="\n";
  @_=split(/\n/,$r_info);
  undef $r_info;
  close readFile;
  
  warn "Constructing index for reads  ...\n";
  my ($r_id,$r_seq)=graphIO::readFas(\@_);
  my $r_count=0;
  my @search_list=@path_Rid;
  print "Searching reads for paths: ".@search_list."\n";
  foreach my $i(0..$#{$r_id}){
    my $r_label=${$r_id}[$i];
    $r_label=(split(/\s+/,$r_label))[0];
    my @r_hits=grep {$r_label eq $search_list[$_]} 0..$#search_list;
    next if @r_hits < 1;
    die "Error: Sequnce labels in the read file and gfa file are not matched!\n" if @r_hits > 1;
    $path_Rseq{$search_list[$r_hits[0]]}=${$r_seq}[$i];
    splice(@search_list,$r_hits[0],1);
    
    $r_count ++;
  }
  splice(@{$r_id});
  splice(@{$r_seq});
  print "Obtained reads for paths: ".$r_count."\n";
  close readFile;
  
  #如果只矫正ctg的序列
  #my @ctg_var;#hash-'alt'=>'count'
  my %ctg_var;
  my %ctg_dep;
  foreach (keys %{$gfa_S}){
    $ctg_var{$_}=[({}) x ${$gfa_S}{$_}{'len'}];
    $ctg_dep{$_}=[(0) x ${$gfa_S}{$_}{'len'}];
  }
  foreach my $path_no(0..$#path_Rid){
    foreach my $this_ctg_no(0..$#{$mapped_ctg[$path_no]}){
      my $this_ctg=$mapped_ctg[$path_no]->[$this_ctg_no];
      my ($r_s,$r_e)=split(/,/,$mapped_Rali[$path_no]->[$this_ctg_no]);
      my ($c_s,$c_e)=split(/,/,$mapped_Cali[$path_no]->[$this_ctg_no]);
      foreach ($c_s-1 .. $c_e-1){
        $ctg_dep{$this_ctg}->[$_] ++;
      }
  
      my $c_seq=substr(${$gfa_S}{$this_ctg}{'seq'},$c_s-1,$c_e-$c_s+1);
      my $r_seq;
      
      if($mapped_ori[$path_no]->[$this_ctg_no] eq '-'){
        $r_seq=comp_rev(substr($path_Rseq{$path_Rid[$path_no]},$r_s-1,$r_e-$r_s+1));
      }
      else{
        $r_seq=substr($path_Rseq{$path_Rid[$path_no]},$r_s-1,$r_e-$r_s+1);
      }
  
      open temp_fas,">${out_prefix}.temp.fas";
      print temp_fas ">ctg\n".$c_seq."\n>read\n".$r_seq."\n";
      close temp_fas;
  
      my $mafft_command="mafft --maxiterate 1000 --globalpair --quiet ${out_prefix}.temp.fas";
      my @mafft_seq=`$mafft_command`;
      die "Fatal: the mafft software is not well installed or the input files are not matched each other. \n" if @mafft_seq < 4;
      my (undef,$ali_seq)=graphIO::readFas(@mafft_seq);
      my $ali_c_seq=uc($ali_seq->[0]);
      my $ali_r_seq=uc($ali_seq->[1]);
      $c_seq=uc($c_seq);
  
      #尝试变异检测和保存
      my @ref=split(//,$c_seq);
      my $ref_no=0;
      my $loc0=$ref_no+$c_s-1;
      my @ali_r_seqs=split(//,$ali_r_seq);
      my $ali_s=-1;
      my $ali_e=-1;
      while($ali_c_seq =~ /(\-*[atcgnATCGN])/g){
        $ali_s ++;
        $ali_e += length($1);
        my $ref_nuc=$ref[$ref_no];
        die"Fatal: The reference and the mafft results is not matched [$ref_nuc : $1].\n" if $1 !~ /$ref_nuc/;
  
        my $alt_nuc=join('',@ali_r_seqs[$ali_s .. $ali_e]);
        if($alt_nuc ne $ref_nuc){
          if(not exists ${$ctg_var{$this_ctg}->[$loc0]}{$alt_nuc}){
           $ctg_var{$this_ctg}->[$loc0]={$alt_nuc => '1'};
          }
          else{
            ${$ctg_var{$this_ctg}->[$loc0]}{$alt_nuc} ++;
          }
        }
        
        $ref_no ++;
        $loc0 ++;
        $ali_s=$ali_e;
      }
    }
  }
  
  #开始矫正；0，记录无变异；>true_count_min,直接矫正；>count_min,报警组装错误；否则仅记录稀有变异。输出变异记录文件，可以包括深度信息
  #记录文件：ctg\tpos\tref\tread_count\tvar_count\tstat[corrected|detected|identical]\talt_nuc[;]\talt_count[;]
  open out_gfa,">${out_prefix}.corrCtg.gfa";
  open out_record,">${out_prefix}.CtgVar.info";
  print out_record "ctg\tpos\tref\tread_count\tvar_count\tstat\talt_nuc\talt_count\n";
  my %corr_c_seq;
  foreach my $ctg(keys %ctg_dep){
    my @corr_seq;
    my @ctg_ref=split(//,${$gfa_S}{$ctg}{'seq'});
    my $corr_seq;
    foreach my $loc0(0..$#{$ctg_dep{$ctg}}){
      my @alt_type=keys %{$ctg_var{$ctg}->[$loc0]};
      my @alt_count=map {${$ctg_var{$ctg}->[$loc0]}{$_}} @alt_type;
      my $dep_count=$ctg_dep{$ctg}->[$loc0];
      my @true_no=grep {$alt_count[$_]/$dep_count >= $true_ratio} 0..$#alt_count;
      my $var_sum=int(sum(@alt_count));
      my $stat;
  
      if(@true_no >1){
        die "Error: Found two or more genotypes which should be corrected in the same locus. Please check your min_ratio_corr option!\n";
      }
      elsif(@true_no == 1){
        $stat = 'corrected';
        $corr_seq .= $alt_type[$true_no[0]];
      }
      elsif($var_sum == 0){
        $stat = 'identical';
        $corr_seq .= $ctg_ref[$loc0];
      }
      else{
        $stat = 'detected';
        $corr_seq .= $ctg_ref[$loc0];
      }
      print out_record "$ctg\t$loc0\t".$ctg_ref[$loc0]."\t$dep_count\t$var_sum\t$stat\t".join(';',@alt_type)."\t".join(';',@alt_count)."\n";
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