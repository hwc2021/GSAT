#!perl
#Generated at Dec 23, 2025
#updated at Jan 03, 2026
#updated at Jul 16, 2026
#updated at Jul 21, 2026

package graphMatrix;

use strict;
use warnings;
#use List::MoreUtils qw(any uniq);
use List::Util qw(max min sum);
#use List::Util 'min';
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/lib";
use GSAT::graphIO;

my @o_x=(['+','+'],['+','-'],['-','+'],['-','-']);
my %ori_r=('+','-','-','+');
my $minEdgeRatio=0.5;#如果是针对core-segment，应启用该设置
my %ctg_info;

sub map2matrix{
  my ($pathFile,$gfa)=@_;
  my %ctg_path;
  my @map_lk;
  my %matrix;
  
  open(gFile,"$gfa") || die $!;
  chomp(my @gfa_content=<gFile>);
  close gFile;
  my ($gfa_S,undef,$gfa_L)=graphIO::readGfa(\@gfa_content,'SL');
  my %t1_lk=graphIO::nodeNeighbors($gfa_L);

  open(pFILE,"$pathFile") || die $!;
  my $head=<pFILE>;
  foreach my $line(<pFILE>){
    chomp($line);
    my ($ctg,undef,undef,undef,$path,$align1,$align2)=split(/\t/,$line);
    my @t_paths=map {die "Error: Alternative path is not allowed!\n" if /\{/;s/[\(\[][^\)\]]+[\)\]]//g;my $n=$_;my $o=substr($n,-1,1,'');[$n,$o]} (split(/,/,$path));
    my ($a1_s,$a1_e)=split(/\-/,$align1);
    shift @t_paths if ($a1_e - $a1_s) / $ctg_info{$t_paths[0]->[0]} < $minEdgeRatio;
    my ($a2_s,$a2_e)=split(/\-/,$align2);
    pop @t_paths if @t_paths > 0 && ($a2_e - $a2_s) / $ctg_info{$t_paths[-1]->[0]} < $minEdgeRatio;
    $ctg_path{$ctg}= [map{$_->[0].$_->[1]} @t_paths];

    my @rc_lks=map {[@{$t_paths[$_]},@{$t_paths[$_ + 1]},0]} 0..$#t_paths-1 if @t_paths > 1;
    push @map_lk,@rc_lks;
  }
  close pFILE;

  #检查被漏掉的ctg，并重构link矩阵
  my @lost=grep {not exists $ctg_path{$_}} (keys %{$gfa_S});
  my @lk_new;
  my %t_lk;
  if(@lost > 0){
    foreach my $tt(@{$gfa_L}){
      my $tn1=grep {$tt->[0] eq $_ || $tt->[2] eq $_} @lost;
      push @lk_new,$tt if $tn1 < 1;
    }

    foreach my $tx(@lost){
      foreach my $ii(@{$t1_lk{$tx}{'+'}{'left'}}){
        my $tln=$ii;
        my $tlo=substr($tln,-1,1,'');
        foreach my $jj(@{$t1_lk{$tx}{'+'}{'right'}}){
          my $trn=$jj;
          my $tro=substr($trn,-1,1,'');
          push @lk_new,[$tln,$tlo,$trn,$tro,0];
        }
      }
    }
    %t_lk=graphIO::nodeNeighbors(\@lk_new);
  }
  else{
    %t_lk=%t1_lk;
  }
  undef %t1_lk;

  foreach my $node(keys %t_lk){
    foreach my $li(@{$t_lk{$node}{'+'}{'left'}}){
      my $ln=$li;
      my $lo=substr($ln,-1,1,'');
      my $ln2;
      my $lo2;
      if($lo eq '+'){
        $ln2=$ctg_path{$ln}->[-1];
        $lo2=substr($ln2,-1,1,'');
      }
      elsif($lo eq '-'){
        $ln2=$ctg_path{$ln}->[0];
        $lo2=$ori_r{substr($ln2,-1,1,'')};
      }
      else{
        die"Error: Fatal error!\n";
      }

      my $nn=$ctg_path{$node}->[0];
      my $no=substr($nn,-1,1,'');
      push @map_lk,[$ln2,$lo2,$nn,$no,0];
    }
    foreach my $ri(@{$t_lk{$node}{'+'}{'right'}}){
      my $rn=$ri;
      my $ro=substr($rn,-1,1,'');
      my $rn2;
      my $ro2;
      if($ro eq '+'){
        $rn2=$ctg_path{$rn}->[0];
        $ro2=substr($rn2,-1,1,'');
      }
      elsif($ro eq '-'){
        $rn2=$ctg_path{$rn}->[-1];
        $ro2=$ori_r{substr($rn2,-1,1,'')};
      }
      else{
        die"Error: Fatal error!\n";
      }

      my $nn=$ctg_path{$node}->[-1];
      my $no=substr($nn,-1,1,'');
      push @map_lk,[$nn,$no,$rn2,$ro2,0];
    }
  }

  my %all_lks=graphIO::nodeNeighbors(\@map_lk);

  my @t_ctg=sort keys %all_lks;

  foreach my $ti(0..$#t_ctg){
    my $ictg=$t_ctg[$ti];
    $matrix{$ictg.'-'}{$ictg.'-'} = 1;
    foreach my $tj($ti..$#t_ctg){
      my $jctg=$t_ctg[$tj];
      foreach my $oo(@o_x){
        foreach my $tx(@{$all_lks{$ictg}{$oo->[0]}{'right'}}){
          my $testj=$jctg.$oo->[1];
          $matrix{$ictg.$oo->[0]}{$testj} = 1 if $tx eq $testj;
        }
      }
    }
  }

  return %matrix;
}

sub cmpMatrix{
  my ($fileList,$ctgList,$outp)=@_;

  open(cList,"$ctgList") || die $!;
  foreach (<cList>){
    chomp;
    my ($n,$l)=split(/\s+/);
    die "Error: The contigList file should contain two colums: contig_name and contig_length. \n" if not defined($l);
    $ctg_info{$n}=$l;
  }
  close cList;

  open(fList,"$fileList") || die $!;
  my @snames;
  my @smatrixs;

  foreach my $fPair(<fList>){
    chomp($fPair);
    my ($nn,$pfile,$gfile)=split(/\s+/,$fPair);
    die 'Error: The sample name, names of paths-file (*mapping.paths) and gfa-file(*.gfa) were required for each sample: '."$nn -- $pfile -- $gfile\n" if $pfile !~ /mapping\.paths$/ || $gfile !~ /\.gfa$/;
    push @snames,$nn;
    my %t_m=map2matrix($pfile,$gfile);
    push @smatrixs,{%t_m};
  }
  close fList;

  my @ctgs=sort keys %ctg_info;

  #my @wholeM_name=@snames;
  my @wholeM_loci;
  my @wholeM_value;
  my @m_sum;
  my @ox2=@o_x;
  #push @ox2,['-','-'] if $pav_code > 1;

  foreach my $i(0..$#ctgs){
    my $ic=$ctgs[$i];
    foreach my $j($i..$#ctgs){
      my $jc=$ctgs[$j];        
      foreach my $z(@ox2){
        push @wholeM_loci,$ic.$z->[0].'#'.$jc.$z->[1];
        my @t_v;
        my $sum=0;
        foreach my $ss(0..$#snames){
          my $ssm=$smatrixs[$ss];
          if(exists ${$ssm}{$ic.$z->[0]}{$jc.$z->[1]}){
            push @t_v,${$ssm}{$ic.$z->[0]}{$jc.$z->[1]};
          }
          else{
            push @t_v,0;
          }
          $sum += $t_v[-1];
        }
        push @wholeM_value,[@t_v];
        push @m_sum,$sum;
      }
    }
  }

  open(outFile,">${outp}.diff.matrix");
  open(alterFile,">${outp}.common.links");
  print outFile join("\t",('Link',@snames))."\n";
  foreach my $lnn(0..$#wholeM_loci){ 
    if($m_sum[$lnn] == @snames){
      print alterFile $wholeM_loci[$lnn]."\n";
    }
    elsif($m_sum[$lnn] != 0){
      print outFile join("\t",($wholeM_loci[$lnn],@{$wholeM_value[$lnn]}))."\n";
    }    
    else{
      ;
    }
  }
  close outFile;
  close alterFile;
}

1;
