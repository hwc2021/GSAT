#!perl
#Generated at Dec 23, 2025
#updated at Jan 03, 2026

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

sub map2matrix{
  my ($pathFile,$gfa)=@_;
  my %ctg_path;
  my @map_lk;
  my %matrix;
  open(pFILE,"$pathFile") || die $!;
  my $head=<pFILE>;
  foreach my $line(<pFILE>){
    chomp($line);
    my ($ctg,undef,undef,undef,$path)=split(/\t/,$line,6);
    my @t_paths=map {die "Error: Alternative path is not allowed!\n" if /\{/;s/[\(\[][^\)\]]+[\)\]]//g;$_} (split(/,/,$path));
    $ctg_path{$ctg}=[@t_paths];

    my @rc_lks=map {my $n1=$t_paths[$_];my $o1=substr($n1,-1,1,'');my $n2=$t_paths[$_+1];my $o2=substr($n2,-1,1,'');[$n1,$o1,$n2,$o2,0]} 0..$#t_paths-1 if @t_paths > 1;
    push @map_lk,@rc_lks;
  }
  close pFILE;

  open(gFile,"$gfa") || die $!;
  chomp(my @gfa_content=<gFile>);
  close gFile;
  my ($gfa_S,undef,$gfa_L)=graphIO::readGfa(\@gfa_content,'SL');
  my %t1_lk=graphIO::nodeNeighbors($gfa_L);

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

  open(cList,"$ctgList") || die $!;
  chomp(my @ctgs=sort <cList>);
  close cList;

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
  print outFile join("\t",('Link',@snames))."\n";
  foreach my $lnn(0..$#wholeM_loci){
    print outFile join("\t",($wholeM_loci[$lnn],@{$wholeM_value[$lnn]}))."\n" if $m_sum[$lnn] != 0 && $m_sum[$lnn] != @snames;
  }
  close outFile;
}

1;