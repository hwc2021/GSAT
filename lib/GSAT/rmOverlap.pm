package rmOverlap;
use strict;
use GSAT::graphIO;

#Adjust the overlap length to 0 to generate a valid genome graph, e.g., MMG, for annoation and further utilization
my @side=("right","left");
my %ori2side;
$ori2side{'left'}{'+'}='right';
$ori2side{'left'}{'-'}='left';
$ori2side{'right'}{'+'}='left';
$ori2side{'right'}{'-'}='right';
my 	%links;
my %m0_stat;#1-remain;-1-m0
sub m0proc{
	my ($in,$out_p)=@_;
	open (gfafile,$in) or die"Error: No gfa file detected: ".$ARGV[0]." !";
	chomp(my @gfa_content=<gfafile>); 
	my ($gfaS,undef,$gfaL)=graphIO::readGfa(\@gfa_content,"SL");
	
	my %mlen;
	my @ctgs=sort keys %{$gfaS};
	foreach my $line(@{$gfaL}){
		#print join("\t",@{$line})."\n";#test
		#第一个ctg
		if($line->[1] eq '+'){
			if(exists $links{$line->[0]}{'right'}){
				push @{$links{$line->[0]}{'right'}},$line->[2].$line->[3];
			}
			else{
				$links{$line->[0]}{'right'}=[($line->[2].$line->[3])];
				#print "test1: ".join(",",@{$links{$line->[0]}{'right'}})."\n";#test
				$mlen{$line->[0]}{'right'}=$line->[4];
			}
		}
		elsif($line->[1] eq '-'){
			my $l_ori=$line->[3];
			$l_ori =~ tr/\+\-/\-\+/;
			if(exists $links{$line->[0]}{'left'}){
				push @{$links{$line->[0]}{'left'}},$line->[2].$l_ori;
			}
			else{
				$links{$line->[0]}{'left'}=[$line->[2].$l_ori];
				$mlen{$line->[0]}{'left'}=$line->[4];
			}
		}
		else{
			die"Error: Format error in the GFA file! [ ".$line->[1]." ] \n";
		}
	
		#第二个ctg
		if($line->[3] eq '+'){
			if(exists $links{$line->[2]}{'left'}){
				push @{$links{$line->[2]}{'left'}},$line->[0].$line->[1];
			}
			else{
				$links{$line->[2]}{'left'}=[$line->[0].$line->[1]];
				$mlen{$line->[2]}{'left'}=$line->[4];
			}
		}
		elsif($line->[3] eq '-'){
			my $l_ori=$line->[1];
			$l_ori =~ tr/\+\-/\-\+/;
			if(exists $links{$line->[2]}{'right'}){
				push @{$links{$line->[2]}{'right'}},$line->[0].$l_ori;
			}
			else{
				$links{$line->[2]}{'right'}=[$line->[0].$l_ori];
				$mlen{$line->[2]}{'right'}=$line->[4];
			}
		}
		else{
			die"Error: Format error in the GFA file!\n";
		}
	}
	
	my @ctgs_sort=sort {${$gfaS}{$b}{'len'} <=> ${$gfaS}{$a}{'len'}} @ctgs;
	#print "sorted ctg: @ctgs_sort \n";
	
	foreach my $ctg(@ctgs_sort){
		foreach my $side1(@side){
			if(exists $m0_stat{$ctg}{$side1}){
				process_side_ctg($ctg,$side1);
			}
			else{
				if(@{$links{$ctg}{$side1}} > 1){
					$m0_stat{$ctg}{$side1}=1;
					process_side_ctg($ctg,$side1);
				}
				elsif(@{$links{$ctg}{$side1}} == 1){
					my $ctg1=$links{$ctg}{$side1}->[0];
					my $ctg1o=substr($ctg1,-1,1,"");
					#print "test2: ".join(",",@{$links{"ctg3"}{'right'}})."\n";
					#print "ctg: $ctg ori: $side1 ; linked: ".join(",",@{$links{"ctg3"}{'right'}})."\n";#test
					if(@{$links{$ctg1}{$ori2side{$side1}{$ctg1o}}} > 1){
						$m0_stat{$ctg}{$side1}=-1;
						process_side_ctg($ctg,$side1);
					}
					else{
						$m0_stat{$ctg}{$side1}=1;
						process_side_ctg($ctg,$side1);
					}
				}
				else{
					$m0_stat{$ctg}{$side1}=1;
				}
			}
		}
	}
	
	#m0操作并输出
	foreach my $pro_ctg(@ctgs_sort){
		print "$pro_ctg : left : ".$m0_stat{$pro_ctg}{'left'}."\n";
		print "$pro_ctg : right : ".$m0_stat{$pro_ctg}{'right'}."\n";
		if($m0_stat{$pro_ctg}{'left'} == -1){
			substr(${$gfaS}{$pro_ctg}{'seq'},0,$mlen{$pro_ctg}{'left'},'');
		}
		if($m0_stat{$pro_ctg}{'right'} == -1){
			substr(${$gfaS}{$pro_ctg}{'seq'},0-$mlen{$pro_ctg}{'right'},$mlen{$pro_ctg}{'right'},'');
		}
	}
	
	$out_p=$in if length($out_p) == 0;
	open outgfa,">${out_p}.m0.gfa";
	my @gfaH=grep {/^H/} @gfa_content;
	print outgfa join("\n",@gfaH)."\n" if @gfaH > 0;
	
	foreach (@ctgs_sort){
		print outgfa "S\t".$_."\t".${$gfaS}{$_}{'seq'}."\tDP:f:".${$gfaS}{$_}{'dep'}."\n";
	}
	
	foreach (@{$gfaL}){
		print outgfa "L\t".join("\t",@{$_}[0..3])."\t0M\n";
	}
	
	foreach (@gfa_content){
		print outgfa $_."\n" if /^[^HSL]/;
	}
	close outgfa;
}

sub process_side_ctg{#ctg,side
	my ($pctg_n,$pctg_side)=@_;

	foreach my $lk_ctg(@{$links{$pctg_n}{$pctg_side}}){
		my $ctgn=$lk_ctg;
		my $ctgo=substr($ctgn,-1,1,"");
		if(exists $m0_stat{$ctgn}{$ori2side{$pctg_side}{$ctgo}} && $m0_stat{$ctgn}{$ori2side{$pctg_side}{$ctgo}} ne $m0_stat{$pctg_n}{$pctg_side}){
			;
		}
		elsif(not exists $m0_stat{$ctgn}{$ori2side{$pctg_side}{$ctgo}}){
			$m0_stat{$ctgn}{$ori2side{$pctg_side}{$ctgo}} = $m0_stat{$pctg_n}{$pctg_side} * -1;
		}
		else{
			die"Error: Maybe some mistakes existed in the gfa-links or the gfa-links is too complex to be further processed!\n";
		}
	}
}

1;