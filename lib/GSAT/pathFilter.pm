#updated on Jul 8, 2022
package pathFilter;

sub filterPath{
    my ($lines,$seqname,$seqlen,$min_end_length)=@_;
    my @bandage;

    foreach (@{$lines}){
        chomp;
        my @p_info=split(/\t/);
        my $this_path=$p_info[4];
        my $f_stat=-1;
        my $l_stat=-1;
        $this_path =~ s/^\([^\(\)]+\)|\([^\(\)]+\)$//g;
        
        while($this_path =~ /^\{/){
            $this_path =~ s/^\{[^,\}]+\},//;
            $f_stat=1;
        }
        while($this_path =~ /\}$/){
            $this_path =~ s/,\{[^,\}]+\}$//;
            $l_stat=1;
        }
  
        $this_path =~ s/\{\@[^\}]+\},*|\([^\)]+\),*|\[[^\]]+\],*//g;#去除距离等附加信息
        my @this_paths=split(/,/,$this_path);
        if($f_stat == -1){
            my $f_ctg=$this_paths[0];
            $f_ctg =~ s/[\+\-]$//;
            die "Error: wrong ctg names were detected: first ctg($f_ctg)\n" if length($f_ctg) ==0;
            my ($f_length)=map {${$seqlen}[$_]} (grep {${$seqname}[$_] eq $f_ctg} 0..$#{$seqlen});
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
            my ($l_length)=map {${$seqlen}[$_]} (grep {${$seqname}[$_] eq $l_ctg} 0..$#{$seqlen});
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

    return @bandage;
}

1;