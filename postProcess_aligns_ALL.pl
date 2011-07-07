#!/usr/bin/perl -w
### 07.25.08
my ($file) = @ARGV;
die "\n USAGE: postProcess.pl <.aligns file>\n\n" unless (scalar(@ARGV)==1);

my %species1_pos;
my %species2_pos;

### now reformat the following file
my %aln_blocks;
open (OUT,">$file.postProcess.csv") || die "ERROR: Can't create $file.postProcess: $!\n";
print OUT 'block_id,pos_in_block,id_1,id_2,chr_1,start_1,end_1,chr_2,start_2,end_2' . "\n";

open (FILE,$file);
while (<FILE>) { chomp $_; 
   # ## Alignment 0: score=517.0 e_value=8.9e-18 N=13 Dmel2L&DsimNODE_1092_length_59670_cov_31.675951 plus
   #   0-  0:        2L.11267474.11267836    NODE_1092_length_59670_cov_31.675951.56945.57307          3e+02
   #   0-  1:        2L.11267967.11268056    NODE_1092_length_59670_cov_31.675951.57414.57503          1e+02
   #   0-  2:        2L.11268116.11268648    NODE_1092_length_59670_cov_31.675951.57560.58096          1e+03
   #   0-  3:        2L.11268656.11268689    NODE_1092_length_59670_cov_31.675951.58099.58132          6e+01
   next if ($_ =~ /^\#/) || (! $_);
   my ($block_id,$pos_in_block,$used_gene1,$used_gene2) = ($1,$2,$3,$4) if ($_ =~ /(\d+)\-\s*(\d+)\:\s+(\S+)\s+(\S+)/);
	my ($thal_id,$lyr_id) = ($used_gene1,$used_gene2);
   #($thal_id,$lyr_id) = ($used_gene2,$used_gene1) if (($used_gene1 =~ /^lyr\./) && ($used_gene2 =~ /^thal\./));

   ($thal_id,$lyr_id) = ($used_gene2,$used_gene1) if ($thal_id =~ /^NODE/); ### reverse thal with lyr
#next if ($thal_id =~ /^NODE/);

	my ($chr_1,$start_1,$end_1) = split(/\./,$thal_id);
	my ($chr_2_1,$chr_2_2,$start_2,$end_2) = split(/\./,$lyr_id);
   my $chr_2 = "$chr_2_1.$chr_2_2";
#die "start_2 $_\n" unless ($start_2);
#die "end_2   ($chr_2_1,$chr_2_2,$start_2,$end_2) $_\n" unless ($end_2);

	### block_id,pos_in_block,id_1,id_2,chr_1,start_1,end_1,chr_2,start_2,end_2
	print OUT "$block_id,$block_id.$pos_in_block,$thal_id,$lyr_id,$chr_1,$start_1,$end_1,$chr_2,$start_2,$end_2\n";

	$aln_blocks{$block_id}{$pos_in_block} = "$thal_id,$lyr_id"; ### keep track of the block's boundaries
	$species1_pos{$thal_id} = &getMid($start_1,$end_1);
	$species2_pos{$lyr_id}  = &getMid($start_2,$end_2);
} close FILE;
close OUT;


open (BOUND,">$file.postProcess.boundaries.csv") || die "ERROR: Can't create $file.postProcess.boundaries: $!\n";
print BOUND 'block_id,num_in_block,',
				'chr_1,id_1_start,id_1_end,start_1,end_1,',
				'chr_2,id_2_start,id_2_end,start_2,end_2',"\n";

foreach my $block_id (sort {$a<=>$b} keys %aln_blocks) {
	my @members = sort {$a<=>$b} keys %{$aln_blocks{$block_id}};
	my ($member1_gene_start,$member2_gene_start) = split(/,/,$aln_blocks{$block_id}{$members[0]});
	my ($member1_gene_end,$member2_gene_end) = split(/,/,$aln_blocks{$block_id}{$members[$#members]});
#print "$block_id: ",scalar(@members),": $members[0] - $members[$#members]\n";

	my ($member1_gene_start_pos,$member1_gene_end_pos) = ($species1_pos{$member1_gene_start},$species1_pos{$member1_gene_end});
	my ($member2_gene_start_pos,$member2_gene_end_pos) = ($species2_pos{$member2_gene_start},$species2_pos{$member2_gene_end});
#print "\tFIRST 1=$member1_gene_start - 2=$member2_gene_start ($member1_gene_start_pos,$member2_gene_start_pos)\n";
#print "\tLAST  1=$member1_gene_end - 2=$member2_gene_end ($member1_gene_end_pos,$member2_gene_end_pos)\n";

	my $chr_1 = $1 if ($member1_gene_start =~ /^(\w+)\./);
	my ($chr_2_1,$chr_2_2) = ($1,$2) if ($member2_gene_start =~ /^(\w+)\.(\w+)\./);
   my $chr_2 = "$chr_2_1.$chr_2_2";

	print BOUND "$block_id,",scalar(@members),',',#,$orient,",
		"$chr_1,$member1_gene_start,$member1_gene_end,$species1_pos{$member1_gene_start},$species1_pos{$member1_gene_end},",
		"$chr_2,$member2_gene_start,$member2_gene_end,$species2_pos{$member2_gene_start},$species2_pos{$member2_gene_end}\n";
}

exit;


sub getMid {
	my ($start,$end) = @_;
	return ($start+$end)/2;
}
