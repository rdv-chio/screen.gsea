#!/usr/bin/perl -w
open (GMT,"/home/rvidana/WestbrookLab/Bioinformatics/tingting/c2.all.v4.0.symbols.gmt");
@gmt=<GMT>; close GMT;

#open (POOL,"/home/rvidana/WestbrookLab/Bioinformatics/datasets/shRNA/shRNA.Index.V3.genes.txt");
open (POOL,"/home/rvidana/WestbrookLab/Bioinformatics/tiffany/rna_seq_genes.txt");
@focus_pool=<POOL>; close POOL;

#open (NEWGMT,">/home/rvidana/WestbrookLab/Bioinformatics/datasets/GSEA/c2.all.v4.0.symbols.pool.4.2.gmt");
open (NEWGMT,">/home/rvidana/WestbrookLab/Bioinformatics/datasets/GSEA/c2.all.v4.0.symbols.ngs.gmt");

foreach $pathway(@gmt)
{
	chomp $pathway;	
	@one_pathway=split(/\t/,$pathway);
	@focus_pathway=();
	push(@focus_pathway,$one_pathway[0]);
	push(@focus_pathway,$one_pathway[1]);
	for($i=2;$i< scalar @one_pathway;$i++)
	{
		foreach $genes_screen(@focus_pool)
		{	
			chomp $genes_screen;
			chomp $one_pathway[$i];
			if($one_pathway[$i] eq $genes_screen)
			{push(@focus_pathway,$one_pathway[$i]);}
		}
	}
	if (scalar @focus_pathway>12)
	{
		$final_pathway=join("\t",@focus_pathway);
		chomp $final_pathway;
		print "$pathway\n";
		print "$final_pathway\n";
		print NEWGMT "$final_pathway\n";
	}
}
close NEWGMT;
