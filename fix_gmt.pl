#!/usr/bin/perl -w

#original GMT file from GSEA website, this iteration of this script does not need to use it, since it is used by the previous script (annotate Nick's script name)
#open (GMT,"/home/rvidana/WestbrookLab/scripts/screen.gsea/datasets/c2.cp.v5.0.symbols.gmt");
#@gmt=<GMT>; close GMT;

#need to work on folder allocation , for now it is a manual operation
$file=$ARGV[0];
$folder="/home/rvidana/WestbrookLab/Bioinformatics/sarah/gsea/";
open (GMT,"$folder$file");
@gmt=<GMT>; close GMT;

#list of valid genes one column, can become an argument at a later time
open (POOL,"/home/rvidana/WestbrookLab/scripts/screen.gsea/datasets/shRNA.Index.V3.genes.txt");
@focus_pool=<POOL>; close POOL;

$length_file=length($file)-4;
$no_extension=substr($file,0,$length_file);
$newfile="$no_extension.pool.gmt";

#new files generated
open (NEWGMT,">$folder$newfile");

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
