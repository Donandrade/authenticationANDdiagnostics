#! /user/bin/perl -w

my $input= $ARGV[0];

open (ARQUIVO_ENTRADA,"<$input");
@linhas= <ARQUIVO_ENTRADA>;
close ARQUIVO_ENTRADA;

for ($i=0; $i <= $#linhas; $i++){
	chomp $linhas[$i];
	@colunas = split (/\s+/, $linhas[$i]);
	$numero_genes = (@colunas - 1);
	$tcca = 0;
	$tcar = 0;
	$cca = 0;
	$car = 0;
	for($j=1; $j <= $#colunas; $j++){
		@elementos = split (/\|/, $colunas[$j]);
		if ($elementos[0] eq "ccan"){
			$cca = $cca + 1;
			$tcca = 1;
		}
		elsif ($elementos[0] eq "cara"){
			$car = $car + 1;
			$tcar = 1;
	}
	$numero_genes= $cca + $car;
	$numero_taxons = $tcca + $tcar;
	}
#$numero_taxons = $tcca + $tcar;
print "$colunas[0]\ttaxons:\t$numero_taxons\tgenes:\t$numero_genes\t@colunas";
print "\n";
}
