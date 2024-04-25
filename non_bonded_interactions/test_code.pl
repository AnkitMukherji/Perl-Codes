#test code
#Backbone H bond
open(F,"101m.pdb") or die "File not found";
while($l=<F>)
{
	for($i=0;$i<=2;$i++)
	{
		if($l=~/ATOM\s+\d+\s+([ON])\s+(\w+)\s+(\w+)\s+$i\s+(.{7})\s+(.{7})\s+(.{7})/)
			{
				if($1 eq 'O')
				{
					push @amino_acid,$2;
					push @chain,$3;
					push @xo,$4;
					push @yo,$5;
					push @zo,$6;
					
				}
				else
				{
					push @amino_acid,$2;
					push @chain,$3;
					push @xn,$4;
					push @yn,$5;
					push @zn,$6;
				}
			}
	}
}
close F;
#print "@amino_acid\n\n@chain\n\n@xo\n\n@yo\n\n@zo\n\n@xn\n\n@yn\n\n@zn\n";
$length_xo=@xo;
$length_xn=@xn;
for($j=0;$j<$length_xo;$j++)
{	
	for($k=0;$k<$length_xn;$k++)
	{
	$d_xo_xn=$xo[$j]-$xn[$k];
	$d_yo_yn=$xo[$j]-$xn[$k];
	$d_zo_zn=$xo[$j]-$xn[$k];
	push @d_xo_xn,($d_xo_xn)**2;
	push @d_yo_yn,($d_yo_yn)**2;
	push @d_zo_zn,($d_zo_zn)**2;
	}
}
#print "@d_xo_xn\n\n@d_yo_yn\n\n@d_zo_zn\n";
$length_xo_xn=@d_xo_xn;
for($z=0;$z<$length_xo_xn;$z++)
{
$diff_o_n=($d_xo_xn[$z]+$d_yo_yn[$z]+$d_zo_zn[$z])**0.5;
$c++;
		if($diff_o_n >= 2.7 and $diff_o_n <= 3.3)
		{
		print "Chain $chain[$c-1]: There is a probability of formation of Hydrogen bond as the distance between donor atom $amino_acid[$c-1] and acceptor atom $amino_acid[$c+1] is $diff_o_n\n";
		}
	}
	