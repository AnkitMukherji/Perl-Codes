#To find which atoms are bonded to which atoms
open(F,"diala.pdb") or die "File not found";
while(chomp($l=<F>))
{
if($l=~/(ATOM\s+(\d+).*)/) #Pattern to extract information out of the PDB file
{
$total_atom_no=$2; #Finding the last atom number
push @file,$1; #Storing all the lines containing ATOM properties into an array
}
}
foreach $i(@file)
{
	if($i=~/ATOM\s+(\d+)\s+(\S+)\s+(.{5})\s+(\d+)\s+(.{7})\s+(.{7})\s+(.{7})\s+\S+\s+\S+\s+(\S+)/) #Pattern to extract different ATOM properties
	{
	$atom_no=$1; #To store the atom numbers
	$atom_name=$2; #To store the atom names
	$aa_name=$3; #To store the amino acid names
	$aa_no=$4; #To store the amino acid numbers
	$x=$5; #To store the x coordinates of the respective atoms
	$y=$6; #To store the y coordinates of the respective atoms
	$z=$7; #To store the z coordinates of the respective atoms
	$atom=$8; #To store the atom names in terms of N, C, O, H or S
		foreach $j($atom_no+1..$total_atom_no) #Loop to check the distance of the atom stored in $atom_no from the consecutive atoms 
		{
			if($file[$j-1]=~/ATOM\s+($j)\s+(\S+)\s+.{5}\s+($aa_no)\s+(.{7})\s+(.{7})\s+(.{7})\s+\S+\s+\S+\s+(\S+)/) #Array $file[$j-1] used to capture the atom properties from the next position of the query atom 
			{
			$d=((($x-$4)**2)+(($y-$5)**2)+(($z-$6)**2))**0.5; #Formula to find the distance
				if($d > 1.33 and $d < 1.82) #Distance range for single bond
				{
					if($atom eq 'H' and $7 eq 'H') #Condition not to consider the bond distance between H atoms
					{
					break;
					}
					else
					{
					print "$atom_no is connected to $1 by single bond\n";
					}
				}
				elsif($d > 1.18 and $d < 1.25) #Distance range for double bond
				{
					if($atom eq 'H' and $7 eq 'H') #Condition not to consider the bond distance between H atoms
					{
					break;
					}
					else
					{
					print "$atom_no is connected to $1 by double bond\n";
					}
				}
				elsif($d > 1.06 and $d < 1.10) #Distance range for single bond, considered differently as the range is different
				{
					if($atom eq 'H' and $7 eq 'H')
					{
					break;
					}
					else
					{
					print "$atom_no is connected to $1 by single bond\n";
					}
				}
				else
				{
				break;
				}
			}
			elsif($file[$j-1]=~/ATOM\s+($j)\s+(\S+)\s+.{5}\s+(\d+)\s+(.{7})\s+(.{7})\s+(.{7})\s+\S+\s+\S+\s+(\S+)/) #Condition to compare atoms involved in peptide bond
			{
			$d=((($x-$4)**2)+(($y-$5)**2)+(($z-$6)**2))**0.5;
				if($d > 1.32 and $d < 1.33) #Distance range for peptide bond
				{
				print "$atom_no is connected to $1 by peptide bond\n";
				}
			}
			else
			{
			break;
			}
		}
	}
}