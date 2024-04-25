open(F,"1abe.pdb") or die "File not found\n";
while($a=<F>)
{
    if($a=~/HEADER.{56}(.{4})/) #Pattern to extract PDB ID
    {
    print "The PDB ID is: $1\n";
    }
    if($a=~/TITLE\s+(.*)/) #Pattern to extract Title
    {
    print "The Title is: $1\n";
    }
    if($a=~/MODEL\s+(\d+)/) #Pattern to extract the total number of models
    {
        $models=$1;
    }
    if($a=~/(^ATOM.*)/) #Pattern to extract all the lines containing ATOMS
    {
    push @atom,$a;
    }
}
if($models > 0)
{
    print "There are $models Models in this PDB file so repetitive chain breaks between same atom numbers and same amino acids indicate they come from different models\n";
}
foreach $i(@atom) #Loop to extract the N and O coordinates
{
    if($i=~/ATOM\s+\d+ {0,}(.{1,4}) {0,}(.{6}) {0,}(.{1,3})\s+(.{7})\s+(.{7})\s+(.{7})/)
    {
        if($1 eq 'O   ') #Condition to extract the properties of the atom O in every amino acid
        {
            push @aa_name_o,$2; #To find the amino acid name where atom O is present   
            push @aa_no_o,$3; #To find the amino acid number where atom O is present
            push @xo,$4; #To find the x coordinate of atom O
            push @yo,$5; #To find the y coordinate of atom O
            push @zo,$6; #To find the z coordinate of atom O
        }
        elsif($1 eq 'N   ') #Condition to extract the properties of the atom N in every amino acid
        {
            push @aa_name_n,$2; #To find the amino acid name where atom N is present
            push @aa_no_n,$3; #To find the amino acid number where atom N is present
            push @xn,$4; #To find the x coordinate of atom N
            push @yn,$5; #To find the y coordinate of atom N
            push @zn,$6; #To find the z coordinate of atom N
        }
        else
        {
            break;
        }
    }
}
$length_o_atom=@xo; #Counting the length of any of the coordinates of O atom to run the outer loop 
$length_n_atom=@xn; #Counting the length of any of the coordinates of N atom to run the inner loop 
for($j=0;$j<$length_o_atom;$j++) #Outer loop
{
	for($k=$j+1;$k<$length_n_atom;$k++) #Inner loop('$j+1' is taken not to consider the hydrogen bond between N and O of the same amino acid)
	{
		$h_bond_dist=((($xo[$j]-$xn[$k])**2)+(($yo[$j]-$yn[$k])**2)+(($zo[$j]-$zn[$k])**2))**0.5; #Distance calculated by taking into account one of the O atom coordinates and then comparing it with all the coordinates of N atom at each round
		if($h_bond_dist < 3.6) #Printing the distance if it is within the range of formation of H bond i.e. 3.6 Angstrom
		{
			print "There exists a Hydrogen Bond between N".$aa_no_n[$k]." of amino acid ".$aa_name_n[$k]." and O".$aa_no_o[$j]." of amino acid ".$aa_name_o[$j]." as the distance between them is ".$h_bond_dist." Angstrom\n";
		}
	}
} 
