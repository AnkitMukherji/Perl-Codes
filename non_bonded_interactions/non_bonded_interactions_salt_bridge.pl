open(F,"1acd.pdb") or die "File not found\n";
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
        if($a=~/ATOM\s+\d+\s+(N |O )/) #Pattern to exclude atoms containing main chain N and O as salt bridge is formed by charged atoms of side chain
        {
            break;
        }
        else
        {
            push @atom,$a; #Pushing the lines containing ATOMS into an array
        }
    }
}
if($models > 0)
{
    print "There are $models Models in this PDB file so repetitive chain breaks between same atom numbers and same amino acids indicate they come from different models\n";
}
foreach $i(@atom) #Loop to extract the N and O coordinates
{
    if($i=~/ATOM\s+\d+ {0,}(.{1,4}) {0,}(.{6}) {0,}(.{1,3})\s+(.{7})\s+(.{7})\s+(.{7})\s+\S+\s+\S+\s+(.)/)
    {
        if($7 eq 'O') #Condition to extract the properties of the atom O in every amino acid
        {
            push @atom_name_o,$1; #To find the atom names of different O atoms
            push @aa_name_o,$2; #To find the amino acid name where atom O is present   
            push @aa_no_o,$3; #To find the amino acid number where atom O is present
            push @xo,$4; #To find the x coordinate of atom O
            push @yo,$5; #To find the y coordinate of atom O
            push @zo,$6; #To find the z coordinate of atom O
        }
        elsif($7 eq 'N') #Condition to extract the properties of the atom N in every amino acid
        {
            push @atom_name_n,$1; #To find the atom names of different N atoms
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
for($j=0;$j<$length_o_atom;$j++)
{
	for($k=0;$k<$length_n_atom;$k++)
	{
		$h_bond_dist=((($xo[$j]-$xn[$k])**2)+(($yo[$j]-$yn[$k])**2)+(($zo[$j]-$zn[$k])**2))**0.5; #Distance calculated by taking into account one of the O atom coordinates and then comparing it with all the coordinates of N atom at each round
		if($h_bond_dist < 4) #Printing the distance if it is within the range of formation of salt bridge i.e. 4 Angstrom
		{
            if($aa_no_n[$k] != $aa_no_o[$j]) #Condition given not to consider the salt bridge formation within the same amino acid
            {
                print "There exists a Salt Bridge between ".$atom_name_n[$k]." of amino acid ".$aa_no_n[$k]." ".$aa_name_n[$k]." and ".$atom_name_o[$j]." of amino acid ".$aa_no_o[$j]." ".$aa_name_o[$j]." as the distance between them is ".$h_bond_dist." Angstrom\n";
            }
		}
	}
} 
