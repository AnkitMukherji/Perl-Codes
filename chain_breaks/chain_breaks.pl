open(F,"8hca.pdb") or die "File not found\n";
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
foreach $i(@atom) #Loop to extract the N and C coordinates
{
    if($i=~/ATOM\s+(\d+) {0,}(.{1,4}) {0,}(.{6}) {0,}(.{1,3})\s+(.{7})\s+(.{7})\s+(.{7})/)
    {
        if($2 eq 'N   ') #Condition to extract the properties of the atom N in every amino acid
        {
            push @aa_name,$3; #To find the amino acid name where atom N is present  
            push @aa_no,$4; #To find the amino acid number where atom N is present 
            push @xn,$5; #To find the x coordinate of atom N
            push @yn,$6; #To find the y coordinate of atom N
            push @zn,$7; #To find the z coordinate of atom N
        }
        elsif($2 eq 'C   ') #Condition to extract the properties of the atom C in every amino acid
        {
            push @aa_name,$3; #To find the amino acid name where atom C is present
            push @aa_no,$4; #To find the amino acid number where atom C is present
            push @xc,$5; #To find the x coordinate of atom C
            push @yc,$6; #To find the y coordinate of atom C
            push @zc,$7; #To find the z coordinate of atom C
        }
        else
        {
            break;
        }
    }
}  
$length_coordinate=@xc; #Finding the length of any of the coordinates to run a loop that can find bond distance which will be further required for calculation of chain breaks
close F;
for($j=0;$j<$length_coordinate;$j++) #Loop for finding Bond Lengths
{
    if($j != 0) #Condition to skip the first N atom and count from the second N atom
    {
    $bond_length_c_n=((($xc[$j-1]-$xn[$j])**2)+(($yc[$j-1]-$yn[$j])**2)+(($zc[$j-1]-$zn[$j])**2))**0.5; #To find the bond length between C(from first amino acid) and N(from next amino acid)
    push @bond_length_c_n,$bond_length_c_n; #Creating an array to store bond length between C(from first amino acid) and N(from next amino acid)
    }
    else
    {
        break;
    }
}
$length_bond=@bond_length_c_n;
for($k=0;$k<$length_bond;$k++) #Loop for printing Bond Lengths
{
    if($bond_length_c_n[$k] > 1.35) #Condition for finding chain break
    {
        print "There is a chain break between atoms C".$aa_no[$count]." of amino acid ".$aa_name[$count]." and N".$aa_no[$count+2]." of amino acid ".$aa_name[$count+2]." as the distance between them is ".$bond_length_c_n[$k]."\n";
    }
    else
    {
        break;
    }
    $count+=2; #Count +2 is considered as two atoms are considered
}