# Calculate the number of bonded interactions(bonds, angles, torsions) and non-bonded interactions
use List::MoreUtils qw(uniq);
use Storable qw/freeze/;
open(F,"diala_mol.mol2") or die "File not found\n";
open(G,">Modeling.txt");
print G"The bonded atoms are:\nBond lengths:\n";
@file = <F>;
$length = @file;
foreach $i(0..$length - 1)
{
    if($file[$i]=~/^@<TRIPOS>BOND/)
    {
        while($file[$i + 1] != '@<TRIPOS>SUBSTRUCTURE')
        {
            if($file[$i + 1]=~/\d+\s+(\d+)\s+(\d+)/)
            {
                push @atom,"$1 $2";
                print G"$1 - $2\n";
                $bond_count++;
            }
            $i++;
        }
    }
    else
    {
        break;
    }
}
print G"Bond angles:\n";
$length_1 = @atom;
foreach $k(0..$length_1 - 1) # For the bond angles, the bonded atoms are used. For eg, 5 - 2 and 2 - 3 are bonded as 5 - 2 - 3. This is checked in four ways
{
    @angle = split/ /,$atom[$k];
    while(($k+1) <= ($length_1 - 1)) # Loop to compare one bonded set with all others
    {
        @check_1 = split/ /,$atom[$k + 1]; 
        if($angle[1] == $check_1[0]) # Forward bonded and forward bonded eg 5 - 2 and 2 - 3
        {
            push @angle_2,"$angle[0] $angle[1] $check_1[1]";
            print G"$angle[0] - $angle[1] - $check_1[1]\n";
            $angle_count++; # Counter for counting the total number of bond angles
        }
        elsif($angle[1] == $check_1[1]) # Forward bonded and reverse bonded eg 5 - 2 and 3 - 2
        {
            push @angle_2,"$angle[0] $angle[1] $check_1[0]";
            print G"$angle[0] - $angle[1] - $check_1[0]\n";
            $angle_count++;
        }
        elsif($angle[0] == $check_1[0]) # Reverse bonded and forward bonded eg 2 - 5 and 2 - 3
        {
            push @angle_2,"$angle[1] $angle[0] $check_1[1]";
            print G"$angle[1] - $angle[0] - $check_1[1]\n";
            $angle_count++;
        }
        elsif($angle[0] == $check_1[1]) # Reverse bonded and reverse bonded eg 2 - 5 and 3 - 2 
        {
            push @angle_2,"$angle[1] $angle[0] $check_1[0]";
            print G"$angle[1] - $angle[0] - $check_1[0]\n";
            $angle_count++;
        }
        $k++;
    }
    print "It's running.....Checking for bond angles\n";
}
print G"Bond torsions:\n";
$length_2 = @angle_2;
foreach $z(0..$length_2 - 1) # The bonded angles are checked with itself to check for torsion angles. This is the same as before but here 2 atoms are checked
{
    @torsion = split/ /,$angle_2[$z];
    while(($z+1) <= ($length_2 - 1)) # eg 5 - 2 - 3 and 2 - 3 - 9 has the torsion 5 - 2 - 3 - 9. This is also checked in four ways
    {
        @check_2 = split/ /,$angle_2[$z + 1];
        if(($torsion[1] == $check_2[0]) and ($torsion[2] == $check_2[1]))
        {
            push @torsion_2,"$torsion[0] $torsion[1] $check_2[1] $check_2[2]";
            print G"$torsion[0] - $torsion[1] - $check_2[1] - $check_2[2]\n";
            $torsion_count++; # Counter for torsion count
        }
        elsif(($torsion[1] == $check_2[2]) and ($torsion[2] == $check_2[1]))
        {
            push @torsion_2,"$torsion[0] $torsion[1] $check_2[1] $check_2[0]";
            print G"$torsion[0] - $torsion[1] - $check_2[1] - $check_2[0]\n";
            $torsion_count++; 
        }
        elsif(($torsion[0] == $check_2[1]) and ($torsion[1] == $check_2[0]))
        {
            push @torsion_2,"$torsion[2] $torsion[1] $check_2[1] $check_2[2]";
            print G"$torsion[2] - $torsion[1] - $check_2[1] - $check_2[2]\n";
            $torsion_count++; 
        }
        elsif(($torsion[0] == $check_2[1]) and ($torsion[1] == $check_2[2]))
        {
            push @torsion_2,"$torsion[2] $torsion[1] $check_2[1] $check_2[0]";
            print G"$torsion[2] - $torsion[1] - $check_2[1] - $check_2[0]\n";
            $torsion_count++; 
        }
        $z++;
    }
    print "It's running.....Checking for bond torsions\n";
}
print G"There are $bond_count bonds, $angle_count angles, $torsion_count torsions\n";
foreach $g(0..10000) # Just running approx max number of iterations to check for all the non-bonded atoms
{
    @longest = split/ /,$torsion_2[$g]; 
    $length_3 = @longest;
    for($a = ($length_3 - 1); $a >= 0; $a--) # This loop reverses the torsion angles
    {
        push @rev,$longest[$a];
    }
    foreach $h(0..$length_1 - 1) # Recursively using @torsion_2 to add new bonds to find all the non bonded interactions
    {   #There are four conditions, same as before. Every condition is checking end and start points for match. Those values are eliminated where: 5 - 2 - 3 - 7 and 7 - 3
        @bond = split/ /,$atom[$h];
        if(($longest[-1] == $bond[0]) and ($longest[-2] != $bond[1]) and ($longest[0] != $bond[1]))
        {
            push @torsion_2,"@longest $bond[1]";
        }
        elsif(($longest[-1] == $bond[1]) and ($longest[-2] != $bond[0]) and ($longest[0] != $bond[0]))
        {
            push @torsion_2,"@longest $bond[0]";
        }
        elsif(($longest[0] == $bond[0]) and ($longest[1] != $bond[1]) and ($longest[-1] != $bond[1]))
        {
            push @torsion_2,"@rev $bond[1]";
        }
        elsif(($longest[0] == $bond[1]) and ($longest[1] != $bond[0]) and ($longest[-1] != $bond[0]))
        {
            push @torsion_2,"@rev $bond[0]";
        }
    }
    @rev = ();
    print "It's running.....Checking for non-bonded interactions\n";
}
@unique_1 = uniq(@torsion_2); # All the redundant values are removed
$length_4 = @unique_1;
foreach $c(0..$length_4 - 1) #This loop checks for redundant elements which are present in reversed form and eliminates them
{
    @check_3 = split/ /,$unique_1[$c];
    foreach $d($c + 1..$length_4 - 1) #This loop checks the considered element with all the forward elements
    {
        @check_4 = split/ /,$unique_1[$d];
        $length_5 = @check_4;
        for($f = $length_5 - 1; $f >= 0; $f--) #This loop reverses the elements
        {
            push @reverse, "$check_4[$f]";
        }
        if(freeze(\@check_3) eq freeze(\@reverse))
        {
            push @store, "@check_4\n"; # Array @store stores all the set of bonded elements whose reverse is present
        }
        elsif(freeze(\@check_3) ne freeze(\@reverse))
        {
            foreach $o(@store) # This loop checks whether the encountered reverse element(that was encountered before as forward) does not get printed and counter increases
            {
                if($o eq "@check_3\n")
                {
                    $count++;
                }
            }
            if($count > 0)
            {
                break;
            }
            else
            {
                push @nBond, "@check_3"; #If the counter is 0, i.e. no reverse and forward elements are printed store the element to @nBond
            }
        }
        @reverse = ();
        $count = 0;
    }
    print "It's running.....Checking and removing redundant elements $c\n";
}
@unique_2 = uniq(@nBond); #Again remove the redundant elements
print G"The Non-Bonded atoms are:\n";
foreach $d(@unique_2) #To print the elements
{
    @split = split/ /,$d;
    foreach $e(@split)
    {
        if($e == $split[-1]) # Condition not to print '-' for the last atom in an element
        {
            print G"$e";
        }
        else
        {
            print G"$e - ";
        }
    }
    print G"\n";
    $non_bonded_count++;
}
print G"The total number of non-bonded interactions are: $non_bonded_count\n";
print "You are set to go. Check in the Modeling file\n";
