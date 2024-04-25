#Calculating bond length, bond angle and torsion angle
use Math::Trig; #Trigonometric Module of Perl used to find Inverse
open(F,"2msa.pdb") or die "File not found\n";
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
    if($a=~/ATOM\s+1\s+\w+\s+.{5}\s+(\d+)/) #Pattern to find what is the first atom number
    {
    $first_atom_no=$1;
    }
    if($a=~/ATOM\s+\d+\s+\w+\s+.{5}\s+(\d+)/) #Pattern to find what is the last atom number
    {
    $last_atom_no=$1;
    }
    for($i=$first_atom_no;$i<=$last_atom_no;$i++) #Loop to extract the N, CA and C coordinates
    {
        if($a=~/ATOM\s+\d+\s+(..)\s+(.{5})\s+($i)\s+(.{7})\s+(.{7})\s+(.{7})/)
        {
            if($1 eq 'N ')
            {
                push @aa_name,$2; #To find the amino acid name where atom N is present  
                push @aa_no,$3; #To find the amino acid number where atom N is present 
                push @xn,$4; #To find the x coordinate of atom N
                push @yn,$5; #To find the y coordinate of atom N
                push @zn,$6; #To find the z coordinate of atom N
            }
            elsif($1 eq 'CA')
            {
                push @aa_name,$2; #To find the amino acid name where atom CA is present
                push @aa_no,$3; #To find the amino acid number where atom CA is present
                push @xca,$4; #To find the x coordinate of atom CA
                push @yca,$5; #To find the y coordinate of atom CA
                push @zca,$6; #To find the z coordinate of atom N
            }
            elsif($1 eq 'C ')
            {
                push @aa_name,$2; #To find the amino acid name where atom C is present
                push @aa_no,$3; #To find the amino acid number where atom C is present
                push @xc,$4; #To find the x coordinate of atom C
                push @yc,$5; #To find the y coordinate of atom C
                push @zc,$6; #To find the z coordinate of atom C
            }
            else
            {
                break;
            }
        }
    }
}
$length=@xn; #Finding the length of any of the coordinates to run a loop that can find bond distance, bond angle and torsion angle
close F;
for($j=0;$j<$length;$j++) #Loop for finding Bond Lengths
{
    $bond_length_n_ca=((($xn[$j]-$xca[$j])**2)+(($yn[$j]-$yca[$j])**2)+(($zn[$j]-$zca[$j])**2))**0.5; #To find the bond length between N and CA
    $bond_length_ca_c=((($xca[$j]-$xc[$j])**2)+(($yca[$j]-$yc[$j])**2)+(($zca[$j]-$zc[$j])**2))**0.5; #To find te bond length between CA and C
    $bond_length_c_n=((($xc[$j]-$xn[$j+1])**2)+(($yc[$j]-$yn[$j+1])**2)+(($zc[$j]-$zn[$j+1])**2))**0.5; #To find the bond length between C and N
    push @bond_length_n_ca,$bond_length_n_ca; #Creating an array to store bond length between N and CA
    push @bond_length_ca_c,$bond_length_ca_c; #Creating an array to store bond length between CA and C
    push @bond_length_c_n,$bond_length_c_n; #Creating an array to store bond length between C and N
}
for($j=0;$j<$length;$j++) #Loop for finding Bond Angles
{
    $bond_angle_n_ca_c=((($xn[$j]-$xca[$j])/$bond_length_n_ca[$j])*(($xc[$j]-$xca[$j])/$bond_length_ca_c[$j]))+((($yn[$j]-$yca[$j])/$bond_length_n_ca[$j])*(($yc[$j]-$yca[$j])/$bond_length_ca_c[$j]))+((($zn[$j]-$zca[$j])/$bond_length_n_ca[$j])*(($zc[$j]-$zca[$j])/$bond_length_ca_c[$j])); #To find the bond angle between N, CA and C by finding the dot product between direction cosines N -> CA and CA -> C
    if(($j+1) != $length) #Condition not to calculate the bond angle between CA, C, N and C, N, CA as the last bond angle in any protein sequence is between N, CA and C
    {
        $bond_angle_ca_c_n=((($xc[$j]-$xca[$j])/$bond_length_ca_c[$j])*(($xc[$j]-$xn[$j+1])/$bond_length_c_n[$j]))+((($yc[$j]-$yca[$j])/$bond_length_ca_c[$j])*(($yc[$j]-$yn[$j+1])/$bond_length_c_n[$j]))+((($zc[$j]-$zca[$j])/$bond_length_ca_c[$j])*(($zc[$j]-$zn[$j+1])/$bond_length_c_n[$j])); #To find the bond angle between CA, C and N by finding the dot product between direction cosines CA -> C and C -> N
        $bond_angle_c_n_ca=((($xc[$j]-$xn[$j+1])/$bond_length_c_n[$j])*(($xca[$j+1]-$xn[$j+1])/$bond_length_n_ca[$k+1]))+((($yc[$j]-$yn[$j+1])/$bond_length_c_n[$j])*(($yca[$j+1]-$yn[$j+1])/$bond_length_n_ca[$k+1]))+((($zc[$j]-$zn[$j+1])/$bond_length_c_n[$j])*(($zca[$j+1]-$zn[$j+1])/$bond_length_n_ca[$k+1])); #To find the bond angle between C, N and CA by finding the dot product between direction cosines C -> N and N -> CA
    }
    push @bond_angle_n_ca_c,$bond_angle_n_ca_c; #Creating an array to store bond angle between N, CA and C
    push @bond_angle_ca_c_n,$bond_angle_ca_c_n; #Creating an array to store bond angle between CA, C and N
    push @bond_angle_c_n_ca,$bond_angle_c_n_ca; #Creating an array to store bond angle between C, N and CA
}
for($j=0;$j<($length-1);$j++) #Loop for finding Torsion Angles
{
    $bond_vector_x_n_ca=$xca[$j]-$xn[$j]; #To find the difference between the x coordinates of the vector N and CA
    $bond_vector_y_n_ca=$yca[$j]-$yn[$j]; #To find the difference between the y coordinates of the vector N and CA
    $bond_vector_z_n_ca=$zca[$j]-$zn[$j]; #To find the difference between the z coordinates of the vector N and CA
    $bond_vector_x_ca_c=$xc[$j]-$xca[$j]; #To find the difference between the x coordinates of the vector CA and C
    $bond_vector_y_ca_c=$yc[$j]-$yca[$j]; #To find the difference between the y coordinates of the vector CA and C
    $bond_vector_z_ca_c=$zc[$j]-$zca[$j]; #To find the difference between the z coordinates of the vector CA and C
    $bond_vector_x_c_n=$xn[$j+1]-$xc[$j]; #To find the difference between the x coordinates of the vector C and N(from next amino acid)
    $bond_vector_y_c_n=$yn[$j+1]-$yc[$j]; #To find the difference between the y coordinates of the vector C and N(from next amino acid)
    $bond_vector_z_c_n=$zn[$j+1]-$zc[$j]; #To find the difference between the z coordinates of the vector C and N(from next amino acid)
    $bond_vector_x_n_ca_next=$xca[$j+1]-$xn[$j+1]; #To find the difference between the x coordinates of the vector N(from next amino acid) and CA(from next amino acid)
    $bond_vector_y_n_ca_next=$yca[$j+1]-$yn[$j+1]; #To find the difference between the y coordinates of the vector N(from next amino acid) and CA(from next amino acid)
    $bond_vector_z_n_ca_next=$zca[$j+1]-$zn[$j+1]; #To find the difference between the z coordinates of the vector N(from next amino acid) and CA(from next amino acid)
    $bond_vector_x_ca_c_next=$xc[$j+1]-$xca[$j+1]; #To find the difference between the x coordinates of the vector CA(from next amino acid) and C(from next amino acid)
    $bond_vector_y_ca_c_next=$yc[$j+1]-$yca[$j+1]; #To find the difference between the y coordinates of the vector CA(from next amino acid) and C(from next amino acid)
    $bond_vector_z_ca_c_next=$zc[$j+1]-$zca[$j+1]; #To find the difference between the z coordinates of the vector CA(from next amino acid) and C(from next amino acid)
    $cross_pdt_n_ca_c_i=($bond_vector_y_n_ca*$bond_vector_z_ca_c)-($bond_vector_z_n_ca*$bond_vector_y_ca_c); #To find the x coordinates of the normal from the plane formed by atoms N, CA and C by calculating the vector product
    $cross_pdt_n_ca_c_j=($bond_vector_z_n_ca*$bond_vector_x_ca_c)-($bond_vector_x_n_ca*$bond_vector_z_ca_c); #To find the y coordinates of the normal from the plane formed by atoms N, CA and C by calculating the vector product
    $cross_pdt_n_ca_c_k=($bond_vector_x_n_ca*$bond_vector_y_ca_c)-($bond_vector_y_n_ca*$bond_vector_x_ca_c); #To find the z coordinates of the normal from the plane formed by atoms N, CA and C by calculating the vector product
    $cross_pdt_ca_c_n_i=($bond_vector_y_ca_c*$bond_vector_z_c_n)-($bond_vector_z_ca_c*$bond_vector_y_c_n); #To find the x coordinates of the normal from the plane formed by atoms CA, C and N(from next amino acid) by calculating the vector product
    $cross_pdt_ca_c_n_j=($bond_vector_z_ca_c*$bond_vector_x_c_n)-($bond_vector_x_ca_c*$bond_vector_z_c_n); #To find the y coordinates of the normal from the plane formed by atoms CA, C and N(from next amino acid) by calculating the vector product
    $cross_pdt_ca_c_n_k=($bond_vector_x_ca_c*$bond_vector_y_c_n)-($bond_vector_y_ca_c*$bond_vector_x_c_n); #To find the z coordinates of the normal from the plane formed by atoms CA, C and N(from next amino acid) by calculating the vector product
    $cross_pdt_c_n_ca_i=($bond_vector_y_c_n*$bond_vector_z_n_ca_next)-($bond_vector_z_c_n*$bond_vector_y_n_ca_next); #To find the x coordinates of the normal from the plane formed by atoms C, N(from next amino acid) and CA(from next amino acid) by calculating the vector product
    $cross_pdt_c_n_ca_j=($bond_vector_z_c_n*$bond_vector_x_n_ca_next)-($bond_vector_x_c_n*$bond_vector_z_n_ca_next); #To find the y coordinates of the normal from the plane formed by atoms C, N(from next amino acid) and CA(from next amino acid) by calculating the vector product
    $cross_pdt_c_n_ca_k=($bond_vector_x_c_n*$bond_vector_y_n_ca_next)-($bond_vector_y_c_n*$bond_vector_x_n_ca_next); #To find the z coordinates of the normal from the plane formed by atoms C, N(from next amino acid) and CA(from next amino acid) by calculating the vector product
    $cross_pdt_n_ca_c_next_i=($bond_vector_y_n_ca_next*$bond_vector_z_ca_c_next)-($bond_vector_z_n_ca_next*$bond_vector_y_ca_c_next); #To find the x coordinates of the normal from the plane formed by atoms N(from next amino acid), CA(from next amino acid) and C(from next amino acid) by calculating the vector product
    $cross_pdt_n_ca_c_next_j=($bond_vector_z_n_ca_next*$bond_vector_x_ca_c_next)-($bond_vector_x_n_ca_next*$bond_vector_z_ca_c_next); #To find the y coordinates of the normal from the plane formed by atoms N(from next amino acid), CA(from next amino acid) and C(from next amino acid) by calculating the vector product
    $cross_pdt_n_ca_c_next_k=($bond_vector_x_n_ca_next*$bond_vector_y_ca_c_next)-($bond_vector_y_n_ca_next*$bond_vector_x_ca_c_next); #To find the z coordinates of the normal from the plane formed by atoms N(from next amino acid), CA(from next amino acid) and C(from next amino acid) by calculating the vector product
    $dot_pdt_n_ca_c_n=(($cross_pdt_n_ca_c_i*$cross_pdt_ca_c_n_i)+($cross_pdt_n_ca_c_j*$cross_pdt_ca_c_n_j)+($cross_pdt_n_ca_c_k*$cross_pdt_ca_c_n_k))/(((($cross_pdt_n_ca_c_i**2)+($cross_pdt_n_ca_c_j**2)+($cross_pdt_n_ca_c_k**2))**0.5)*((($cross_pdt_ca_c_n_i**2)+($cross_pdt_ca_c_n_j**2)+($cross_pdt_ca_c_n_k**2))**0.5)); #To find the angle between the two normals from the planes N, CA, C and CA, C, N(from next amino acid) by calculating the dot product
    $dot_pdt_ca_c_n_ca=(($cross_pdt_c_n_ca_i*$cross_pdt_ca_c_n_i)+($cross_pdt_c_n_ca_j*$cross_pdt_ca_c_n_j)+($cross_pdt_c_n_ca_k*$cross_pdt_ca_c_n_k))/(((($cross_pdt_c_n_ca_i**2)+($cross_pdt_c_n_ca_j**2)+($cross_pdt_c_n_ca_k**2))**0.5)*((($cross_pdt_ca_c_n_i**2)+($cross_pdt_ca_c_n_j**2)+($cross_pdt_ca_c_n_k**2))**0.5)); #To find the angle between the two normals from the planes CA, C, N(from next amino acid) and C, N(from next amino acid), CA(from next amino acid) by calculating the dot product
    $dot_pdt_c_n_ca_c=(($cross_pdt_c_n_ca_i*$cross_pdt_n_ca_c_next_i)+($cross_pdt_c_n_ca_j*$cross_pdt_n_ca_c_next_j)+($cross_pdt_c_n_ca_k*$cross_pdt_n_ca_c_next_k))/(((($cross_pdt_c_n_ca_i**2)+($cross_pdt_c_n_ca_j**2)+($cross_pdt_c_n_ca_k**2))**0.5)*((($cross_pdt_n_ca_c_next_i**2)+($cross_pdt_n_ca_c_next_j**2)+($cross_pdt_n_ca_c_next_k**2))**0.5)); #To find the angle between the two normals from the planes C, N(from next amino acid), CA(from next amino acid) and N(from next amino acid), CA(from next amino acid), C(from next amino acid) by calculating the dot product
    push @dot_pdt_n_ca_c_n,$dot_pdt_n_ca_c_n; #Creating an array to store the torsion angle between N, CA, C and N(from next amino acid)
    push @dot_pdt_ca_c_n_ca,$dot_pdt_ca_c_n_ca; #Creating an array to store the torsion angle between CA, C, N(from next amino acid) and CA(from next amino acid)
    push @dot_pdt_c_n_ca_c,$dot_pdt_c_n_ca_c; #Creating an array to store the torsion angle between C, N(from next amino acid), CA(from next amino acid) and C(from next amino acid)
}
print "Press 1 for calculating Bond Length, 2 for calculating Bond Angle, 3 for calculating Torsion Angle and any other number to exit\n"; #Choice to select whether to calculate Bond Length, Bond Angle or Torsion Angle
chomp($n=<>);
$count=0;
if($n == 1)
{
    for($k=0;$k<$length;$k++) #Loop to print Bond Lengths
    {
        print "The bond length between N".$aa_no[$count]." and CA".$aa_no[$count]." of $aa_name[$count] is $bond_length_n_ca[$k] Angstrom\n"; #Printing the bond length between N and CA
        print "The bond length between CA".$aa_no[$count]." and C".$aa_no[$count]." of $aa_name[$count] is $bond_length_ca_c[$k] Angstrom\n"; #Printing the bond length between CA and C
        if($k != ($length-1)) #Condition not to calculate the bond length of C and N at the end as the last bond length is between CA and C
        {
            print "The bond length between C".$aa_no[$count]." of $aa_name[$count] and N".$aa_no[$count+3]." of $aa_name[$count+3] is $bond_length_c_n[$k] Angstrom\n"; #Printing the bond length betwen C and N(from next amino acid)
        }
        $count+=3; #As atom names and atom numbers are stored 3 times, the next atom can be extracted by adding 3 to the previous count
    }
}
elsif($n == 2) 
{
    for($j=0;$j<$length;$j++) #Loop to print Bond Angles
    {
        print "The bond angle of atoms N".$aa_no[$count].", CA".$aa_no[$count]." and C".$aa_no[$count]." is: ".rad2deg(acos($bond_angle_n_ca_c[$j]))." degrees\n"; #Printing the bond length between N, CA and C
        if($j != ($length-1)) #Condition not to calculate the bond angle between CA, C, N(from next amino acid) and C, N(from next amino acid), CA(from next amino acid) at the end as the last bond angle is between N, CA, C
        {
            print "The bond angle of atoms CA".$aa_no[$count].", C".$aa_no[$count]." and N".$aa_no[$count+3]." is: ".rad2deg(acos($bond_angle_ca_c_n[$j]))." degrees\n"; #Printing the bond length between CA, C and N(from next amino acid)
            print "The bond angle of atoms C".$aa_no[$count].", N".$aa_no[$count+3]." and CA".$aa_no[$count+3]." is: ".rad2deg(acos($bond_angle_c_n_ca[$j]))." degrees\n"; #Printing the bond length between C, N(from next amino acid) and CA(from next amino acid)
        }
        $count+=3; #As atom names and atom numbers are stored 3 times, the next atom can be extracted by adding 3 to the previous count
    }
}
elsif($n == 3)
{
    for($k=0;$k<($length-1);$k++) #Loop to print Torsion Angles and the loop runs 2 less than the length as torsion angles are less than the total number of backbone atoms
    {
        print "The torsion angle of atoms N".$aa_no[$count].", CA".$aa_no[$count].", C".$aa_no[$count]." and N".$aa_no[$count+3]." is :".rad2deg(acos($dot_pdt_n_ca_c_n[$k]))." degrees\n"; #Printing the torsion angle between N, CA, C and N(from next amino acid)
        print "The torsion angle of atoms CA".$aa_no[$count].", C".$aa_no[$count].", N".$aa_no[$count+3]." and CA".$aa_no[$count+3]." is :".rad2deg(acos($dot_pdt_ca_c_n_ca[$k]))." degrees\n"; #Printing the torsion angle between CA, C, N(from next amino acid) and CA(from next amino acid)
        print "The torsion angle of atoms C".$aa_no[$count].", N".$aa_no[$count+3].", CA".$aa_no[$count+3]." and C".$aa_no[$count+3]." is :".rad2deg(acos($dot_pdt_c_n_ca_c[$k]))." degrees\n"; #Printing the torsion angle between C, N(from next amino acid), CA(from next amino acid) and C(from next amino acid)
        $count+=3; #As atom names and atom numbers are stored 3 times, the next atom can be extracted by adding 3 to the previous count
    }
}
else
{
    print "Press the correct number\n"; #If wrong choice is chosen
}
#close G;