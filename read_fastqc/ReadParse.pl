#Calculating the average Quality and GC% per base of the reads

%ascii=('@'=>0,'A'=>1,'B'=> 2,'C'=> 3,'D'=> 4,'E'=> 5,'F'=> 6,'G'=> 7,'H'=> 8,'I'=> 9,'J'=> 10,'K'=> 11,'L'=> 12,'M'=>13,
'N'=> 14,'O'=>15,'P'=> 16,'Q'=> 17,'R'=> 18,'S'=> 19,'T'=> 20,'U'=> 21,'V'=> 22,'W'=>23,'X'=> 24,'Y'=> 25,'Z'=> 26,'['=> 27,
'\\'=> 28,']'=> 29,'^'=>30,'_'=>31,'`'=> 32,'a'=> 33,'b'=> 34,'c'=> 35, 'd'=>36,'e'=>37,'f'=> 38,'g'=>39,'h'=>40,'i'=> 41,'j'=>42);
#foreach $i(keys%ascii)
#{print"$i:$ascii{$i}";}
open (F,"fastq_code.fa") or die "File not found\n";
open (G,">result.txt");
while (chomp($line= <F>))
{push @wholefile,$line;
$l=@wholefile;}
foreach $a(0..$l-1)
{
if (@wholefile[$a] =~ /^\@HWUSI.*/i)
		{$sequence= $sequence.$wholefile[($a+1)];}
if (@wholefile[$a] =~/^\+HWUSI.*/i)
		{$basecall= $basecall.$wholefile[($a+1)];}
}
close (F);
#$len=length$sequence;
print"$l\n";
=head
#print"@wholefile\n\n";
#print G"$sequence\n";
#print G"$basecall\n";
@split_seq=split//,$sequence;
@split_base=split//,$basecall;
$length=@split_seq;
#print"@split_base";
#print"$length\n";
#print"@split_seq\n";
$column=$l/4;
$row=$length/$column;
#print"row=$row,col=$column";
for ($c=0;$c<$column;$c++)
	{$quality+= $ascii{$split_base[$c]};
	}
	$avg_quality = ($quality/$column);
#print"$quality";
print "The Average Quality of the reads is $avg_quality\n";
print G"The Average Quality of the reads is $avg_quality\n";
for ($d=0;$d<$row-1;$d++)
{
	for($e=0;$e<$column-1;$e++)
		{
			if(@split_seq[$e]=~/[G|C]/i)
			{#print"@split_seq[$e]\n";
				$count++;
			#print"$count\n";		
			$gc=($count/$column)*100;			
			}
	
	
		}
print "The GC% of bases in position $d is $gc\n";
						print G"The GC% of bases in position $d is $gc\n"; 

}
close (G);
