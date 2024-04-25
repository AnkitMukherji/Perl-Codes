#Finding average quality score and %GC content
open(F,"fastq_code.fa") or die "File not found"; #Opening a fasta file to parse the average genome sequence quality and the per base GC content
open(G,">file.txt");
@fastqc=<F>; #storing the file in an array
close F;
$length=@fastqc; #Finding the length of the array
%average_quality=('@',0,'A',1,'B',2,'C',3,'D',4,'E',5,'F',6,'G',7,'H',8,'I',9,'J',10,'K',11,'L',12,'M',13,'N',14,'O',15,'P',16,'Q',17,'R',18,'S',19,'T',20,'U',21,'V',22,'W',23,'X',24,'Y',25,'Z',26,'[',27,'\/',28,']',29,'^',30,'_',31,'`',32,'a',33,'b',34,'c',35,'d',36,'e',37,'f',38,'g',39,'h',40,'i',41,'j',42); #Creating a dictionary using hash for the quality of each base
for($i=0;$i<$length;$i++) #Running a loop upto the length of @fastq to extract the sequence and the quality of each nucleotide
{
    if($fastqc[$i]=~/^\@HWUSI.*/) #Pattern to extract the sequence
    {
        $seq.=$fastqc[($i+1)]; #Appending the second line from each read containing 4 lines
        if($i == 0)
        {
            $seq_length=length($seq); #Finding the length of each read
        }
    }
    if($fastqc[$i+2]=~/^\+HWUSI.*/) #Pattern to extract the quality
    {
        $quality.=$fastqc[($i+3)]; #Appending the fourth line from each read containing 4 lines
    }
}
print G"The length of each read is: ",$seq_length-1,"\n";
@seq=split//,$seq; #Splitting the sequence into an array
@quality=split//,$quality; #Splitting the quality into an array
$total_read_length=@seq; #Finding the total number of nucleotides in the seq variable including the \n value
$total_read=$total_read_length/$seq_length; #finding the total number of reads from the illumina data
for($j=0;$j<$total_read_length;$j++) #Running a loop upto the total number of nucleotides to find the average quality score of the sequence
{
        $sum+=$average_quality{$quality[$j]}; #Finding the sum of the qualities of the respective nucleotides and piping it to the dictionary to extract the respective quality scores
        $avg=$sum/($total_read_length-$total_read); #Finding the average quality score by dividing the total quality score by the total number of nucleotides excluding the \n
}
print G"The total number of reads are: $total_read\n";
print G"The average quality score of the sequence is: $avg\n";  
for($k=0;$k<($seq_length-1);$k++) #Running a loop upto the length of each read to run the loop horizontally across a single read
{
    for($l=$k;$l<($total_read_length-($seq_length-($k+1)));$l+=$seq_length) #Running a loop upto the last nucleotide position of each vertical lane to run the loop vertically considering a single nucleotide from each read
    {
        if($seq[$l]=~/[GC]/) #Pattern to check whether the considered nucleotide is either G or C or not
        {
            $gc++; #Variable to count the total GC content of each vertical lane
        }
    }
    print G"The GC percentage for bases at position ",$k+1," of the sequence is: ",($gc/$total_read)*100," %\n";
    $gc=0; #Reinitializing the value of $gc to 0 to start counting the GC content of a new vertical lane from the beginning
}
