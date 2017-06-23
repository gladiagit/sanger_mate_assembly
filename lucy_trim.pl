#!/usr/bin/perl

use strict;
use warnings; 
use diagnostics;
#a program that attempts to trim the sequence regions identified by lucy as being low quality
my ($ffile,$qfile,$out)=&param_check(3,"perl lucy_trim.pl fastafile  quailty_file out_filename_template");
&main($ffile,$out,$qfile);

#checks for the correct number of transmitted parameters v3
sub param_check
{
   my $nparam=$_[0];
   my $msg=$_[1];          chomp($msg);
   my @arguments=@ARGV;    chomp(@arguments);
  
   if ( ($msg) && ($#arguments<($nparam-1)))
   {  print "------------------------------------------------------------------------------------------\n";
      print $msg."\n";
      print "------------------------------------------------------------------------------------------\n";
   }
   for(my $i=0;$i<$nparam;$i++)
   {
      if (!$arguments[$i])
         {print "please give the ".($i+1)."st paramater \n";
          $arguments[$i]=<STDIN>;  chomp(@arguments);}
      
   }
   return @arguments;
}

# function that takes a filename and reads in and returns a the data in the file as an array of lines 
sub file_reader
{
   my $filename=$_[0];  chomp($filename);
   open(INP,"<$filename")||die "could not open $filename for reading $!\n";
   my @data=<INP>;
   return @data;

}

# a simple function that goes through the sequences in a given file 
sub fasta_file_parsing
{
   use Bio::SeqIO;

   my $file=$_[0];         chomp($file);
   my %seqs;
   my $seqio_object = Bio::SeqIO->new(-file => $file, -format=>'fasta');
   while (my $seq_object = $seqio_object->next_seq)
   {
      if (!exists $seqs{$seq_object->display_id})
      {
         if ($seq_object->description)
            {$seqs{$seq_object->display_id}{'desc'}=$seq_object->description;}
        $seqs{$seq_object->display_id}{'seq'}=$seq_object->seq;
      }
      else
         {print $seq_object->display_id." sequence is already in the hash\n";}
   }
   return %seqs;
}
# file that parses a quality file
sub qual_file_parsing
{
   use Bio::SeqIO;
   my $sfile=$_[0];      chomp($sfile);   #the sequence file
   my %qual;
   my $seqio_object = Bio::SeqIO->new(-file => $sfile , -format=> 'qual');
   while (my $seq_object = $seqio_object->next_seq)
   {
      if (!exists $qual{$seq_object->display_id})
      {
         $qual{$seq_object->display_id}{'ref_qual'}=\@{$seq_object->qual()};
         if ($seq_object->desc)
            {$qual{$seq_object->display_id}{'qual_desc'}=$seq_object->desc;}
      }
      else {print $seq_object->display_id." is present more than once in the file.\n";}
   }
  return %qual;
}
#a function that trims the sequences based on the values in the description
sub  trim
{
   my $hash_ref=$_[0];     chomp($hash_ref);
   my %seqs=%{$hash_ref};
   my $hash_ref2=$_[1];     chomp($hash_ref2);
   my %quals=%{$hash_ref2};   
   my $outt=$_[2];         chomp($outt);
   my $outf=$outt.".fas";
   my $outq=$outt.".qual";
   
   use Bio::SeqIO;
   
   open(OUTF,">$outf") || die "could not open $outf due to  $!\n";
   open(OUTQ,">$outq") || die "could not open $outq due to  $!\n";
   my $tlength=0;
   foreach my $seq (sort keys %seqs)
   {
      my @temp=split/\s/, $seqs{$seq}{'desc'};
      my $start=$temp[$#temp-1]+1;
      my $stop=$temp[$#temp];
      my $length=$stop-$start;
      $tlength+=$length;
#      print $seq."\t".$length."\n";
#      print "the given coordinates are : ".$temp[$#temp-1]." - ".$temp[$#temp]."\n";
    my $seqf1=Bio::Seq->new(-id =>$seq,
                              -seq =>$seqs{$seq}{'seq'},
                              -desc =>$seqs{$seq}{'desc'});
    print OUTF ">".$seq."\n";
    print OUTF $seqf1->subseq($start,$stop)."\n";
    print OUTQ ">",$seq,"\n";
    my @scores=@{$quals{$seq}{'ref_qual'}};
    
    my $sline;
   for(my $i=$start-1;$i<$stop;$i++)
   {#print $scores[$i]." ";
      $sline.=$scores[$i]." ";
   }
      chop($sline);
#     print ">".$seq."\n";
#     print $sline."\n";
    print OUTQ $sline."\n";
   }
   print "avrage length trim : ".($tlength/(keys %seqs))."\n";
}
#the main function
sub main
{
   my $fasta=$_[0];        chomp($fasta);
   my $output=$_[1];       chomp($output);
   my $qual=$_[2];         chomp($qual);
   my %seqs=&fasta_file_parsing($fasta);
   my %quals=&qual_file_parsing($qual);
   &trim(\%seqs,\%quals,$output);
}
