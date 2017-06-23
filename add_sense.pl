#!/usr/bin/perl

use strict;
use warnings; 
use diagnostics;
# a program that takes in a set of sequences and a list of sense and attaches the sense information to the sequences 
my ($sfile,$ffile,$ofile,$opt)=&param_check(3,"perl add_sense.pl  sense_file fasta_file output_file [all] [list]");
&main($sfile,$ffile,$ofile,$opt);
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

#a function that takes an array as refernce and converts it to a hash
sub array_to_hash
{
   my $array_ref=$_[0];                chomp($array_ref);
   my @array=@{$array_ref};
   my %hash;
   
   foreach my $data (@array)
   {
      if ($data=~ /\s|\t/gi)
       {
         my @temp=split/\s|\t/, $data;
         if (!exists $hash{$temp[0]})
            {$hash{$temp[0]}=$temp[1];}
         else {print "we have a double : $temp[0] - ".$hash{$temp[0]}."\n";}
       }
      else
      {
         if (!exists $hash{$data})
            {$hash{$data}=1;}
         else {print "we have a double : $data - ".$hash{$data}."\n";}
         
      }
   }
return %hash;
}
# a simple function that goes through the sequences in a given file 
sub fasta_file_parsing
{
   use Bio::SeqIO;

   my $file=$_[0];         chomp($file); 
   my $out=$_[1];        chomp($out);
   my $hash_ref=$_[2];  chomp($hash_ref);
   my $option=$_[3];    chomp($option);
   my %senses=%{$hash_ref};
   open(OUT,">$out")||die "could not open $out due to $!\n";
   print "in the fastafile function\n";
   if ($opt ne 'list')
   {
      my $seqio_object = Bio::SeqIO->new(-file => $file, -format=>'fasta');
      while (my $seq_object = $seqio_object->next_seq)
      {
         if (!exists $senses{$seq_object->display_id})
         {
            print "no sense information available for ".$seq_object->display_id."\n";
            if ($option eq 'all')
               {print OUT ">".$seq_object->display_id."\n";
               print OUT $seq_object->seq."\n";}
         }
         else
         {
            print OUT ">".$seq_object->display_id.".".$senses{$seq_object->display_id}."\n";
            print OUT $seq_object->seq."\n";
         }
      }
   }
   else
   {
      print "down the list option...\n";
      my @ids=&file_reader($file);        chomp(@ids);
      print "we have ".($#ids+1)." seqids\n";
      foreach my $id (@ids )
      {
         if ($senses {$id})
         {
            print OUT $id.".".$senses{$id}."\n";
#            print $id.".".$senses{$id}."\n";
         }
         else
         {
            print " no sense for  : ".$id."\n";
            print OUT $id."\n";
         }
      }
   }
}
# the main subfunction
sub main
{
   my $sense=$_[0];                           chomp($sense);
   my $fasta=$_[1];                             chomp($fasta);
   my $out=$_[2];                               chomp($out);
   my $opt=$_[3];                         chomp($opt);
   my @data=&file_reader($sense);         chomp(@data);
   my %hsense=&array_to_hash(\@data);
   &fasta_file_parsing($fasta,$out,\%hsense,$opt);
}
