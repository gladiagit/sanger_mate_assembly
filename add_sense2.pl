#!/usr/bin/perl

use strict;
use warnings; 
use diagnostics;
# a program to change the name of the sequence into a sense containing format
my ($in,$out,$opt)=&param_check(2,"perl add_sense2.pl input_fasta/qual_file -fasta/-qual");
&main($in,$out,$opt);

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

# a function that goes through the sequences in a given file 
sub fasta_file_parsing
{
   use Bio::SeqIO;

   my $file=$_[0];         chomp($file);
   my $seqio_object = Bio::SeqIO->new(-file => $file, -format=>'fasta');
   my %seqs;
   while (my $seq_object = $seqio_object->next_seq)
   {
      if (!exists $seqs{$seq_object->display_id})
      {
         $seqs{$seq_object->display_id}{'seq'}=$seq_object->seq;
      }
      else
      {
         print "we have a double ! ".$seq_object->display_id."\n";
      }
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
#a function that renames the sequences based on their original name
sub add_sense
{
    my $hash_ref=$_[0];             chomp($hash_ref);
    my %seqs=%{$hash_ref};
    my %nseqs;
    
  #  print "inside the add_sense function\n";
    foreach my $seq (keys %seqs)
    {
        my @temp=split/_|\./, $seq;
 #       print "seqID : ".$seq." ".($#temp+1). " ".$temp[$#temp-1]."\n";
        my $id=$seq;

#        print "initial ID : ".$id."\n";
        $id=~s/(_F)(\.)*([A-Za-z]+)*(\.seq)|(_F)(\.)*([A-Za-z]+)*(\.qual)/\.b1/;
        $id=~s/(_R)(\.)*([A-Za-z]+)*(\.seq)|(_R)(\.)*([A-Za-z]+)*(\.qual)/\.g1/;
        if (exists $seqs{$seq}{'seq'})
        {
            $nseqs{$id}{'seq'}=$seqs{$seq}{'seq'};
 #           print "----seq : ".$seq." \n".$seqs{$seq}{'seq'}."\n";
  #          print "seq : ".$id." \n".$nseqs{$id}{'seq'}."\n";
        }
        elsif (exists $seqs{$seq}{'ref_qual'})
        {
            $nseqs{$id}{'ref_qual'}=$seqs{$seq}{'ref_qual'};           
        }
        
    }

    
    return %nseqs;
}

# a function that prints either fasta or qual sequnecess
sub seq_print
{
    my $ref_hash=$_[0];             chomp($ref_hash);
    my $out=$_[1];                  chomp($out);
    my %nseqs=%{$ref_hash};
    
    open(OUT,">$out")||die "could not open $out due to $!\n";

    foreach my $seq (keys %nseqs)
    {
        print OUT ">".$seq."\n";
   #     print  ">".$seq."\n";
        if (exists $nseqs{$seq}{'seq'})
        {
            print OUT $nseqs{$seq}{'seq'}."\n";
       #     print  $nseqs{$seq}{'seq'}."\n";
            
        }
        elsif (exists $nseqs{$seq}{'ref_qual'})
        {
            my @temp=@{$nseqs{$seq}{'ref_qual'}};
            foreach my $t (@temp)
                {print OUT $t." ";}
                print OUT "\n";
        }
    }
    close(OUT);
}
# the main subfunction
sub main
{
    my $file=$_[0];                 chomp($file);
    my $ofile=$_[1];                chomp($ofile);
    my $option=$_[2];               chomp($option);
    my %seqs;
    
    if ($option eq '-fasta')
    {
        %seqs=&fasta_file_parsing($file);
    }
    elsif ($option eq '-qual')
    {
        %seqs=&qual_file_parsing($file);
    }
    %seqs=&add_sense(\%seqs);
    &seq_print(\%seqs,$out);
}

