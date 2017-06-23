#!/usr/bin/perl

use strict;
use warnings; 
use diagnostics;
# aprogram that parses a .ace file
my ($infile)=&param_check(1,"perl ace_parser.pl input_ace_file");

&main($infile);

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
#the main subfunction
sub main
{
   my $acefile=$_[0];            chomp($acefile);
   my ($ncont,$nreads,$arr_ref)=&ace_parser($acefile);
   my %ace=%{$arr_ref};
   print "we have $ncont contigs in the file composed of $nreads reads in total\n";
}
# the ace parser
sub ace_parser
{
   my $infile=$_[0];                      chomp($infile);
   my @data=&file_reader($infile);   chomp(@data);
   my %ace;
   my $ncontigs;
   my $nreads;
   my $ccont;
   my $cseq;
   my $cqual;
   my $frag=0;
   my $isrd=0;
   my $cread;
   for (my $i=0;$i<=$#data;$i++)
 #  for (my $i=0;$i<=350;$i++)
   {
      my @line=split/\s/,$data[$i];
      my $nwords=$#line+1;
 #     print "line $i : ",$data[$i],"\n";
#      print " we have ".$nwords." words on line ".$i."\n";
      if ($nwords>0)
      {
         if ($line[0] eq 'AS')
         {
            $ncontigs=$line[1];
            $nreads=$line[2];
         }
         elsif ($line[0] eq 'CO')
         {
            if ($ccont) # placing the completed seq and qual of the previous contig inot the hash
            {
               $ace{$ccont}{'raw_seq'}=$cseq;
               chop($cqual);
               $ace{$ccont}{'raw_qual'}=$cqual;
            }
            $ccont=$line[1]; #re-initializing the ccont, cseq and cqual variables
            $cseq='';
            $cqual='';
            $frag=0;
            $isrd=0;
            $cread='';
            if (!exists $ace{$ccont})
            {
               $ace{$ccont}{'length'}=$line[2];
               $ace{$ccont}{'nreads'}=$line[3];
               $ace{$ccont}{'base_segmets'}=$line[4];
               $ace{$ccont}{'phrap_compl'}=$line[5];
            }
            else  {print "contig ".$ccont," is already in the ace hash\n";}
         }
         elsif ( ($nwords==1) &&($line[0]!~ /[0-9]|[a-z]/ ) && ($line[0] ne 'BQ') && ($line[0] ne 'DS') && (!$isrd))
         {
            $cseq.=$line[0];
         }
         elsif ($line[0]!~ /[A-Z]|[a-z]/ )
         {
            foreach my $score (@line)
            {
               chomp($score);
               $cqual.=$score." ";
            }
         }
         elsif ($line[0] eq 'AF')
         {
            $ace{$ccont}{'read'}{$line[1]}{'compl'}=$line[2];
            $ace{$ccont}{'read'}{$line[1]}{'start'}=$line[3];
         }
         elsif ($line[0] eq 'BS')
         {
            $ace{$ccont}{'padding'}{$line[3]}{$frag}{'startp'}=$line[1];
            $ace{$ccont}{'padding'}{$line[3]}{$frag}{'endp'}=$line[2];
            $frag++;
         }
         elsif ($line[0] eq 'RD')
         {
            $isrd=1;
            $cread=$line[1];
            $ace{$ccont}{'read'}{$cread}{'np_bp'}=$line[2]; #number of paded bases - i have no idea what this means , but the parser should contain all the data in the ACE file
            $ace{$ccont}{'read'}{$cread}{'n_whole_reads'}=$line[3];
            $ace{$ccont}{'read'}{$cread}{'tags'}=$line[4];
         }
         elsif ($line[0] eq 'QA')
         {
            $ace{$ccont}{'read'}{$cread}{'qcstart'}=$line[1]; #quality cliping start and end
            $ace{$ccont}{'read'}{$cread}{'qcend'}=$line[2];
            $ace{$ccont}{'read'}{$cread}{'acstart'}=$line[3]; #align cliping start and end
            $ace{$ccont}{'read'}{$cread}{'acend'}=$line[4];
         }
         elsif (($nwords>1)&&($line[0] eq 'DS'))
         {
            $ace{$ccont}{'cromat_file'}=$line[2];
            $ace{$ccont}{'phd_file'}=$line[4];
            $ace{$ccont}{'time'}=$line[6];
            $ace{$ccont}{'chem'}=$line[8];
            $ace{$ccont}{'dye'}=$line[10];
            $ace{$ccont}{'template'}=$line[12];
            $ace{$ccont}{'direction'}=$line[14];
            
         }
         
      }
   }
#   print "the total number of contigs in the file : ".$ncontigs."\n";
#   print "the total number of reads in the file : ".$nreads."\n";
#   print "csq :\n".$cseq."\n";
   #chop($cqual);
   #print "cq :\n".$cqual.".\n";
   return $ncontigs,$nreads,\%ace;
}
