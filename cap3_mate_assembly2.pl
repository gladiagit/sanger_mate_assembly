#!/usr/bin/perl

use strict;
use warnings; 
use diagnostics;
use lib "/opt/bioperl";
# a program that uses cap3 to assemble one pair of sequences (b1 & g1) at a time
my ($inf,$inq,$pcap3,$tempp,$outt,@options)=&param_check(5,"perl cap3_mate_assembly2.pl input_fasta_file input_qual_file path_to_cap3 temporary_files_path output_filename_template
                                                                                  [cap3_specific_parameters]");
my %seqs;
&main($inf,$inq,$pcap3,$tempp,$outt,\@options);

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
# a simple function that goes through the sequences in a given file 
sub fasta_file_parsing
{
   use Bio::SeqIO;
   
   my $file=$_[0];         chomp($file);
   my $opt=$_[1];          if ($opt)   {chomp($opt);}
   my %tseqs;
   my $seqio_object = Bio::SeqIO->new(-file => $file, -format=>'fasta');

   while (my $seq_object = $seqio_object->next_seq)
   {
      my ($id,$type)=split/\./,$seq_object->display_id;
      
      if (!$type)
      {$type='none';}
         
         if (($seq_object->desc) && (($opt)&&($opt eq 'gen')))           # the sequences are stored in a global variable called 'seqs'
            {$seqs{$id}{'type'}{$type}{'desc'}=$seq_object->desc;}
         elsif  (($seq_object->desc) && (!$opt))
            {$tseqs{$id}{'type'}{$type}{'desc'}=$seq_object->desc;}
         if (($opt)&&($opt eq 'gen'))
         {
            $seqs{$id}{'type'}{$type}{'seq'}=$seq_object->seq;         #we want the seq as it is
            $seqs{$id}{'type'}{$type}{'length'}=$seq_object->length;
         }
         else
         {
            $tseqs{$id}{'type'}{$type}{'seq'}=$seq_object->seq;         #we want the seq as it is
            $tseqs{$id}{'type'}{$type}{'length'}=$seq_object->length;
         }

   }
   if (!$opt)
      {return %tseqs;}
}

# file that parses a quality file
sub qual_file_parsing
{
   use Bio::SeqIO;
   my $sfile=$_[0];      chomp($sfile);   #the sequence file
   my $opt=$_[1];        if ($opt)  {chomp($opt)}
   my %tseqs;
   my $seqio_object = Bio::SeqIO->new(-file => $sfile , -format=> 'qual');
   
   while (my $seq_object = $seqio_object->next_seq)
   {
      my ($id,$type)=split/\./,$seq_object->display_id;
      
      if (!$type)
         {$type ='none';}
      if ($seq_object->desc)
      {
         if (($opt)&&($opt eq 'gen'))
            {$seqs{$id}{'type'}{$type}{'desc'}=$seq_object->desc;}
         else
            {$tseqs{$id}{'type'}{$type}{'desc'}=$seq_object->desc;}
      }
      if (($opt)&&($opt eq 'gen'))
      {
         $seqs{$id}{'type'}{$type}{'qual'}=\@{$seq_object->qual()};
         $seqs{$id}{'type'}{$type}{'length'}=$seq_object->length;
      }
      else
      {
         $tseqs{$id}{'type'}{$type}{'qual'}=\@{$seq_object->qual()};
         $tseqs{$id}{'type'}{$type}{'length'}=$seq_object->length;
      }
   }
   if (!$opt)
      {return %tseqs;}
}
# function that takes a filename and reads in and returns a the data in the file as an array of lines 
sub file_reader
{
   my $filename=$_[0];  chomp($filename);
   open(INP,"<$filename")||die "could not open $filename for reading $!\n";
   my @data=<INP>;
   return @data;

}


#a function that adds a slash to a path if necesary
sub add_slash
{
   my $word=$_[0];   chomp($word);
   my $npath;
   my @str=split//, $word;
   
   if ($str[$#str] ne '/')
      { $npath=$word."/";}
   else { $npath=$word; }
   return $npath;
}

# a function that creates the temporary pairs files
sub create_tfiles
{
   my $tempf=$_[0];           chomp($tempf);
   
   foreach my $id (keys %seqs)
   {
      my $sfile=&add_slash($tempf).$id;
      my $qfile=&add_slash($tempf).$id.".qual";
         $seqs{$id}{'sfile'}=$sfile;                        #the temporary sequence file
         $seqs{$id}{'qfile'}=$qfile;                        #the temporary quality file
         $seqs{$id}{'afile'}=$sfile.".cap.ace";             #the generated ace file
         $seqs{$id}{'csfile'}=$sfile.".cap.contigs";        #the generated contigs fasta file
         $seqs{$id}{'cqfile'}=$sfile.".cap.contigs.qual";   #the generated contigs quality file
      
      open(OUTS,">$sfile")||die "could not open $sfile for writing";  
      open(OUTQ,">$qfile")||die "could not open $qfile for writing";
      
      foreach my $type (keys %{$seqs{$id}{'type'}})
      {
         print OUTS ">".$id.".".$type."\n";
         print OUTS $seqs{$id}{'type'}{$type}{'seq'}."\n";
         my @temp=@{$seqs{$id}{'type'}{$type}{'qual'}};
         print OUTQ ">".$id.".".$type."\n";
	     foreach my $t (@temp)
            {print OUTQ $t." ";}
            print  OUTQ "\n";
      }

   }
}

#a function that goes through all the sequences and tries to assemble them
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

#a function that deals with the sense of the consensus
sub get_sense
{
	my $acefile=$_[0];			chomp($acefile);
	my ($nc,$nr,$href)=&ace_parser($acefile);
	my %ace=%{$href};
	my $csense;
		foreach my $cont (keys %ace)
	{
		foreach my $read (keys %{$ace{$cont}{'read'}})
		{
			if ($read=~ /\.b1|\.p1/gi) # reverse strand
			{
				if  ($ace{$cont}{'read'}{$read}{'compl'} eq 'C') 
				{
					if (!$csense)
						{$csense='f1';}
					elsif ($csense ne 'f1')
						{print "We have a conflict in read orientation for read $read - sense $csense\n";
                         $csense="null";}
				}
				elsif ( ($ace{$cont}{'read'}{$read}{'compl'} eq 'U') && (!$csense))
					{
					if (!$csense)
						{$csense='r1';}
					elsif ($csense ne 'r1')
						{print "We have a conflict in read orientation for $read - sense + $csense\n";
                         $csense="null";}
					}
			}
			elsif ($read=~ /\.g1|\.q1/gi) #forward strand
			{
				if  ($ace{$cont}{'read'}{$read}{'compl'} eq 'U') 
				{
					if (!$csense)
						{$csense='f1';}
					elsif ($csense ne 'f1')
						{print "We have a conflict in read orientation for read $read -| sense $csense\n";
                         $csense="null";}
				}
				elsif ( ($ace{$cont}{'read'}{$read}{'compl'} eq 'C') && (!$csense))
					{
					if (!$csense)
						{$csense='r1';}
					elsif ($csense ne 'r1')
						{print "We have a conflict in read orientation for $read - sense * $csense\n";
                         $csense="null";}
					}

			}
		}
	}
	return $csense;

}

#the main function
sub main
{
   my $fin=$_[0];             chomp($fin);                           #input fasta file
   my $qin=$_[1];             chomp($qin);                           #input quality file
   my $cap3p=$_[2];           chomp($cap3p);                         #cap3 path
   my $tempf=$_[3];           chomp($tempf);                         #temporary files path
   my $outft=$_[4];           chomp($outft);                         #results output_filename template
   my $arr_ref=$_[5];         if ($arr_ref)  {chomp($arr_ref);}      #cap3 options array reference
   my @cap3o=@{$arr_ref};                          #cap3 options 
   my $options;                                    #cap3 options in string format 
   
 print "outft : ".$outft."\n";   
   my $ofilecs=$outft."_contigs";           #fasta file with assembled sequences
   my $ofilecq=$outft."_contigs.qual";      #qual file with assembled sequences quality scores
   my $ofiless=$outft."_singlets";          #fasta file with unassembled sequences
   my $ofilesq=$outft."_singlets.qual";     #qual file with assembled sequences quality scores
   
   open(OUTCS,">$ofilecs")||die "could not open $ofilecs due to $!\n";
   open(OUTCQ,">$ofilecq")||die "could not open $ofilecq due to $!\n";
   open(OUTSS,">$ofiless")||die "could not open $ofiless due to $!\n";
   open(OUTSQ,">$ofilesq")||die "could not open $ofilesq due to $!\n";
   
   &fasta_file_parsing($fin,'gen');          #placing the sequences and quality scores in an single general hash
   &qual_file_parsing($qin,'gen');
   foreach my $opt (@cap3o)
      {$options.=$opt." ";}         chomp($options);
   
   print "cap3 desired options : ".$options."\n";
   print "we have ".(keys %seqs)." pairs \n";
   
   &create_tfiles($tempf);
   foreach my $id (keys %seqs)
   {
      my $com = &add_slash($cap3p)."cap3 ".$seqs{$id}{'sfile'}." ".$options;
      `$com`;
      my ($ncont,$nreads,$href)=&ace_parser($seqs{$id}{'afile'});
      if ( ($ncont>0) && ($ncont<2))
      {
         my %temps=&fasta_file_parsing($seqs{$id}{'csfile'});
         my %tempq=&qual_file_parsing($seqs{$id}{'cqfile'});
         my $csense=&get_sense($seqs{$id}{'afile'});             #the sense of the newly constructed contig
         use Data::Dumper;
#         print Dumper($csense);

#         print $id." : ".$csense."\n";
         $seqs{$id}{'contig'}{'seq'}=$temps{'Contig1'}{'type'}{'none'}{'seq'};
         $seqs{$id}{'contig'}{'qual'}=$tempq{'Contig1'}{'type'}{'none'}{'qual'};
         $seqs{$id}{'contig'}{'sense'}=$csense;
         print OUTCS ">".$id.".".$seqs{$id}{'contig'}{'sense'}."\n";
         print OUTCS $seqs{$id}{'contig'}{'seq'}."\n";
         print OUTCQ ">".$id.".".$seqs{$id}{'contig'}{'sense'}."\n";
         my @temp=@{$seqs{$id}{'contig'}{'qual'}};
	     foreach my $t (@temp)
            {print OUTCQ $t." ";}
            print OUTCQ "\n";
      }
      elsif ($ncont==0)
      {
         foreach my $type (keys %{$seqs{$id}{'type'}})
         {
            print OUTSS ">".$id.".".$type."\n";
            print OUTSS $seqs{$id}{'type'}{$type}{'seq'}."\n";
         
            print OUTSQ ">".$id.".".$type."\n";
            my @temp=@{$seqs{$id}{'type'}{$type}{'qual'}};
            foreach my $t (@temp)
               {print OUTSQ $t." ";}
               print OUTSQ "\n";
         }
      }
      elsif ($ncont>=2)
      {
         print "!!!!!!!!!!!!!!!!!!!!!!Warning!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
         print "mate assembly step generates more than one assembled Contig for sequence ID : ".$id."\n";
      }
      
   }
   
}