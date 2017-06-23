#!/usr/bin/perl

use strict;
use warnings; 
use diagnostics;
use lib "/opt/bioperl";

#a program that tries to assemble the pair contigs
#my ($cfile,$sfile,$qfile,$tempf)=&param_check(3,"cap3_cpair_assembly.pl contig_pairs_file sequence_file quality_file fullpath_template_name");
my ($infasta,$inqual,$tempp,$tempt,$cp3p,$ffile,$opt,$cfile,@cap3op)=&param_check(2,"perl cap3_cpai_assembly.pl input_fasta_file, input_quality_file [temp folder] [tempory file template name] [path to cap3] [final output file] [option:sense/ignore] [ccontig_cluster_file] [cap3 options]");
&main($infasta,$inqual,$tempp,$tempt,$cp3p,$ffile,$opt,$cfile,\@cap3op);

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


# a simple function that goes through the sequences in a given file 
sub fasta_file_parsing
{
   use Bio::SeqIO;

   my $file=$_[0];         chomp($file);
   my %seqs;
   my $seqio_object = Bio::SeqIO->new(-file => $file, -format=>'fasta');
   while (my $seq_object = $seqio_object->next_seq)
   {
      my ($id,$type)=split/\./,$seq_object->display_id;
      if ($type)
      {
         if ($seq_object->desc)
            {$seqs{$id}{$type}{'desc'}=$seq_object->desc;}    
         $seqs{$id}{$type}{'seq'}=$seq_object->seq;         #we want the seq as it is
         $seqs{$id}{$type}{'length'}=$seq_object->length;
      }
      else
      {
         if ($seq_object->desc)
            {$seqs{$id}{'0'}{'desc'}=$seq_object->desc;}    
         $seqs{$id}{'0'}{'seq'}=$seq_object->seq;         #we want the seq as it is
         $seqs{$id}{'0'}{'length'}=$seq_object->length;
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
      my ($id,$type)=split/\./,$seq_object->display_id;
      if ($type)
      {
      if (!exists $qual{$id}{$type})
      {
      if ($seq_object->desc)
            {$qual{$id}{$type}{'desc'}=$seq_object->desc;}    
      $qual{$id}{$type}{'qual'}=\@{$seq_object->qual()};
      $qual{$id}{$type}{'length'}=$seq_object->length;
      }
      else {print $seq_object->display_id." is present more than once in the file.\n";}
   
      }
      else
      {
           if (!exists $qual{$id}{'0'})
      {
      if ($seq_object->desc)
            {$qual{$id}{'0'}{'desc'}=$seq_object->desc;}    
      $qual{$id}{'0'}{'qual'}=\@{$seq_object->qual()};
      $qual{$id}{'0'}{'length'}=$seq_object->length;
      }
      else {print $seq_object->display_id." is present more than once in the file.\n";}
   
      }
   }
  return %qual;
}


#a function that parses the cluster files
sub cfile_parser
{
   my $clfile=$_[0];       chomp($clfile);
   my @data=&file_reader($clfile);  chomp(@data);
   my %clusters;
   
   foreach my $line (@data)
   {
      my @temp=split/\t/,$line;
      for (my $i=1;$i<=$#temp;$i++)
      {
         if (!exists $clusters{$temp[0]}{$temp[$i]})
         {
            $clusters{$temp[0]}{$temp[$i]}=1;
         }
         else
         {
            print "we already have $temp[$i] in cluster $temp[0]\n";
         }
      }
   }
   return %clusters;
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
			if ($read=~ /\.b1|\.p1|\.r1/gi) # reverse strand
			{
				if  ($ace{$cont}{'read'}{$read}{'compl'} eq 'C') 
				{
					if (!$csense)
						{$csense='f1';}
					elsif ($csense ne 'f1')
						{print "We have a conflict in read orientation for read $read\n";}
				}
				elsif ( ($ace{$cont}{'read'}{$read}{'compl'} eq 'U') && (!$csense))
					{
					if (!$csense)
						{$csense='r1';}
					elsif ($csense ne 'r1')
						{print "We have a conflict in read orientation for $read\n";}
					}
			}
			elsif ($read=~ /\.g1|\.q1|\.f1/gi) #forward strand
			{
				if  ($ace{$cont}{'read'}{$read}{'compl'} eq 'U') 
				{
					if (!$csense)
						{$csense='f1';}
					elsif ($csense ne 'f1')
						{print "We have a conflict in read orientation for read $read\n";}
				}
				elsif ( ($ace{$cont}{'read'}{$read}{'compl'} eq 'C') && (!$csense))
					{
					if (!$csense)
						{$csense='r1';}
					elsif ($csense ne 'r1')
						{print "We have a conflict in read orientation for $read\n";}
					}

			}
		}
	}
	return $csense;

}

#the main function
sub main
{
   my $infas=$_[0];           chomp($infas);             # input fasta file
   my $inqual=$_[1];          chomp($inqual);          #input qualfile  
   my $tempf=$_[2];           chomp($tempf);          #tempory folder
   my $tempt=$_[3];           chomp($tempt);          #temporary filename
   my $pathcp3=$_[4];         chomp($pathcp3);      #path to Cap3
   my $out=$_[5];             chomp($out);                    #outfile name
   my $contf=$_[7];           chomp($contf);   
   my $arr_ref=$_[8];         chomp($arr_ref);            #cap3 parameters
   my @cap3=@{$arr_ref};
   my $options;
   my @singl;
   my %contigs;
   my %singlets;
   my %seqs=&fasta_file_parsing($infas);
   my %quals=&qual_file_parsing($inqual);
   my %clusters=&cfile_parser($contf);
   my $times=0;
   print "we have ".(keys %seqs)." sequence pairs \n";
   print "we have ".(keys %quals)." quality pairs \n";
   
   #the set of cap3 given options 
   foreach my $opt (@cap3)
      {$options.=$opt." ";}         chomp($options);
   print "cap3 desired options : ".$options."\n";
   
   #creating the two files of sequences to be assembled for each cluster separately
   foreach my $cluster (keys %clusters)
   {
      #opening the temporary files for sequences and quality scores for each cluster
      my $sfile=&add_slash($tempf).$tempt;
      my $qfile=&add_slash($tempf).$tempt.".qual";
      open(OUTS,">$sfile")||die "could not open $sfile for writing";  
      open(OUTQ,">$qfile")||die "could not open $qfile for writing";
      
      foreach my $id (keys %{$clusters{$cluster}})
      {
         
         foreach my $type (keys %{$seqs{$id}})
         {
            print OUTS ">".$id.".".$type."\n";
            print OUTS $seqs{$id}{$type}{'seq'}."\n";
            print OUTQ   ">".$id.".".$type."\n";
            my @temp=@{$quals{$id}{$type}{'qual'}};
            foreach my $t (@temp)
               {print OUTQ $t." ";}
            print OUTQ "\n";
         }
      }
      my $com = &add_slash($pathcp3)."cap3 ".$sfile." ".$options;
      print " comand ".$times.": ".$com." out of ".(keys %seqs)."\n";
      `$com`;
      #assinging a sense to the contig	
      my $acefile=$sfile.".cap.ace";	
      my $csense;
      my $res=`ls $acefile`;
      if ($res)
      	{$csense=&get_sense($acefile);}
      #dealing with the contings a.k.a. : going thruogh the resulting conting file and trying to write everything in anew file
      my $contfile=$sfile.".cap.contigs";
      my $qcontfile=$sfile.".cap.contigs.qual";
      my %tempo=&fasta_file_parsing($contfile);
      my %tempq=&qual_file_parsing($qcontfile);   
      #attaching the newly built conting to the cluster of sequences to be assembled 
      foreach my $seq (keys %tempo)
      {
         $contigs{$cluster}{'seq'}=$tempo{$seq}{'0'}{'seq'};
         $contigs{$cluster}{'qual'}=$tempq{$seq}{'0'}{'qual'};
         if ($csense)
         	{$contigs{$cluster}{'sense'}=$csense;}
      }
		#dealing with the singlets
		my $sinfile=$sfile.".cap.singlets";
		my %sino=&fasta_file_parsing($sinfile);
		foreach my $seq (keys %sino)
		{
			foreach my $type (keys %{$sino{$seq}})
			{
				my $id2=$id.".".$type;
				push @singl,$id2;
			}
		}
   }
	my $outc=$out.".contigs";
   open (OUT,">$outc")||die "could not open file $outc for writing\n";
   my $outql=$out.".contigs.qual";
   my $outsing=$out.".singlets";
   my $outsingq=$out.".singlets.qual";
   open (OUTQL,">$outql")||die "could not open file $outql for writing\n";
   open (OUTSL,">$outsing")||die "could not open file $outsing for writing\n";
   open (OUTSQL,">$outsingq")||die "could not open file $outsingq for writing\n";
   foreach my $seq (keys %contigs)
   {
		if ($contigs{$seq}{'sense'})
			{ print OUT ">".$seq.".".$contigs{$seq}{'sense'}."\n";
			print OUTQL ">".$seq.".".$contigs{$seq}{'sense'}."\n";}
		else
			{print OUT ">".$seq."\n";
			print OUTQL ">".$seq."\n";}
      print OUT $contigs{$seq}{'seq'}."\n";
      my @temp=@{$contigs{$seq}{'qual'}};
	   foreach my $t (@temp)
         {print OUTQL $t." ";}
      print OUTQL "\n";
   }
   foreach my $list (@singl)
   {
      my @temp=split/\./,$list;
      print OUTSL ">".$list."\n";
      print OUTSL $seqs{$temp[0]}{$temp[1]}{'seq'}."\n";
      print OUTSQL ">".$list."\n";
      my @temp2=@{$quals{$temp[0]}{$temp[1]}{'qual'}};
	   foreach my $t (@temp2)
         {print OUTSQL $t." ";}
      print OUTSQL "\n";
 
   }
   

   
   

}
