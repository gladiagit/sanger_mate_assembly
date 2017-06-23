#perl functions collection :

#!/usr/bin/perl

use strict;
use warnings; 
use diagnostics;
# an EST pipeline runing program
my ($cfile)=&param_check(1,"perl est_pipeline.pl config_file");
&main($cfile);

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
#a function for reading in the configuration file
#config format
#fasta_file=full_filename_to_fastafile
#qual_file=full_filename_to_quality_score_file
#vector_file=filename_to_vector
#vector_splice_filename_to_vector_splice
sub config_file_reader
{
   my $file=$_[0];      chomp($file);
   my %setings;
   my @data=&file_reader($file); chomp(@data);
   
#   print "inside configfile\n";
#   print" we have setings : ".(keys %setings)."\n";
   foreach my $line (@data)
   {
     my @temp=split/\t/,$line;
     $setings{$temp[0]}=$temp[1];
#     print "line : ".$line."\n";
   }
#   foreach my $seting (keys %setings)
#   {
#      print $seting." : ".$setings{$seting}."\n";
#   }
   return %setings;
}
#a function that takes a file and returns its path
sub get_path
{
   my $filename=$_[0];     chomp($filename);
   my @temp= split /\//, $filename;
   my $path;
   for (my $i=0;$i<$#temp;$i++)
   {#print $temp[$i] ."\n";
    $path.=$temp[$i]."/";}
   #print "path to folder : ". $path."\n";
   return $path;
}
# a file that given a full path returns a filename
#requires get_path function
sub get_filename
{
   my $cfile=$_[0];  chomp($cfile);
   my $ttemp=&get_path($cfile);
   my $temp=$cfile;
   $temp=~ s/$ttemp//;
   return $temp;
}
# a function that creates a filename with full path based on the option/ step in the pipeline and the original filename
sub frename
{
   my $ofname=$_[0];       chomp($ofname);            #the original filename
   my $opt1=$_[1];         chomp($opt1);    
   my $opt2=$_[2];         chomp($opt2);    
   my $temp=&get_filename($ofname);
   
 #  print "before : ".$temp."\n";
   
   if ($opt1 eq 'fasta')
   {
      $temp=~ s/\.fasta|\.fst|\.fas|\.fa|\.seq|\.sqn//;
   }
   elsif ($opt1 eq 'qual')
   {
      $temp=~ s/\.qual|\.qul//;
   }
   elsif ($opt1 eq 'log')
   {
      $temp=~ s/\.log|\.clean|\.cln//;
   }
   
   if ($opt2 eq 'lucy')
   {
      $temp.='_lucy';  
   }
   elsif ($opt2 eq 'seqclean')
   {
      $temp.='_sqcln';
   }
   elsif ($opt2 eq 'trimmed')
   {
      $temp.="_trim";
   }
#   print "after : ".$temp."\n";

}
#sub main
sub main
{
   my $configf=$_[0];      chomp($configf);
#   print " main function\n";
   my %setings=&config_file_reader($configf);
#   print "main files of interest : \n".$setings{'fasta_file'}."\n".$setings{'qual_file'}."\n";
   $setings{'lucy_out_fasta'}=&add_slash($setings{'temp_path'}).&frename($setings{'fasta_file'},'fasta','lucy').".fas";
   $setings{'lucy_out_qual'}=&add_slash($setings{'temp_path'}).&frename($setings{'qual_file'},'qual','lucy').".qual";
   
   my $com=&add_slash($setings{'lucy_path'})."lucy ".$setings{'fasta_file'}." ".$setings{'qual_file'}." -vector ".$setings{'vector_file'}." ".$setings{'splice_file'}.
                                                                                                      " -output ".$setings{'lucy_out_fasta'}." ".$setings{'lucy_out_qual'};
                                                                                                    #." ".$setings{'lucy_extra_params'}; #to be added later (10/25/11)
   $setings{'lucy_trimmed_fasta'}=&add_slash($setings{'temp_path'}).&frename($setings{'lucy_out_fasta'},'fasta','trimmed').".fas";;
   $setings{'lucy_trimmed_qual'}=&add_slash($setings{'temp_path'}).&frename($setings{'lucy_out_qual'},'qual','trimmed').".qual";
   print "\ncom : ".$com."\n";
   `$com`;       # where the lucy quallity score based trimmming is being initiated
  
  $com="perl lucy_trim.pl ".$setings{'lucy_out_fasta'}." ".$setings{'lucy_out_qual'}." ". &add_slash($setings{'temp_path'}).&frename($setings{'lucy_out_fasta'},'fasta','trimmed');
   print "\ncom : ".$com."\n";
   `$com`;      # where the quality based trimming sugested in lucy_triming step is actually being executed

   $setings{'sqcln_fasta_out'}=&add_slash($setings{'temp_path'}).&frename($setings{'lucy_trimmed_fasta'},"fasta","seqclean").".fas";
   print " \t".$setings{'sqcln_fasta_out'}."\n";
   $setings{'sqcln_report_out'}=&add_slash($setings{'temp_path'}).&frename($setings{'lucy_trimmed_fasta'},"fasta","seqclean").".cln";
   $com=&add_slash($setings{'seqclean_path'})."seqclean ".$setings{'lucy_trimmed_fasta'}." -v ".$setings{'UniVec_path'}."UniVec,".$setings{'vector_file'}.",".$setings{'extra_clean'}
                                                                                        ." -o ".$setings{'sqcln_fasta_out'}
                                                                                        ." -r ".$setings{'sqcln_report_out'};
   print "\ncom : ".$com."\n";
   `$com`;
   $setings{'sqcln_qual_out'}=$setings{'lucy_trimmed_qual'}.".clean";
   $com=&add_slash($setings{'seqclean_path'})."cln2qual ".$setings{'sqcln_report_out'}." ". $setings{'lucy_trimmed_qual'};
   print "\ncom : ".$com."\n";
   `$com`;
   
   $com="perl cap3_mate_assembly.pl ".$setings{'fasta_file'}." ".$setings{'qual_file'}." ".$setings{'cap3_mate_params'};
    print "com : ".$com."\n";
#   `$com`;
   
}