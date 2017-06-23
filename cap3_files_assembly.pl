#!/usr/bin/perl

use strict;
use warnings; 
use diagnostics;
#a program that takes in a file templatename a list of cap3 parameters and assembles all the files in the path with template naame using the param
my ($tempn, $capp,@opts)=&param_check(2,'perl cap3_files_assembly.pl template_name path_to_cap3 cap3_opts');
&main($tempn, $capp,\@opts);

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
#the main subfunction
sub main
{
   my $tname=$_[0];                 chomp($tname);
   my $cpath=$_[1];                  chomp($cpath);
   my $arr_ref=$_[2];                 chomp($arr_ref);
   my @cap3=@{$arr_ref};
   my $com='ls '.$tname."*";
   my @files=`$com`;                 chomp(@files);
   my $options;
   
#   print "first file : ".$files[0]."\n";
#   my @temp=split/$tname|\.txt|\.fasta|\.fas|\.fa/, $files[0];
#   print "insertion name : ".$temp[1]."\n";
   foreach my $opt (@cap3)
      {$options.=$opt." ";}         chomp($options);
 #  print "cap3 desired options : ".$options."\n";   
 
   foreach my $file (@files)
   {
      my $com2=&add_slash($cpath)."cap3 ".$file." ".$options;
      print "com2= ".$com2."\n";
      `$com2`;
   }
   
}
