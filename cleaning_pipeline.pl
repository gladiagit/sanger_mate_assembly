#!/usr/bin/perl

use strict;
use warnings; 
use diagnostics;
#a program that attempts to put toghether the lucy and seqclean cleaning algorithm
my ($ffile,$qfile,$vfile,$spfile,$cfile,$tempf,$lpath,$sqpath,$vdbfile,@params)=&param_check(8,"perl cleaning_pipeline.pl input_fastafile input_quality_file vectorfile splicefile contaminantfile{untrimed vector/linker/adapter} temporary_folder path_to_lucy path_to_seqclean vectordb [lucy/seqclean params : -lucy:lucy_param:value / -sqcln:sqcln_param:value]");
&main($ffile,$qfile,$vfile,$spfile,$cfile,$tempf,$lpath,$sqpath,$vdbfile,\@params);

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

# the main function
sub main
{
   my $fastaf=$_[0];              chomp($fastaf);
   my $qualf=$_[1];               chomp($qualf);
   my $vectorf=$_[2];             chomp($vectorf);
   my $splicef=$_[3];             chomp($splicef);
   my $contf=$_[4];               chomp($contf);
   my $tempf=$_[5];               chomp($tempf);
   my $lucyp=$_[6];               chomp($lucyp);
   my $sqclnp=$_[7];              chomp($sqclnp);
   my $vectordb=$_[8];            chomp($vectordb);
   my $arr_ref=$_[9];             chomp($arr_ref);
   my @params=@{$arr_ref};
   
   my $lparam='';
   my $sparam='';
   my $com1='';   
   my $com2='';   
   my $com3='';
   my $com4='';
   my $outt=&add_slash($tempf)."intermediate";
   
   foreach my $set (@params)
   {
      my @temp=split/:/,$set;
      if ($temp[0] eq 'lucy')
      {if ($temp[1]!~ /-/gi)
         {$lparam.="-".$temp[1]." ".$temp[2]." ";} 
      else
         {$lparam.=$temp[1]." ".$temp[2]." ";} 
      }
      elsif ($temp[0] eq 'sqcln')
      {if ($temp[1]!~ /-/gi)
         {$sparam.="-".$temp[1]." ".$temp[2]." ";} 
      else
         {$sparam.=$temp[1]." ".$temp[2]." ";} 
      }
   }
   $com1=&add_slash($lucyp)."lucy ".$fastaf." ".$qualf." -v ".$vectorf." ".$splicef." -o "
               .&add_slash($tempf)."lucy_temp.fas ".&add_slash($tempf)."lucy_temp.qual ".$lparam;
   $com3='perl lucy_trim.pl '.&add_slash($tempf)."lucy_temp.fas ".&add_slash($tempf)."lucy_temp.qual ".$outt;
   $com2=&add_slash($sqclnp)."seqclean ".$outt.".fas -v ".$vectordb.",".$contf." -r ".$outt.".fas.cln"." ".$sparam;
   $com4=&add_slash($sqclnp)."cln2qual ".$outt.".fas.cln ".$outt.".qual";
   print "com1 : ".$com1."\n";
   print "com3 : ".$com3."\n";
   print "com2 : ".$com2."\n";
   print "com4 : ".$com4."\n";
   `$com1`;
   `$com3`;
   `$com2`;
   `$com4`;
}

