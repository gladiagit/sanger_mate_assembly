#!/usr/bin/perl

use strict;
use warnings; 
use diagnostics;
#a blast result parser that can hopefully do anything
# -nhits:kword:e-val should return the number of sequences in the file that have the kword in the hit descrion/name with at least an evalue of e-val 
my ($bfile,@opt)=&param_check(1,"perl blast_parser_complete.pl blast_result_file\n
                              [-bhit_specie_dist:[eval]:outputfile:[specie]]
                              [-nhits:kword:e-val]
                              [-bhit_pident:outputfile]
                              [-ident:second_bres_file:outputfile:{all/max/min}]
                              [-cluster_2nd_file:second_blastfile:keyword:outputfile]
                              [-bhit_list:seqs_nreads_file:hits_nreads_file:output_file]
                              [-split_id]
                              [-scount:template_filename]
                              [-combine_go:GO_annnotationfile:option{yes/no}]
                              [-wordfq_list:e-val]
                              [-bhit_pident_specie:outputfile]
                              [-bhit_cluster:outputfile:[e-val]:[nseq]:[list]]
                              [-bhit_list_clean:outputfile:[e-val]:[id_only]]
                              [-map_to_bhit:[e-value]]
                              [-pident_interval]
                              [-cluster_pident:%identity]
                              [-cluster_bhit_val:clusterfile:outputfile:[ignore_sense]]
                              [-no_hits:[n/ids/id_print]]");
&main($bfile,\@opt);

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
# the full blast parser function
sub blast_results_parser
{
   my $filename=$_[0];        chomp($filename);
   my $tnhits=$_[1];             # only looking at first hit
  
   use Bio::SearchIO;
   my $bres=Bio::SearchIO->new(-format=>"blast",-file=>$filename);
   my %blres;
   my %hits;
   my @queryl;
   while( my $result = $bres->next_result )        #going through each result in the file 
   {
       my $chit=0;
         push @queryl,$result->query_accession;
      while( my $hit = $result->next_hit )         #going through each hit
      {
         $chit++;
         if (($tnhits)&&($chit>=2))          # only looking at the best hit, this was done to optimize running time for large blast files
         {next;}
         my $cfrag=1;
         $blres{$result->query_accession}{'algorithm'}=$result->algorithm;
         $blres{$result->query_accession}{'length'}=$result->query_length;
         $blres{$result->query_accession}{'description'}=$result->query_description;
         $blres{$result->query_accession}{'name'}=$result->query_name;
         $blres{$result->query_accession}{'nhits'}=$result->num_hits;
         $blres{$result->query_accession}{'hits'}{$hit->accession}{'name'}=$hit->name;
         $blres{$result->query_accession}{'hits'}{$hit->accession}{'description'}=$hit->description;
         $blres{$result->query_accession}{'hits'}{$hit->accession}{'length'}=$hit->length;
         $blres{$result->query_accession}{'hits'}{$hit->accession}{'e-value'}=$hit->significance;
         $blres{$result->query_accession}{'hits'}{$hit->accession}{'raw_score'}=$hit->raw_score;
         $blres{$result->query_accession}{'hits'}{$hit->accession}{'rank'}=$chit;
         while (my $hsp = $hit ->next_hsp)   #going through each field for each hit
         {
            $hits{$hit->accession}{$result->query_accession}{'frags'}{$cfrag}{'end'}=$hsp->end('hit');
            $hits{$hit->accession}{$result->query_accession}{'frags'}{$cfrag}{'start'}=$hsp->start('hit');
            $hits{$hit->accession}{$result->query_accession}{'eval'}=$hit->significance;
            $hits{$hit->accession}{$result->query_accession}{'description'}=$hit->description;
            $hits{$hit->accession}{$result->query_accession}{'raw_score'}=$hit->raw_score;
            $hits{$hit->accession}{$result->query_accession}{'rank'}=$chit;
          #  print "hit description : ".$hits{$hit->accession}{$result->query_accession}{'description'}."\n";
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'e-value'}=$hsp->evalue;
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'fconserved'}=$hsp->frac_conserved;
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'fident'}=$hsp->frac_identical;
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'tlength'}=$hsp->length;               #total lenght of alignment, including gaps
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'hlength'}=$hsp->length('hit');        #total lenght of hit, minus gaps
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'qlength'}=$hsp->length('query');      #total lenght of query, minus gaps
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'nrconserved'}=$hsp->num_conserved;    #number of conserved residues
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'nrident'}=$hsp->num_identical;        #number of identical residues
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'score'}=$hsp->score;
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'pident'}=$hsp->percent_identity;
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'qframe'}=$hsp->query->frame;
           # print "query frame : ".$hsp->query->frame."\n";
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'hframe'}=$hsp->hit->frame;
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'qstrand'}=$hsp->strand('query');
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'hstrand'}=$hsp->strand('hit');
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'hstart'}=$hsp->start('hit');          #start position from alignment for hit
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'qstart'}=$hsp->start('query');        #start position from alignment for query
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'hend'}=$hsp->end('hit');              #end position from alignment for hit
            $blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'qend'}=$hsp->end('query');            #end position from alignment for query
            push @{$blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'qident_pos'}},$hsp->seq_inds('query','identical');        #the positions of ifentities in query
            push @{$blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'qcons_pos'}},$hsp->seq_inds('query','conserved-not-identitical');        #the positions of conserved in query
            push @{$blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'hident_pos'}},$hsp->seq_inds('hit','identical');        #the positions of ifentities in hit
            push @{$blres{$result->query_accession}{'hits'}{$hit->accession}{'frags'}{$cfrag}{'hcons_pos'}},$hsp->seq_inds('hit','conserved-not-identitical');        #the positions of conserved in hit
            $cfrag++;
         }   
      }
   }
   return \%blres,\%hits,\@queryl;
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
#the main subfunction
sub main
{
   my $blastf=$_[0];       chomp($blastf);
   my $arr_ref=$_[1];      chomp($arr_ref);
   my @options=@{$arr_ref};
   my $nhits;
   foreach my $option (@options)
   { if (($option=~ /bhit/gi)||($option eq '-no_hits')) {$nhits=1;}}
   my ($bres_ref,$hit_ref,$ql_ref)=&blast_results_parser($blastf,$nhits);
   my %ress=%{$bres_ref};
   my %hress=%{$hit_ref};
   my %query_ids;
   print "total number of results  keys : ".(keys %ress)." number of options : ".($#options+1)."\n";
   foreach my $option (@options)
   {
      my @temp=split/:/,$option;
      if (($temp[0] eq '-nhits') && (($temp[1])&&($temp[2])))
      {
         my ($nres,$gl_ref)=&nresults_count(\%ress,$temp[1],$temp[2]);
         if (($temp[3])&&($temp[3] eq 'list'))
         {
            my %genel=%{$gl_ref};
            foreach my $gene  (keys %genel)
            {
               print $gene."\t".$genel{$gene}."\n";
            }
         }
         print "we have ".$nres." results with ".$temp[1]." keyword at a e-val less then ".$temp[2]."\n";
      }
      elsif (($temp[0] eq '-ident') &&($temp[2]))
      {
         my ($ress2_ref,$hress2)=&blast_results_parser($temp[1]);
         my %ress2=%{$ress2_ref};
         my %hits1=&get_ident(\%ress);
         my %hits2=&get_ident(\%ress2);
         print " have ".(keys %hits1)." results in the first hash and ".(keys %hits2)." in the second\n";
         open(OUT,">$temp[2]")||die "could not open $temp[2] due to $!\n";
         
         foreach my $hit (sort keys %hits1)
         {
            if (exists $hits2{$hit})
            {
               my $min=100;
               my $max=0;
               my $seqmax;
               my $seqmin;
               
               foreach my $query1 (keys %{$hits1{$hit}{'query'}})
               {
                  if ((!$temp[3])||($temp[3] eq 'all'))
                  {
#                     print $hits1{$hit}{'query'}{$query1}{'rpident'}."\t".$query1."\t";
                     print OUT $hits1{$hit}{'query'}{$query1}{'rpident'}."\t".$query1."\t";
                  }
                  elsif ( ($temp[3] eq 'max') &&($max<=$hits1{$query1}{'rpident'}))
                  {
                     $max=$hits1{$query1}{'rpident'};
                     $seqmax=$query1;
                  }
                  elsif ( ($temp[3] eq 'max') &&($min>=$hits1{$query1}{'rpident'}))
                  {
                     $min=$hits1{$query1}{'rpident'};
                     $seqmin=$query1;
                  }
                  
               }
               if (($temp[3]) &&($temp[3] eq 'max'))
               {
  #                print $hits1{$hit}{'query'}{$seqmax}{'rpident'}."\t".$seqmax."\t";
                  print OUT $hits1{$hit}{'query'}{$seqmax}{'rpident'}."\t".$seqmax."\t";
               }
               elsif (($temp[3]) &&($temp[3] eq 'min'))
               {
  #                print $hits1{$hit}{'query'}{$seqmin}{'rpident'}."\t".$seqmin."\t";
                  print OUT $hits1{$hit}{'query'}{$seqmin}{'rpident'}."\t".$seqmin."\t";
               }
#               print $hit."\t";
               print OUT $hit."\t";
               $min=100;
               $max=0;
               $seqmax='';
               $seqmin='';
               foreach my $query2 (keys %{$hits2{$hit}{'query'}})
               {
                  if ((!$temp[3])||($temp[3] eq 'all'))
                  {
      #               print $query2."\t".$hits2{$hit}{'query'}{$query2}{'rpident'}."\t";
                     print OUT $query2."\t".$hits2{$hit}{'query'}{$query2}{'rpident'}."\t";
                  }
                  elsif ( ($temp[3] eq 'max') &&($max<=$hits2{$query2}{'rpident'}))
                  {
                     $max=$hits2{$query2}{'rpident'};
                     $seqmax=$query2;
                  }
                  elsif ( ($temp[3] eq 'max') &&($min>=$hits2{$query2}{'rpident'}))
                  {
                     $min=$hits2{$query2}{'rpident'};
                     $seqmin=$query2;
                  }
               }
               if (($temp[3]) &&($temp[3] eq 'max'))
               {
  #                print $seqmax."\t".$hits2{$hit}{'query'}{$seqmax}{'rpident'}."\t";
                  print OUT $seqmax."\t".$hits2{$hit}{'query'}{$seqmax}{'rpident'}."\t";
               }
               elsif (($temp[3]) &&($temp[3] eq 'min'))
               {
    #              print $seqmin."\t".$hits2{$hit}{'query'}{$seqmin}{'rpident'}."\t";
                  print OUT $seqmin."\t".$hits2{$hit}{'query'}{$seqmin}{'rpident'}."\t";
               }               
    #           print "\n";
               print OUT   "\n";
            }
         }
      }
      
      elsif (($temp[0] eq '-ident2') &&($temp[2])) #dreprecated 02/23/10
      {
         my ($ress2_ref,$hress2)=&blast_results_parser($temp[1]);
         my %nhres1=&get_totals($bres_ref,$hit_ref);
         my %nhres2=&get_totals($ress2_ref,$hress2);
         &get_ident2(\%nhres1,\%nhres2,$bres_ref,$ress2_ref);
      }
      elsif (($temp[0]  eq '-ident3')&&($temp[2])&&($temp[4])&&($temp[7])&&($temp[5]))
      {
         #-ident3:taxa1:cluster_file_tx1:taxa2:blast_file_tx2:cluster_file_tx2:taxa3:cluster_file_tx3
         my ($bres2_ref,$hres2_ref)=&blast_results_parser($temp[4]);
         my %clust1=&get_clusters($temp[2]);
         my %clust2=&get_clusters($temp[5]);
         my %clust3=&get_clusters($temp[7]);
            
#            print "tx1 : ".$temp[1]."\ntx2 : ".$temp[3]."\ntx3 : ".$temp[6]."\n";
         my %master=&associate($temp[1],$bres_ref,\%clust1,$temp[3],$bres2_ref,\%clust2,$temp[6],\%clust3);
         &get_ident3(\%master,$bres_ref,$bres2_ref,$temp[1],$temp[3],$temp[6],$temp[8]);
        # }
      }
      #-bhit_list:seqs_nreads_file:hits_nreads_file:output_file   <------------ for lutzomyia assembly assertion
      elsif (($temp[0] eq '-bhit_list') && ($temp[3]))
      {
         my $out=$temp[3];
         open(OUT,">$out")||die "could not open $out due to $!\n";
         my %seqr=&get_nreads($temp[1]);
         my %hitr=&get_nreads($temp[2]);
         my $times=0; # n different length
         my $times2=0;  #n different nreads
         my $times3=0;  #n mine > paper
         my $times4=0; # n paper > mine
         my $times5=0; #n nreads mine >paper
         my $times6=0; #n nreads mine < paper
         foreach my $seq (keys %ress)
         {
            foreach my $hit (keys %{$ress{$seq}{'hits'}})
            {
              $ress{$seq}{'hits'}{$hit}{'tconserved'}=&get_tident(\%ress,$seq,$hit);
               if ($ress{$seq}{'hits'}{$hit}{'rank'}==1)
               {
#                  print OUT $seqr{$seq}{'nreads'}."\t".$seqr{$seq}{'length'}.$seq."\t".$hit."\t".$hitr{$hit}{'length'}."\t".$hitr{$hit}{'nreads'}."\n";
                     if ($seq=~ /\.f1|\.r1/)
                     {
#                        print "we have mates \n";
                        $seqr{$seq}{'length'}=$ress{$seq}{'length'};
                        $seqr{$seq}{'nreads'}=2;
                     
                     }
                     elsif ($seq=~ /\.p1|\.q1/)
                     {
       #                  print "we have singlets\n";
                        $seqr{$seq}{'length'}=$ress{$seq}{'length'};
                        $seqr{$seq}{'nreads'}=1;
#                          print "length : ".$ress{$seq}{'length'}."\n";
                     }
#                  if ($ress{$seq}{'hits'}{$hit}{'e-value'}>0)
                  if ($seqr{$seq}{'length'}!=$hitr{$hit}{'length'})
                  {
                     my $ldiff=$seqr{$seq}{'length'}-$hitr{$hit}{'length'};
                     
                     if ($ldiff>0)
                     {$times3++;}
                     else {$times4++;}
#                     print $seqr{$seq}{'nreads'}."\t".$seqr{$seq}{'length'}."\t".$seq."\t".$hit."\t".$hitr{$hit}{'length'}."\t".
  #                            $hitr{$hit}{'nreads'}."\t".$ress{$seq}{'hits'}{$hit}{'raw_sc  ore'}."\t".$ress{$seq}{'hits'}{$hit}{'e-value'}."\t".$ress{$seq}{'hits'}{$hit}{'tconserved'}."\n";                     
                     $times++;
                  }
                  if ($seqr{$seq}{'nreads'}!=$hitr{$hit}{'nreads'})
                  {
                     my $nrdiff=$seqr{$seq}{'nreads'}-$hitr{$hit}{'nreads'};
                     if ($nrdiff>0)    {$times5++;}
                     else {$times6++;}
#                     print $seqr{$seq}{'nreads'}."\t".$seqr{$seq}{'length'}."\t".$seq."\t".$hit."\t".$hitr{$hit}{'length'}."\t".
  #                            $hitr{$hit}{'nreads'}."\t".$ress{$seq}{'hits'}{$hit}{'raw_sc  ore'}."\t".$ress{$seq}{'hits'}{$hit}{'e-value'}."\t".$ress{$seq}{'hits'}{$hit}{'tconserved'}."\n";                     
                     $times2++;
                  }
#print OUT $seqr{$seq}{'nreads'}."\t".$seqr{$seq}{'length'}."\t".$seq."\t".$hit."\t".$hitr{$hit}{'length'}."\t".
  #                               $hitr{$hit}{'nreads'}."\t".$ress{$seq}{'hits'}{$hit}{'raw_score'}."\t".$ress{$seq}{'hits'}{$hit}{'e-value'}."\t".$ress{$seq}{'hits'}{$hit}{'tconserved'}."\n";
            
               }
            }
         }
         print " we had different nreads " .$times2." with ".$times5." of mine having more nreads and ".$times6." of theirs\n or diffferent length ".$times." pairs of seuqneces of which ".$times3." of mine lenger compaired to ".$times4."\n";
      }
      elsif ( ($temp[0] eq '-cluster_2nd_file')&&($temp[1] ))
      {
         my ($bres_ref2,$hit_ref2)=&blast_results_parser($temp[1]);
         &two_file_clustering($hit_ref,$hit_ref2,$bres_ref,$bres_ref2,$temp[2],$temp[3]);
      }
      elsif ( ($temp[0] eq '-bhit_pident')&&($temp[1] ))
      {
         my $out=$temp[1];
         my %idents;
         open(OUT,">$out")||die "could not open $out due to $!\n";
         foreach my $seq (keys %ress)
         {
            foreach my $hit (keys %{$ress{$seq}{'hits'}})
            {
               if ($ress{$seq}{'hits'}{$hit}{'rank'}==1)
               {
                  
#               print  OUT $ress{$seq}{'length'}."\t".$seq."\t".$hit."\t".$ress{$seq}{'hits'}{$hit}{'length'}."\t".
               #                  $ress{$seq}{'hits'}{$hit}{'raw_score'}."\t".$ress{$seq}{'hits'}{$hit}{'e-value'}."\t".$ress{$seq}{'hits'}{$hit}{'frags'}{1}{'pident'}."\n";
                  for (my $pid=10;$pid<=100;$pid+=10)
                  {
                     if ((($pid-10)<=$ress{$seq}{'hits'}{$hit}{'frags'}{1}{'pident'}) && ($pid>$ress{$seq}{'hits'}{$hit}{'frags'}{1}{'pident'}))
                     {$idents{$pid}++;}
                  }
               }
            }
         }
         foreach my $pid (sort {$a <=> $b} keys %idents)
         {
            print OUT $pid."\t".$idents{$pid}."\n";
         }
      }
      elsif (($temp[0] eq '-bhit_pident_specie')&&($temp[1]))
      {
         print " we have pident specie option\n";
 #        &bhit_species_pident($bres_ref,$temp[2],$temp[1]);
         &bhit_species_pident($bres_ref,$temp[1]);
      }
      elsif ($temp[0] eq '-split_id')
      {
         %query_ids=&split_id($ql_ref);
      }
      elsif (($temp[0] eq '-scount') && ((keys %query_ids)>=1) && ($temp[1]))
      {
         my %qhits;
         foreach my $rid (keys %query_ids)
         {
            my $hit=0;
            my $temp=$rid."_";
            foreach my $sense (keys %{$query_ids{$rid}})
            {
               $temp.=$sense."_frame_";
               foreach my $frame (keys %{$query_ids{$rid}{$sense}})
               {
                  my $ttemp=$temp.$frame;
                  #print "temp_id : ".$ttemp."\n";
                  if (exists $ress{$ttemp})
                  {
#                     print "query ".$ttemp." has ".$ress{$ttemp}{'nhits'}." hits\n";
                     $hit=1;
                  }
               }
            }
            if ($hit)
            {
               $qhits{$rid}=1;
            }
           }
         print "we have ".(keys %qhits)." queries with hits out of ".(keys %query_ids)."\n";
      }
      elsif (($temp[0] eq '-combine_go') && ($temp[1]))
      {
         my %got=&go_terms_list_b2go($temp[1]);
         &combine_gene_list($temp[2],\%query_ids,\%got);
         
      }
      elsif ($temp[0] eq '-wordfq_list')
      {
         my %freqs=&kword_fhash_builder(\%ress,$temp[1]);
         print "we have kwords to look for\n";
         foreach my $word (sort  keys %freqs)
         {
            print $word."\t".$freqs{$word}{'times'}."\t".$freqs{$word}{'eval_avrg'}."\n";
         }
      }
      elsif (($temp[0] eq '-bhit_specie_dist') && ($temp[2]))
      {
#         my $key=$temp[3];
 #        $key=s/'// ;
##         print "temp[3]=$temp[3]\n";
        &get_bhit_specie($bres_ref,$temp[2],$temp[1],$temp[3]);
      }
#      elsif (($temp[0] eq '-bhit_cluster') && ($temp[1]))
#      {
#         &bhit_cluster($hit_ref,$temp[1],$temp[2]);  
#      }
      elsif (($temp[0] eq '-bhit_cluster') && ($temp[1]))
      {
         &bhit_cluster2($hit_ref,$temp[1],$temp[2],$bres_ref,$temp[3],$temp[4]);  
      }
      elsif (($temp[0] eq '-bhit_list_clean') && ($temp[1]))
      {
         &bhit_list($bres_ref,$temp[1],$temp[2],$temp[3]);
      }
      elsif ($temp[0] eq '-map_to_bhit')
      {
         &mapp($bres_ref,$temp[1]);
      }
      elsif ($temp[0] eq '-pident_interval')
      {
         &get_int_pident($bres_ref);
      }
      elsif (($temp[0] eq '-cluster_pident') && ($temp[1]))
      {
         &cluster_pident($bres_ref,$temp[1]);
      }
      elsif (($temp[0] eq '-cluster_bhit_val') && ($temp[1]))
      {
         &cluster_bhit_validation($bres_ref,$temp[1],$temp[2],$temp[3]);
      }
      elsif ($temp[0] eq '-no_hits')
      {
         &get_no_hits($bres_ref,$temp[1],$temp[2]);  
      }
   }  
}

#a function that goes through the blast hits and counts the number
#of results where there is a hit containing keyword at an e-value of e-val or smaller
sub nresults_count
{
   my $hash_ref=$_[0];        chomp($hash_ref);
   my %bres=%{$hash_ref};
   my $kword=$_[1];           chomp($kword);       # the 'keyword' that we are trying to find
   my $eval=$_[2];                                 #the e-value of limit
   my $nres=0;
   my %genes;
   
   foreach my $query (keys %bres)
   {
      my $found=0;
      foreach my $hit (keys %{$bres{$query}{'hits'}})
      {
#         print "comparing : ".$bres{$query}{'hits'}{$hit}{'e-value'}." vs ".$eval."\n";
#         if ($bres{$query}{'hits'}{$hit}{'e-value'}<=$eval)
         if ( ( ($bres{$query}{'hits'}{$hit}{'description'}=~ /$kword/gi)||($bres{$query}{'hits'}{$hit}{'name'}=~ /$kword/gi))
             && ($bres{$query}{'hits'}{$hit}{'e-value'}<=$eval) &&(!$found))
         {
            $nres++;
            $found++;
            $genes{$query}=$bres{$query}{'hits'}{$hit}{'e-value'};
         }
      }
   }
   return $nres,\%genes;
}
# a function that takes in a blast results hash
sub get_ident
{
   my $hash_ref=$_[0];                       chomp($hash_ref);
   #my $ofile=$_[1];                             chomp($ofile);
   my %bres=%{$hash_ref};
   my %hits;
  # open(OUT,">$ofile")||die "could not open $ofile due to $!\n";   
   print " we are in get_ident_function\n";  
   foreach my $query (keys %bres)
   {
      foreach my $hit (keys %{$bres{$query}{'hits'}})
      {
         my $tqstart=0;
         my $tqend=0;
         my %idents;
         my $tidents=0;
         my $tfrags=0;
         my $hlength=0;
         my $pidents=0;
         my $pident_round;
         my $thlengt;
         my %overl;
         my %nolp;
         
         foreach my $frag (keys %{$bres{$query}{'hits'}{$hit}{'frags'}})
         {
            
            foreach my $frag2 (keys %{$bres{$query}{'hits'}{$hit}{'frags'}})
            {
               if ($frag!=$frag2)
               {
                  if (($bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'qstart'}<$bres{$query}{'hits'}{$hit}{'frags'}{$frag2}{'qstart'}) &&($bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'qend'}>$bres{$query}{'hits'}{$hit}{'frags'}{$frag2}{'qend'}))
                  {
                     $overl{$frag2}=1;
#                     print $query." - ".$hit." ".$frag," | ".$frag2." overlap~!\n";
                  }
                                                                                                                                        
               }
            }
            if (!exists $overl{$frag})
            {$nolp{$frag}=1;}
            $tfrags++;
            if (($tqstart==0) && ($tqend==0))
            {
               $tqstart=$bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'qstart'};
               $tqend=$bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'qend'};
            }
            elsif ($bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'qstart'}<$tqstart)
            {
               $tqstart=$bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'qstart'};
            }
            elsif ($tqend<$bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'qend'})
            {
               $tqend=$bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'qend'};               
            }
            
            foreach my $ident (@{$bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'qident_pos'}})
            {
               if (! exists $idents{$ident})
               {
                  $idents{$ident}=1;
               }
            }
         }
     #    print "tfrag : ".$tfrags."\n";
         
         if ($tfrags==1)
         {
#            print $bres{$query}{'hits'}{$hit}{'frags'}{1}{'qlength'}." - ".$bres{$query}{'hits'}{$hit}{'frags'}{1}{'nrident'}."\n";
            $hlength=$bres{$query}{'hits'}{$hit}{'frags'}{1}{'qlength'};
            $tidents=$bres{$query}{'hits'}{$hit}{'frags'}{1}{'nrident'};
         }
         else
         {
            foreach my $frag (keys %nolp)
            {$hlength+=$bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'qlength'};            
#          $tidents=(keys %idents)
            $tidents+=$bres{$query}{'hits'}{$hit}{'frags'}{$frag}{'nrident'};
            }
         }
         if ($bres{$query}{'algorithm'}!~ /blastn/gi)
            {  $pidents=$tidents*100*3/$hlength; }
         else
            {  $pidents=$tidents*100/$hlength; }
            $pident_round=sprintf("%.2f", $pidents);
#         my $pidents=$nrident*100/$hlength;
#         print $query." - ".$hit." total qstart  - ".$tqstart.", total qend - ".$tqend." nfrags ".(keys %{$bres{$query}{'hits'}{$hit}{'frags'}}).#"\n";
#               " length ratio : ". $bres{$query}{'length'}." / ".$bres{$query}{'hits'}{$hit}{'length'}." idents : ".( (keys %idents)/3).
#               "  %idnetity  ".$pident_round."\n";
#               " tlength ".$hlength." tidents ".$tidents."\n";
#         print OUT $query."\t".$hit. "\t". $bres{$query}{'length'}."\t".$bres{$query}{'hits'}{$hit}{'length'}."\t".$pident_round."\n";
         $hits{$hit}{'query'}{$query}{'length'}=$bres{$query}{'length'};
         $hits{$hit}{'query'}{$query}{'rpident'}=$pident_round;
         $hits{$hit}{'query'}{$query}{'tqstart'}=$tqstart;
         $hits{$hit}{'query'}{$query}{'tqend'}=$tqend;
         $hits{$hit}{'query'}{$query}{'nfrags'}=(keys %{$bres{$query}{'hits'}{$hit}{'frags'}});
         $hits{$hit}{'query'}{$query}{'tlength'}=$hlength;
         $hits{$hit}{'query'}{$query}{'tidents'}=$tidents;
      }
   }
   print " we should have ".(keys %hits)." results\n";
   return %hits;
}
sub get_hits
{
   my $hash_ref=$_[0];                       chomp($hash_ref);
   my %bres1=%{$hash_ref};
   my $hash_ref2=$_[1];                       chomp($hash_ref2);
   my %bres2=%{$hash_ref2};
   my %hits;
   my $times=0;
   print "we are getting the hitz \n";
   foreach my $query1 (keys %bres1)
   {
      foreach my $hit1 (keys %{$bres1{$query1}{'hits'}})
      {
         if ($bres1{$query1}{'hits'}{$hit1}{'rank'}==1)
         {
            foreach my $query2 (keys %bres2)
            {
               if  ( $bres2{$query2}{'hits'}{$hit1}{'rank'}==1)
                  {
                     print $query1." - ".$query2." ".$hit1." rank ".$bres1{$query1}{'hits'}{$hit1}{'rank'}."\n";
                  $times++;
                  }
               
            }
         }
      }
   }
   print "times : ".$times."\n";
   return %hits;
}

sub get_ident2
{
   my $hash_ref1=$_[0];          chomp($hash_ref1);
   my $hash_ref2=$_[1];          chomp($hash_ref2);
   my $hash_ref3=$_[2];          chomp($hash_ref3);
   my $hash_ref4=$_[3];          chomp($hash_ref4);
   my %hres1=%{$hash_ref1};
   my %bres1=%{$hash_ref3};
   my %hres2=%{$hash_ref2};
   my %bres2=%{$hash_ref4};
   my $times=0;
   
   foreach my $hit (keys %hres1)
   {
      my $nqs=(keys %{$hres1{$hit}});
      my $once=0;
      if ($nqs>2)
      {
      print "for hit ".$hit." we have ".$nqs. " queries\n";
      }
      foreach my $query1 (keys %{$hres1{$hit}})
      {
      if ($nqs>2)
      {print " ".$query1."\n";}
         foreach my $query2 (keys %{$hres2{$hit}})
         {
            if ((!$once)&&($nqs>2))
            {$once++;
               print "--".$query2."\n";}
            if ( ( (($hres1{$hit}{$query1}{'thstart'}<=$hres2{$hit}{$query2}{'thstart'}) &&  ($hres1{$hit}{$query1}{'thend'}>=$hres2{$hit}{$query2}{'thstart'})) ||
                  (($hres2{$hit}{$query2}{'thstart'}<=$hres1{$hit}{$query1}{'thstart'}) &&  ($hres2{$hit}{$query2}{'thend'}>=$hres1{$hit}{$query1}{'thstart'})))&&
                ($hres1{$hit}{$query1}{'rank'}==1) && ($hres2{$hit}{$query2}{'rank'}==1) )
            {
               $times++;
               
#               print $hit." ".$query1." - ".$hres1{$hit}{$query1}{'rank'}." ; ".$query2 ." - ".$hres2{$hit}{$query2}{'rank'}." : ".$hres1{$hit}{$query1}{'thstart'}." - ".$hres1{$hit}{$query1}{'thend'} ."\n";
               foreach my $frag1 (keys %{$bres1{$query1}{'hits'}{$hit}{'frags'}})
               {
                  
               }
            }
         }
      }
   }
   print "we have ".$times."\n";
}
# a function that returns the total leght and identities of a hit
sub get_totals
{
   my $hash_ref=$_[0];        chomp($hash_ref);
   my %bres=%{$hash_ref};
   my $hash_ref2=$_[1];          chomp($hash_ref2);
   my %hres=%{$hash_ref2};
   
   foreach my $hit (keys %hres)
   {
      foreach my $query (keys %{$hres{$hit}})
      {
         my %overl;
         my %nolp;
         my $tfrags=0;
         my $tqstart=0;
         my $tqend=0;
         my $thstart=0;
         my $thend=0;
         
         foreach my $frag1 (keys %{$bres{$query}{'hits'}{$hit}{'frags'}})
         {
            foreach my $frag2 (keys %{$bres{$query}{'hits'}{$hit}{'frags'}})
            {
                if ($frag1!=$frag2)
               {
                  if (($bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'qstart'}<$bres{$query}{'hits'}{$hit}{'frags'}{$frag2}{'qstart'}) &&
                      ($bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'qend'}>$bres{$query}{'hits'}{$hit}{'frags'}{$frag2}{'qend'}))
                  {
                     $overl{$frag2}=1;
                  }
               }
            }
            if (!exists $overl{$frag1})
            {$nolp{$frag1}=1;}
            $tfrags++;
            if (($tqstart==0) && ($tqend==0))
            {
               $tqstart=$bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'qstart'};
               $tqend=$bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'qend'};
               $thstart=$bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'hstart'};
               $thend=$bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'hend'};
               
            }
            elsif ($bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'qstart'}<$tqstart)
            {
               $tqstart=$bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'qstart'};
               $thstart=$bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'hstart'};
            }
            elsif ($tqend<$bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'qend'})
            {
               $tqend=$bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'qend'};               
               $thend=$bres{$query}{'hits'}{$hit}{'frags'}{$frag1}{'hend'};               
            }
         }
         $hres{$hit}{$query}{'thstart'}=$thstart;
         $hres{$hit}{$query}{'thend'}=$thend;
      }
   }
   return %hres;
}
#a function that  gets all the clusterr information
sub get_clusters
{
   my $file=$_[0];         chomp($file);
   my @data=&file_reader($file);       chomp(@data);
   my %clusters;
   
   foreach my $line (@data)
   {
      my @temp=split/\t/,$line;
      for (my $i=1;$i<$#temp;$i++)
      {
         $clusters{$temp[0]}{$temp[$i]}=1;
      }
   }
   return %clusters;
}
#a function that calculates the identity of sequences in a cluster
sub associate
{
   my $tx1=$_[0];                   chomp($tx1);               #the first taxa, the blast results hash & cluster file
   my $hash_ref1=$_[1];        chomp($hash_ref1);
   my %bres1=%{$hash_ref1};
   my $hash_ref2=$_[2];         chomp($hash_ref2);
   my %clust1=%{$hash_ref2};
   my $tx2=$_[3];                     chomp($tx2);              #the second taxa, the blast results hash & cluster file
   my $hash_ref3=$_[4];          chomp($hash_ref3);
   my %bres2=%{$hash_ref3};
   my $hash_ref4=$_[5];          chomp($hash_ref4);
   my %clust2=%{$hash_ref4};
   my $tx3=$_[6];                       chomp($tx3);              #the third taxa - taxa of reference
   my $hash_ref5=$_[7];          chomp($hash_ref5);
   my %clust3=%{$hash_ref5};
   my %master;
   
  foreach my $cluster (keys %clust3)               #going through each cluster in the reference sequence hash
  {
      foreach my $hit (keys %{$clust3{$cluster}})
      {
         if (! exists $master{$cluster}{$tx3}{$hit})
         {
            $master{$cluster}{$tx3}{$hit}=1;
         }
         foreach my $seq1 (keys %bres1)
         {
            if (exists $bres1{$seq1}{'hits'}{$hit})
            {
               if (! exists $master{$cluster}{$tx1}{$seq1})
               {
                  $master{$cluster}{$tx1}{$seq1}=1;                  
               }
            }
         }
         foreach my $seq2 (keys %bres2)
         {
            if (exists $bres2{$seq2}{'hits'}{$hit})
            {
               if (! exists $master{$cluster}{$tx2}{$seq2})
               {
                  $master{$cluster}{$tx2}{$seq2}=1;                  
               }
            }
         }
      }
  }
   return %master;
}
#a function that goes through the master hash
sub get_ident3
{
   my $hash_ref1=$_[0];             chomp($hash_ref1);         #the master hash
   my %master=%{$hash_ref1};
   my $hash_ref2=$_[1];             chomp($hash_ref2);         #the first blast results hash
   my %bres1=%{$hash_ref2};
   my $hash_ref3=$_[2];             chomp($hash_ref3);         #the second blast results hash
   my %bres2=%{$hash_ref3};
   my $tx1=$_[3];                       chomp($tx1);
   my $tx2=$_[4];                       chomp($tx2);
   my $tx3=$_[5];                       chomp($tx3);                    #reference taxa
   my $minlen=$_[6];
   my $singles=0;
#  print "tx1 : ".$tx1."\ntx2 : ".$tx2."\ntx3 : ".$tx3."\n";
#   print "keys : ".(keys %master)."\n";
      my $zero=0;
#   my $cluster='ORTHOMCL1861';
   foreach my $cluster (keys %master)
   {
      my @vals1=();
      my @vals2=();
      my $maxpident1=0;
      my $maxpident2=0;
      my $maxseq1='';
      my $maxseq2='';
      my $maxhit='';
      foreach my $hit (keys %{$master{$cluster}{$tx3}})
      {
        # ($seqm1,$hitm1,$minp1,$maxp1,$pident1,$algo1,$arr_ref1)=&get_max_ident(\%master,$cluster,$tx1,$hit,\%bres1,$minlen);
         @vals1=&get_max_ident(\%master,$cluster,$tx1,$hit,\%bres1,$minlen);
        # ($seqm2,$hitm2,$minp2,$maxp2,$pident2,$algo2,$arr_ref2)=&get_max_ident(\%master,$cluster,$tx2,$hit,\%bres2,$minlen);
         @vals2=&get_max_ident(\%master,$cluster,$tx2,$hit,\%bres2,$minlen);
      
         if ( ($vals1[0]) && ($vals2[0]))   #we have sequences hitting to the same hit for both taxa of interest
         {
            my $minp=0;
            my $maxp=0;
            my $length=0;
            if ( ($vals1[2]<=$vals2[2]) && ($vals2[2]<=$vals1[3]))         #the first sequence hit start earlier but they still overlap
               {$minp=$vals2[2];}
            elsif ( ($vals1[2]>$vals2[2]) && ($vals1[2]<=$vals2[3]))       #the second taxa hits early but still overlap;
               {$minp=$vals1[2];}
            if (($vals1[3]<=$vals2[3])&&($vals2[2]<=$vals1[3]))            #the first seq terminates early but some overlap occurs
               {$maxp=$vals1[3];}
            elsif (($vals1[3]>$vals2[3])&&($vals1[2]<=$vals2[3]))            #the first seq terminates early but some overlap occurs
               {$maxp=$vals2[3];}
            $length=$maxp-$minp+1;
            if ($length>=$minlen)
            {
               my %hidents1=&array_to_hash($vals1[6]);
               my %hidents2=&array_to_hash($vals2[6]);
               my $nidents1=0;
               my $nidents2=0;
               foreach my $ident (keys %hidents1)              #getting the number of idents in new interval for taxa1
              {
                  if (($ident>=$minp) && ($ident<=$maxp))
                     {$nidents1++;}
               }
               foreach my $ident (keys %hidents2)                 #getting the number of idents in new interval for taxa2
               {
                  if (($ident>=$minp) && ($ident<=$maxp))
                     {$nidents2++;}
               }
               my $pident1=$nidents1*100/$length;
               if ($pident1>100)
               {print $vals1[0]." - ".$vals1[2]." pident :".$pident1."=".$nidents1."*100/".$length."\n";}
               my $pident2=$nidents2*100/$length;
               if ($pident2>100)
               {print $vals2[0]." - ".$vals2[2]." pident :".$pident2."=".$nidents2."*100/".$length."\n";}
               if ( ($pident1>$maxpident1) && ($pident2>$maxpident2))
               {
                  $maxpident1=$pident1;
                  $maxpident2=$pident2;
                  $maxseq1=$vals1[0];
                  $maxseq2=$vals2[0];
                  $maxhit=$hit;
               }
            }
      }
      else {$zero++;}
      }
#     if (($maxpident1<1)||($maxpident1>100)||($maxpident2<1)||($maxpident2>100))
      if (($maxpident1)&&($maxpident2))
      {
         my $round1=sprintf("%.2f", $maxpident1);
         my $round2=sprintf("%.2f", $maxpident2);
         print $round1."\t".$maxseq1."\t".$maxhit."\t".$maxseq2."\t".$round2."\n";
      }
   }
   
   #print "we have ".$zero." clusters missing one sequence hit\n";
   
}
#an auxilary function
sub get_max_ident
{
   my $hash_ref1=$_[0];       chomp($hash_ref1); #master hash
   my %master=%{$hash_ref1};
   my $cluster=$_[1];          chomp($cluster);
   my $tx1=$_[2];                chomp($tx1);
   my $hit=$_[3];                chomp($hit);
   my $hash_ref2=$_[4];     chomp($hash_ref2);
   my %bres1=%{$hash_ref2};                              #blast results hash
   my $minlen=$_[5];
   my $minp=0;
   my $maxp=0;
   my $pident=0;
   my $hitm;
   my $seqm;
   my @identsm=();
   my $algorithm;
   
foreach my $seq1 (keys %{$master{$cluster}{$tx1}})
{
  
   if ((keys %{$bres1{$seq1}{'hits'}{$hit}{'frags'}})==1)
   {
      my $pident_round=sprintf("%.2f", $bres1{$seq1}{'hits'}{$hit}{'frags'}{1}{'pident'});
      if ($pident_round<1)
{      print "\t ".$seq1."\t".$hit."\t".$bres1{$seq1}{'hits'}{$hit}{'frags'}{1}{'hstart'}." - ".$bres1{$seq1}{'hits'}{$hit}{'frags'}{1}{'hend'}."\t"
            .$pident_round."\n";}
      if ( ($pident<$pident_round)&&($bres1{$seq1}{'hits'}{$hit}{'frags'}{1}{'qlength'}>=$minlen))
      {
         $minp=$bres1{$seq1}{'hits'}{$hit}{'frags'}{1}{'hstart'};
         $maxp=$bres1{$seq1}{'hits'}{$hit}{'frags'}{1}{'hend'};
         $pident=$pident_round;
         $hitm=$hit;
         $seqm=$seq1;
          @identsm=();
#         print " before adding : ".($#identsm+1)."\n";
         @identsm=@{$bres1{$seq1}{'hits'}{$hit}{'frags'}{1}{'hident_pos'}};
#          push @identsm,  @{$bres1{$seq1}{'hits'}{$hit}{'frags'}{1}{'hident_pos'}};
         $algorithm=$bres1{$seq1}{'algorithm'};
      }
   }
   else
   {
      my %overl;
      my %added;
      my $extra=0;
      do
      {
         $extra=0;
         foreach my $frag1 (sort  keys %{$bres1{$seq1}{'hits'}{$hit}{'frags'}})
         {
            foreach my $frag2 (sort  keys %{$bres1{$seq1}{'hits'}{$hit}{'frags'}})
            {
#               print "------------".$seq1."-".$hit." ".$frag1." ".$frag2."\n";
               if ((!$frag1)||(!$frag2))
               {print "----------------------".$seq1,"-----------------------------\n";}
                elsif (( $frag1 ne $frag2)&&
                     (! exists $overl{$frag1}) &&
                     ($bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'hlength'}>=$minlen)&&($bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag2}{'hlength'}>=$minlen))
                  {
                     if (($bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'qstart'}<$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag2}{'qstart'}) &&
                        ($bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'qend'}>$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag2}{'qend'}))
                    {
                           $overl{$frag2}{$frag1}=2;
                    }
                    elsif ( (! exists $overl{$frag2})&&
                           ((! exists $added{$frag1})&&(!exists $added{$frag2}))  &&
                           (($bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'qstart'}<$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag2}{'qstart'}) &&
                           ($bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'qend'}>$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag2}{'qstart'}) &&
                           ($bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'qend'}<$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag2}{'qend'})))
                    {
                        my $nfrag=$frag1."_".$frag2;
                        $added{$frag1}=$frag2;
                        $added{$frag2}=$frag1;
                        $overl{$frag2}{$nfrag}=2;
                        $overl{$frag1}{$nfrag}=2;                              
                        $extra=1;
                        $bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'qstart'}=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'qstart'};
                        $bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'qend'}=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag2}{'qend'};
                        $bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hstart'}=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'hstart'};
                        $bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hend'}=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag2}{'hend'};
                        $bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'qlength'}=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'qend'}-$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'qstart'};
                        $bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'length'}=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hend'}-$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hstart'};
                        foreach my $ident (@{$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'hident_pos'}})
                        {
                           my $found=0;
                           foreach my $ident2 (@{$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hident_pos'}})
                           {
                              if ($ident==$ident2)
                              {$found=1;}
                           }
                           if ((!$found)&& ($ident>=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hstart'}) && ($ident<=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hend'}))
                           {
                              push @{$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hident_pos'}},$ident;
                           }
                        }
                        foreach my $ident (@{$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag2}{'hident_pos'}})
                        {
                           my $found=0;
                           foreach my $ident2 (@{$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hident_pos'}})
                           {
                              if ($ident==$ident2)
                              {$found=1;}
                           }
                                 
                           if ( (!$found)&&($ident>=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hstart'}) && ($ident<=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hend'}))
                           {
                              push @{$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hident_pos'}},$ident;
                           }
                        }
                        my $npident=0;
                        if ($bres1{$seq1}{'algorithm'} eq 'blastn')
                        {
                           $npident=($#{$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hident_pos'}}+1)*100/$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'qlength'};
                        }
                        else
                        {
                           $npident=($#{$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'hident_pos'}}+1)*100*3/$bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'qlength'};
                        }
                        $bres1{$seq1}{'hits'}{$hit}{'frags'}{$nfrag}{'pident'}=$npident;
                    }
                  }
               }
               }
              } while ($extra==1);
               foreach my $frag1 (sort  keys %{$bres1{$seq1}{'hits'}{$hit}{'frags'}})
               {
                   my $pident_round=sprintf("%.2f", $bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'pident'});
                   if (  (! exists $overl{$frag1})&&($bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'length'}>=$minlen) &&($pident<$pident_round))
                   {
#                         print $seq1."\t".$hit."\t".$frag1."\t".$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'hstart'}." - ".$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'hend'}."\t"
#                            .$pident_round."\n";#."olap :".$overl{$frag1}."\n";
                  @identsm=();
                  @identsm=@{$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'hident_pos'}};
                  $minp=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'hstart'};
                   $maxp=$bres1{$seq1}{'hits'}{$hit}{'frags'}{$frag1}{'hend'};
                   $pident=$pident_round;
                   $hitm=$hit;
                   $seqm=$seq1;
                   $algorithm=$bres1{$seq1}{'algorithm'};
                  }
               }
   }
}
    return $seqm,$hitm,$minp,$maxp,$pident,$algorithm,\@identsm; #returning the sequence hit pair with the highest   
}
#a function that returns the number of reads for each sequence from a table file
sub get_nreads
{
   my $file=$_[0];         chomp($file);
   my @data=&file_reader($file);    chomp(@data);
   my %seqs;
   foreach my $line (@data)
   {
      my @temp=split/\t/,$line;
 #    print "id : ".$temp[0]." length : ".$temp[1]." nreads : ".$temp[2]."\n";
      $seqs{$temp[0]}{'length'}=$temp[1];
      $seqs{$temp[0]}{'nreads'}=$temp[2];
   }
   return %seqs;
}
#a function that computes the % identity per total hit length
sub get_tident
{
   my $hash_ref=$_[0];                 chomp($hash_ref);
   my $seq=$_[1];                         chomp($seq);
   my $hit=$_[2];                          chomp($hit);
   my %bres=%{$hash_ref};
   my %bfinal;
   my $pident=0;
   
   my %cons;
   foreach my $frag (keys %{$bres{$seq}{'hits'}{$hit}{'frags'}})
   {
       foreach my $pos (@{$bres{$seq}{'hits'}{$hit}{'frags'}{$frag}{'hcons_pos'}})
      {
         if (!exists $cons{$pos})
         {
            $cons{$pos}=1;
         }
      }
   }
#   print $seq." - ".$hit." we have ".(keys %cons)." conserved pos \n"; #<---- just temporarerly blocked 04.12.10
   if ($bres{$seq}{'algorithm'} eq 'blastn')
   {
#      $pident=(keys %cons)*100/3*$bres{$seq}{'length'};
      $pident=(keys %cons)*100/$bres{$seq}{'length'};
   }
   else
   {
      $pident=(keys %cons)*100/$bres{$seq}{'length'};         
   }
   return $pident;
}
# a function that clusters queries based on their best hit
sub two_file_clustering
{
   my $hash_ref1=$_[0];          chomp($hash_ref1);
   my $hash_ref2=$_[1];          chomp($hash_ref2);
   my %queryh1=%{$hash_ref1};
   my %queryh2=%{$hash_ref2};
   my $hash_ref12=$_[2];         chomp($hash_ref12);
   my $hash_ref22=$_[3];         chomp($hash_ref22);
   my %bres1=%{$hash_ref12};
   my %bres2=%{$hash_ref22};
   my $keyw=$_[4];                  chomp($keyw);
   my $out=$_[5];                    chomp($out);
   my %clusters;
   
   open(OUT,">$out")||die "could not open $out due to $!\n";
   foreach my $hit (keys %queryh1)
   {
      foreach my $query1 (keys %{$queryh1{$hit}})
      {
         if ($queryh1{$hit}{$query1}{'rank'}==1)
         {
#            print $query1."\t";
            foreach my $query2 (keys %{$queryh2{$hit}})
            {
               if ($queryh2{$hit}{$query2}{'rank'}==1)
               {
                  foreach my $frame1 (keys %{$bres1{$query1}{'hits'}{$hit}{'frags'}})
                  {
                     
  #                print "frames: ".$bres1{$query1}{'hits'}{$hit}{'frags'}{$frame1}{'qframe'}."\t";
                  #.$queryh1{$hit}{$query1}{'rank'}." ".$query1." - ".$hit." - ".$query2." ".$queryh1{$hit}{$query1}{'rank'}."\n";
                     if (! exists $clusters{$hit}{$query1}{$bres1{$query1}{'hits'}{$hit}{'frags'}{$frame1}{'qframe'}})
                        {$clusters{$hit}{$query1}{$bres1{$query1}{'hits'}{$hit}{'frags'}{$frame1}{'qframe'}}=11;}
                  }
                  foreach my $frame2 (keys %{$bres2{$query2}{'hits'}{$hit}{'frags'}})
                  {
                     if (! exists $clusters{$hit}{$query2}{$bres2{$query2}{'hits'}{$hit}{'frags'}{$frame2}{'qframe'}})
                        {$clusters{$hit}{$query2}{$bres2{$query2}{'hits'}{$hit}{'frags'}{$frame2}{'qframe'}}=11;}
                  }
                  if (!exists $clusters{$hit}{$query1})
                     {$clusters{$hit}{$query1}=100}
                  if (!exists $clusters{$hit}{$query2})
                     {$clusters{$hit}{$query2}=100}
               }
            }
#            print"\n";
         }
      }
   }
   my $inc=0;
   foreach my $hit (keys %clusters)
   {
      $inc++;
      print OUT $keyw.$inc."\t".$hit."\t";
      foreach  my $query (keys %{$clusters{$hit}})
      {
         print OUT $query."\t";
      }
      print OUT "\n";
   }
}
# a function that splits the gene id by "_" character
sub split_id
{
   my $arr_ref=$_[0];        chomp($arr_ref);
   my @qarray=@{$arr_ref};
   my %querys;
   foreach my $id (@qarray)
   {
      my @temp=split/_/, $id;
      #print "our id has ".($#temp+1)." parts\n";
      for(my $i=0;$i<=$#temp;$i++)
      {
#         print "temp[$i]=".$temp[$i]."\n";
         $querys{$temp[0]}{$temp[1]}{$temp[$#temp]}=10;
      }
   }
   return %querys;
}
#a function that reads through a .annot file
#blast2go annot parsers
#a function that takes in an .annot file from blast2GO and returns the list of genes and GO terms associated with them
sub go_terms_list_b2go
{
   my $file=$_[0];                  chomp($file);
   my @data=&file_reader($file);    chomp(@data);
   my %b2go;
   
   for (my $i=0;$i<=$#data;$i++)
   {
      my @temp=split/\t/,$data[$i];
      if (!exists $b2go{$temp[0]}{$temp[1]})
      {
         if (exists $temp[2])
            {$b2go{$temp[0]}{$temp[1]}{'name'}=$temp[2];
      #       print $temp[0]."\t".$temp[1]."\t".$temp[2]."\n";
             }
         else
            {$b2go{$temp[0]}{$temp[1]}{'name'}=2;}      
      }
      else
      {
         print " we already have the ".$temp[1]." GO term for gene  ".$temp[0]."\n";
      }
   }
   return %b2go;
}
#a function that goes through two lists of genes and outputs either the comon names or the missing
sub combine_gene_list
{
   my $opt=$_[0];
   if (!$opt)
      {$opt ="yes";}                  chomp($opt);
   my $hash_ref1=$_[1];          chomp($hash_ref1);   #list of genes in the blast_results hts
   my $hash_ref2=$_[2];          chomp($hash_ref2);   #list of genes with GO annotation
   my %gblist=%{$hash_ref1};
   my %golist=%{$hash_ref2};
   my %listoi;
   print "we have ".(keys %gblist)." genes whith hits\n";
   print "we have ".(keys %golist)." genes whith GO terms\n";
   
   foreach my $id (keys %gblist)
   {
      my $ori=$id;
      if ($id=~ /Contig/gi)
      {
         my ($temp)=split/\./,$id;
         $id=$temp;
         print "id : ".$id."\n";
      }
      if (($opt eq "present") || ($opt eq "yes"))
      {
         if (exists $golist{$id})
         {
            $listoi{$ori}=10;
         }
      }
      elsif ($opt eq 'no')
      {
         if (!exists $golist{$id})
         {
            $listoi{$ori}=10;
         }
      }
  
   }
    print "we have ".(keys %listoi)." genes that satify the condition ".$opt."\n";
}
# a sub that computes a list of description keywords frequencies for hits with an e-value
sub kword_fhash_builder
{
   my $hash_ref=$_[0];           chomp($hash_ref);
   my %bres=%{$hash_ref};
   my $eval=$_[1];
   my %freqs;
   
   if (!$eval)
      {$eval=1e-05;}
   foreach my $query (keys %bres)
   {
#      print "query : ". $query."\n";
      foreach my $result (keys %{$bres{$query}{'hits'}})
      {
         my %temph;
  #       print "e-val : ".$bres{$query}{'hits'}{$result}{'e-value'}."\n";
         if (($bres{$query}{'hits'}{$result}{'e-value'}!=0)
             &&($bres{$query}{'hits'}{$result}{'e-value'}<=$eval))
         {
            my @words=split/,|\s|\(|\)|\.|\t|"|'/,$bres{$query}{'hits'}{$result}{'description'}; # looking at each word in the description individualy
  #          my @words;
#            push @words, $bres{$query}{'hits'}{$result}{'description'};
            foreach my $word (@words)        # each word is considered once for each hit
            {
#               print "word : ".$word."\n";
               $temph{$word}=1;  
            }
         }
         elsif ($bres{$query}{'hits'}{$result}{'eval'}==0)
         {
            my @words=split/\s|\t/,$bres{$query}{'hits'}{$result}{'description'}; # looking at each word in the description individualy
            foreach my $word (@words)        # each word is considered once for each hit
            {
               $temph{$word}=1;  
            }
         }
         foreach my $word (keys %temph)
         {
            $freqs{$word}{'teval'}+=$bres{$query}{'hits'}{$result}{'e-value'};
            $freqs{$word}{'times'}++;
         }
      }
   }
   foreach my $word (keys %freqs)
   {
      $freqs{$word}{'eval_avrg'}=$freqs{$word}{'teval'}/$freqs{$word}{'times'};  #the eval avrage for each hit word for hits with e-val less than a given eval
   }
   return %freqs;
}
#a function that tries to computes 
#a function that counts the identity of a sequence for best hit according to their species
sub bhit_species_pident
{
   my $hash_ref=$_[0];        chomp($hash_ref);
   my %bres=%{$hash_ref};
   my $out=$_[1];             chomp($out);
   my %spident;
   open(OUT,">$out")||die "could not open $out due to $!\n";

   print " we are here  :in the species ident function\n";
   foreach my $query (keys %bres)
   {
      foreach my $hit (keys %{$bres{$query}{'hits'}})
      {
         if ($bres{$query}{'hits'}{$hit}{'rank'}==1)
         {
            my @temp=split/\[|\]/,$bres{$query}{'hits'}{$hit}{'description'};
            my $specie=$temp[$#temp];
            my ($ident,$plen,$pcons)=&get_tident2($hash_ref,$query,$hit);
            $spident{$specie}{$ident}{'times'}++;
#            if ($ident>100)
#            {
               print OUT $query."\t".
               #$bres{$query}{'hits'}{$hit}{'description'}."\t".
               $specie."\t".$ident."\t".$pcons."\t".$plen."\n";
               print $query."\t".
               $specie."\t".$ident."\t".$pcons."\t".$plen."\n";
#           }
         }
      }
   }
   return %spident;
}
# aseond function
#a function that computes the % identity per total hit length # second /05/10/10
sub get_tident2
{
   my $hash_ref=$_[0];                 chomp($hash_ref);
   my $seq=$_[1];                         chomp($seq);
   my $hit=$_[2];                          chomp($hit);
   my %bres=%{$hash_ref};
   my %cons;
   my %ident;
   my $pident=0;
   my $pcons=0;
   my @cpos;
   my @ipos;
   foreach my $frag (keys %{$bres{$seq}{'hits'}{$hit}{'frags'}})
   {
       foreach my $pos (@{$bres{$seq}{'hits'}{$hit}{'frags'}{$frag}{'hident_pos'}})
      {
         if (!exists $ident{$pos})
         {
            $ident{$pos}=1;
            push @ipos, $pos;
         }
      }
       foreach my $pos (@{$bres{$seq}{'hits'}{$hit}{'frags'}{$frag}{'hcons_pos'}})
      {
         if (!exists $cons{$pos})
         {
            $cons{$pos}=1;
            push @cpos, $pos;
         }
      }
   }
   @ipos=sort {$a <=> $b} @ipos;
   @cpos=sort {$a <=> $b} @cpos;
   my $length=abs($cpos[0]-$cpos[$#cpos]);
#   my $plength=$length*100/$bres{$seq}{'length'};
#   my $plength=$bres{$seq}{'hits'}{$hit}{'length'}*100/$bres{$seq}{'length'};
   my $plength=$length*100/$bres{$seq}{'length'};
   
    if ($bres{$seq}{'algorithm'} eq 'blastx')
      {$pcons=(keys %cons)*100*3/$length;
       $pident=(keys %ident)*100*3/$length;
       }
   else 
      {$pcons=(keys %cons)*100/$length;
      $pident=(keys %ident)*100/$length;
      }
   if ($plength>100)
      {
         print $seq,"\t"
 #        .$pident.
         .$plength
          ."\t".$cpos[0]." - ".$cpos[$#cpos]."\t".($#cpos+1)."\n";
         }
return $pident,$plength,$pcons;
}
#another function to calculate the total % identity of a alignment, regardless of the length of either query or hit sequence
#  input : hash_ref, seq, hit
# output : pident
sub get_tident4
{
   my $hash_ref=$_[0];                    chomp($hash_ref);
   my $seq=$_[1];                         chomp($seq);
   my $hit=$_[2];                         chomp($hit);
   my %bres=%{$hash_ref};
   my %ident;
   my $minpos=1e+200;
   my $maxpos=-1;
   my $pident=0;
   my $tlength=0;
   
   foreach my $frag (keys %{$bres{$seq}{'hits'}{$hit}{'frags'}})
   {
       foreach my $pos (@{$bres{$seq}{'hits'}{$hit}{'frags'}{$frag}{'hident_pos'}}) # am counting only the unique positions on the hit, not the query this would eliminate tandem duplications errors 
      {
         if (!exists $ident{$pos})
         {
            $ident{$pos}=1;
            if ($pos<=$minpos)   {$minpos=$pos;}
            if ($pos>=$maxpos)   {$maxpos=$pos;}    
         }
      }
   }
   $tlength=abs($maxpos-$minpos);
   $pident=(keys %ident)*100/$tlength;
  print "tlength= ".$tlength." | pident= ".$pident."\n";
   print $seq."\t".$hit."\ttlength= ".$tlength." | tident= ".(keys %ident)." tcons ".($#{$bres{$seq}{'hits'}{$hit}{'frags'}{1}{'hident_pos'}}+1). "\n";
   #as of 03/25/2011 it seems that bioperl might no longer make distinction between identical adn conserved positions 
   return $pident;
}
 #a function that returns a list of the number of sequences having a specie as best hit for a given e-value
 sub get_bhit_specie
{
      my $hash_ref=$_[0];        chomp($hash_ref);
   my %bres=%{$hash_ref};
   my $out=$_[1];             chomp($out);
   my $eval=$_[2];
   my $spoi=$_[3];      if ($spoi) {chomp($spoi);} #specie of interest only 
   my %spoil;
   
   my %species;
   open(OUT,">$out")||die "could not open $out due to $!\n";
   foreach my $query (keys %bres)
   {
      foreach my $hit (keys %{$bres{$query}{'hits'}})
      {
         if (($eval)&&($bres{$query}{'hits'}{$hit}{'e-value'}<=$eval)&&($bres{$query}{'hits'}{$hit}{'rank'}==1))
         {
#            print $query,"\t".$hit."\t".$bres{$query}{'hits'}{$hit}{'e-value'}." / ".$eval."\n";
            my @temp=split/\[|\]/,$bres{$query}{'hits'}{$hit}{'description'};
            my $specie=$temp[$#temp];
            $species{$specie}++;
            my ($ident,$plen,$pcons)=&get_tident2($hash_ref,$query,$hit);
            if ( ($spoi)&&($specie=~/$spoi/gi))
            {
               $spoil{$query}{'pident'}=$ident;
               $spoil{$query}{'plength'}=$plen;
               $spoil{$query}{'pcons'}=$pcons;
           }
            
         }
         elsif ((!$eval)&&($bres{$query}{'hits'}{$hit}{'rank'}==1))
         {
#           print "no ".$query,"\t".$hit."\t".$bres{$query}{'hits'}{$hit}{'e-value'}."\n";
            my @temp=split/\[|\]/,$bres{$query}{'hits'}{$hit}{'description'};
            my $specie=$temp[$#temp];
            $species{$specie}++;
         }
      }
   }
   
   foreach my $specie (sort {$species {$a} <=> $species {$b}} keys %species)
   {
    #  print $specie."\t",$species{$specie}."\n";
      if (!$spoi)
         {print OUT $specie."\t",$species{$specie}."\n";}
      else
         {
            foreach my $seq (keys %spoil)
            {print OUT$seq."\t".$spoil{$seq}{'pident'}."\t".$spoil{$seq}{'pcons'}."\t".$spoil{$seq}{'plength'}."\n";}
         }
   }
 }
 
# a function that takes in a blast file and outputs the query sequences based on their common blast hit
sub bhit_cluster
{
   my $hash_ref=$_[0];        chomp($hash_ref);       #the blast results hash depending on their blast hit
   my %bhits=%{$hash_ref};
   my $out=$_[1];                   chomp($out);
   my $eval=$_[2];
   open(OUT,">$out")||die "could not open $out due to $!\n";

   foreach my $hit (keys %bhits)
   {
      print $hit."\n";
      foreach my $query (keys %{$bhits{$hit}})
      {
         if ($bhits{$hit}{$query}{'rank'}==1)
         {
            if (($eval)&&($bhits{$hit}{$query}{'eval'}<=$eval))
            {
               print "\t".$query."\n";
               print OUT $hit."\t".$query."\t".$bhits{$hit}{$query}{'description'}."\n";
            }
            elsif (!$eval)
            {
               print "\t".$query."\n";
               print OUT $hit."\t".$query."\t".$bhits{$hit}{$query}{'description'}."\n";
            }
         }
      }
   }
}
#a function that clusters the hits
sub bhit_cluster2
{
   my $hash_ref=$_[0];        chomp($hash_ref);       #the blast results hash depending on their blast hit
   my %bhits=%{$hash_ref};
   my $out=$_[1];                   chomp($out);
   my $eval=$_[2];
   my $hash_ref2=$_[3];       chomp($hash_ref2);
   my %bres=%{$hash_ref2};
   my $ns=$_[4];
   my $lfile=$_[5];           if ($lfile) {chomp($lfile);}
   my %clusters;
   my %list;
   
   if (!$ns)
      {$ns=1;}
   foreach my $hit (keys %bhits)
   {
      foreach my $query (keys %{$bhits{$hit}})
      {
         if (($lfile)&&(exists $list{$query}))
         {
            my @lista=&file_reader($lfile);
            %list=&array_to_hash(\@lista);
 #           print" we should a list file\n"
            ;
            if ($bhits{$hit}{$query}{'rank'}==1)
            {
               if (($eval)&&($bhits{$hit}{$query}{'eval'}<=$eval))
               {
                  $clusters{$hit}{'seqs'}{$query}=1;
                  $clusters{$hit}{'nseqs'}++;
                  
               }
               elsif (!$eval)
               {
                  $clusters{$hit}{'seqs'}{$query}=1;
                  $clusters{$hit}{'nseqs'}++;
               }
            }
         }
         elsif (!$lfile)
         {
#            print "no list given~\n";
            if ($bhits{$hit}{$query}{'rank'}==1)
            {
               if (($eval)&&($bhits{$hit}{$query}{'eval'}<=$eval))
               {
                  $clusters{$hit}{'seqs'}{$query}=1;
                  $clusters{$hit}{'nseqs'}++;
                  
               }
               elsif (!$eval)
               {
                  $clusters{$hit}{'seqs'}{$query}=1;
                  $clusters{$hit}{'nseqs'}++;
               }
            }
         }
         
      }
   }
   foreach my $hit (keys %clusters)
   {
      if (($ns) && ($clusters{$hit}{'nseqs'}>=$ns))
      {
         print $hit."\n";
         foreach my $seq (keys %{$clusters{$hit}{'seqs'}})
         {
            print "\t".$seq."\n";
            my $filename=$out."_".$hit.".txt";
            open(OUT,">>$filename")||die "could not open $filename due to $!\n";
            print OUT $seq."\n";
            print $seq."\n";
         }
      }
   }
}
#a function that prints out the best hits list in a cleaned fashion
sub bhit_list
{
   my $hash_ref=$_[0];        chomp($hash_ref);
   my %bres=%{$hash_ref};
   my $out=$_[1];             chomp($out);
   my $eval=$_[2];
   my $opt=$_[3];             chomp($out);
   
   if ($eval)
      {print "we care about e-val smaller than ".$eval."\n";}
   open (OUT,">$out")|| die "could not open $out for writing due to  : $!\n";
   foreach my $query (keys %bres)
   {
      foreach my $hit (keys %{$bres{$query}{'hits'}})
      {
         if ($bres{$query}{'hits'}{$hit}{'rank'}==1)
         {
            my @temp=split/\[|\]/,$bres{$query}{'hits'}{$hit}{'description'};
            my $species=$temp[$#temp];
            my $seqtype=$temp[0];
            #obtaining species reduced name 
            my $spred;
            my @temp2=split/ /,$species;
            foreach my $name (@temp2)
            {
               my @temp3=split//,$name;
               $spred.=$temp3[0];
            }
#            print "sequence type : ".$seqtype."\n";
#            print " the full species name : ".$species." - reduced to  : ".$spred."\n";
#            print "  score : ".$bres{$query}{'hits'}{$hit}{'raw_score'}."\n";
            if (!$opt)
            {
               if (($eval)&&($eval>=$bres{$query}{'hits'}{$hit}{'e-value'}))
                   {print OUT $query."\t".$seqtype."\t".$bres{$query}{'hits'}{$hit}{'raw_score'}."\t".$bres{$query}{'hits'}{$hit}{'e-value'}."\t".$hit."\t".$spred."\t".$bres{$query}{'length'}." / ".$bres{$query}{'hits'}{$hit}{'length'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{'1'}{'pident'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{'1'}{'qstart'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{'1'}{'hstart'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{'1'}{'qend'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{'1'}{'hend'}."\n";}
               elsif (!$eval)
                   {print OUT $query."\t".$seqtype."\t".$bres{$query}{'hits'}{$hit}{'raw_score'}."\t".$bres{$query}{'hits'}{$hit}{'e-value'}."\t".$hit."\t".$spred."\t".$bres{$query}{'length'}." / ".$bres{$query}{'hits'}{$hit}{'length'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{'1'}{'pident'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{'1'}{'qstart'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{'1'}{'hstart'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{'1'}{'qend'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{'1'}{'hend'}."\n";}
            }
            elsif (($opt)&&($opt eq 'id_only'))
            {
               if (($eval)&&($eval>=$bres{$query}{'hits'}{$hit}{'e-value'}))
                   {print OUT $query."\t"."\n";}
               elsif (!$eval)
                   {print OUT $query."\t"."\n";}
            }
         }
         else
         {next;}
      }
   }
}
# function that "mapps a query back to its hit - basically returns a list with ID, smallest start of match and largest end of match
sub mapp
{
   my $hash_ref=$_[0];        chomp($hash_ref);
   my %bres=%{$hash_ref};
   my $eval=$_[1];
   
   foreach my $query (keys %bres)
   {
      foreach my $hit (keys %{$bres{$query}{'hits'}})
      {
         if (($eval)&&($bres{$query}{'hits'}{$hit}{'e-value'}<=$eval)&&($bres{$query}{'hits'}{$hit}{'rank'}==1))
         {
            print $query."\t".$hit."\t".$bres{$query}{'hits'}{$hit}{'frags'}{0}{'hstart'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{0}{'hend'}."\n";
         }
         elsif (!$eval)
         {
            print $query."\t".$hit."\t".$bres{$query}{'hits'}{$hit}{'frags'}{0}{'hstart'}."\t".$bres{$query}{'hits'}{$hit}{'frags'}{0}{'hend'}."\n";
         }
      }
   }
   
}
# a function for outputting the identity interval
sub get_int_pident
{
   my $hash_ref=$_[0];     chomp($hash_ref);
   my %bres=%{$hash_ref};
   my $ominpi=100;
   my $omaxpi=0;
   my $shit;                  #the ID of the hit with the smallest ovelall identity in the results
   my $squery;                #the ID of the query with the smallest overall identity in the results
   my $hhit;                  #the ID of the hit with the highest overlall identity
   my $hquery;                #the ID of the query with the highest overall identity
   
   foreach my $query (keys %bres)
   {
      my $minpi=100;
      my $maxpi=0;
      foreach my $hit ( keys %{$bres{$query}{'hits'}})
      {
         if ($hit ne $query)
         {
            my $ident=&get_tident4($hash_ref,$query,$hit);
            if ($ident<=$minpi) {$minpi=$ident}
            if ($ident>=$maxpi) {$maxpi=$ident}
            if ($ident<=$ominpi)
            {
               $ominpi=$ident;
               $shit=$hit;
               $squery=$query;
            }
            if ($ident>=$omaxpi)
            {
               $omaxpi=$ident;
               $hhit=$hit;
               $hquery=$query;
            }
         }
      }
#      print $query ." ".$minpi." ".$maxpi."\n";
   }
   print "overall minimum identity ".$ominpi."\t".$squery."\t".$shit."\n";
   print "overall maximum identity ".$omaxpi."\t".$hquery."\t".$hhit."\n";
}
# a function that clusters all the sequences with an identity greater than X a 03/25/2011
#  input: hash_ref, min_ident
sub cluster_pident
{
   my $hash_ref=$_[0];                 chomp($hash_ref); #blast results hash reference
   my %hash=%{$hash_ref};
   my $mident=$_[1];                            #the minimal identity required to be placed in an cluster
   my %clusters;
   my %placed;
   my $cclust=1; #current cluster id
   
   foreach my $qseq (keys %hash)
   {
  #    print " qseq :".$qseq."\n";
      foreach my $hseq (keys %{$hash{$qseq}{'hits'}})
      {
         if ($qseq ne $hseq)
#         if (($qseq eq "Contig1538.f1")&& ($qseq ne $hseq))
         {
#            my $ident=&get_tident4($hash_ref,$qseq,$hseq);
            my $ident=$hash{$qseq}{'hits'}{$hseq}{'frags'}{1}{'fident'}*100; # only looking at the first match fragment - this was done as a fast way to get info and avoid the complication of overlapping matched fragments, however it would not feaseble for possiblity of insertion between to mathced fragments 03/25/2011
            if ((!exists $placed{$hseq})&&($ident >=$mident))
            {
            $ident=sprintf("%.2f", $ident);
               if (!exists $placed{$qseq})         # we have potentially a new cluster
               {
#                  print " current cluster ".$cclust."\t".$qseq."\t".$hseq."\n";
                  $clusters{$cclust}{'seq'}{$qseq}=$ident;
                  $clusters{$cclust}{'seq'}{$hseq}=$ident;
                  $clusters{$cclust}{'first'}=$qseq;           # the first/reference sequence of the cluster
                  $placed{$qseq}=1;
                  $placed{$hseq}=1;
                  $cclust++;
                #  print "cluster ".$cclust."\t".$qseq."\t".$hseq."\t".$ident."\n";
               }
               else
               {
                  foreach my $clust (keys %clusters)
                  {
                     if (exists $clusters{$clust}{'seq'}{$qseq})
                     {
                        $clusters{$clust}{'seq'}{$hseq}=$ident;
                        $placed{$hseq}=1;
               #         print "-present cluster ".$clust."\t".$qseq."\t".$hseq."\t".$ident."\n";
                     }
                  }
               }
            }
#           print $qseq."\t".$hseq."\t".$ident."\n";
         }
      }
   }
   #calculating the average % identity
   foreach my $cluster (keys %clusters)
   {
      my $allident=0;
      foreach my $seq (keys %{$clusters{$cluster}{'seq'}})
      {
         $allident+=$clusters{$cluster}{'seq'}{$seq};
      }
      $clusters{$cluster}{'pident_av'}=sprintf("%.2f",$allident/(keys %{$clusters{$cluster}{'seq'}}));
 #     print "cluster ".$cluster." : ".$clusters{$cluster}{'pident_av'}." ; ".$clusters{$cluster}{'first'}."\n";
      foreach my $seq (keys %{$clusters{$cluster}{'seq'}})
      {
         print $cluster."\t".$seq."\t".$clusters{$cluster}{'seq'}{$seq}."\t".$clusters{$cluster}{'pident_av'}."\t".$clusters{$cluster}{'first'}."\n";
      }
   }
   my $notp=0;
   if ((keys %hash)!=(keys %placed))
   {
#      print "sequences with closest hit % identity less than ".$mident."%\n";
     foreach my $seq (keys %hash)
     {
         if (!exists $placed{$seq})
         {
            print"\t".$seq."\n";
         $notp++;
         }
     }
   }
   print "we have ".(keys %placed )." clusters with sequences with %ident >= ".$mident." representing ".((keys %hash)-$notp)." placed sequences  and " .$notp." sequences not in clusters"."\n";
   return %clusters;
}

# a function that takes in a cluster file and for each cluster looks to see if all the sequences have the same blast best match
sub cluster_bhit_validation
{
   my $hash_ref=$_[0];              chomp($hash_ref);       #the blast results hash depending on their blast hit
   my %bhits=%{$hash_ref};
   my $clfile=$_[1];                chomp($clfile);
   my $out=$_[2];                   chomp($out);
   my $option=$_[3];                chomp($option); 
   my %clusters;
   my $nseqs=0;
   
   my @data=&file_reader($clfile);  chomp(@data);
   foreach my $line (@data)
   {
      if ($line=~ /\s|\t/gi)
       {
         my @temp=split/\s|\t/, $line;
         if (!exists $clusters{$temp[0]})
         {
            for (my $i=1;$i<=$#temp;$i++)
            {
               $clusters{$temp[0]}{$temp[$i]}=1;
               $nseqs++;
            }
         }
       }
   }
   foreach my $cluster (keys %clusters)
   {
      my $bhit;
      foreach my $seq (keys %{$clusters{$cluster}})
      {
         if ((!$option )|| ($option ne 'ignore_sense'))
         {if (exists $bhits{$seq})
            {
               print $seq." blast hit foun! \n";
            }
            else
            {
              print "NO HITS for ".$seq."\n";
            }
         }
         else
         {
            
            foreach my $seq2 (keys %bhits)
            {
               my ($temp)=split/\./, $seq2;
               if ($seq eq $temp)
               {
#                  print "cluster information found for sequence : ".$seq."\n";
#                  $clusters{$cluster}{$seq}{'bhit'}=
#                  print " we have ".(keys %{$bhits{$seq2}{'hits'}})." hits for sequence of interest  :" .$seq."\n";
                  foreach my $hit (keys %{$bhits{$seq2}{'hits'}})
                  {
                     if (!$bhit)
                        {$bhit=$hit;}
                     elsif ($bhit ne $hit)
                     {
#                        print "bhit - ".$bhit," hit - ".$hit."\n";
                        print "sequences in cluster ".$cluster." have different best hits\n";
                        next;
                     }
                  }
               }
            }
         }
      }
   }
   print "we have ".(keys %clusters)." clusters coressponding to ".$nseqs." sequences\n";
}
# a function that reads in a blast_Results and outputs some data regarding the no hits
# incomplete as of (08/30/2011)
sub get_no_hits
{
   my $hash_ref=$_[0];     chomp($hash_ref);
   my %hits=%{$hash_ref};
   my $opt=$_[1];          if ($opt) {chomp($opt);}
   my $out=$_[2];          if ($out) {chomp($out);}

   foreach my $seq (keys %hits)
   {
      print $seq ," : ".$hits{$seq}{'nhits'}."\n";
   }
# something fishy going on (08/30/2011) komodo won't show the last } in the file 

   
   
}