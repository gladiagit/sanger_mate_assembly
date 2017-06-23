# Sanger sequencing mate-pair based assembly pipeline
 A set of programs for mate-based assembly for Sanger type sequences. 
## Summary
This assembly pipeline uses several existing tools for calling ((phrap/phred),  cleaning (Lucy/Seqclean) and assembling (CAP3) of long EST reads with mate pair information. The reason for this pipeline is that, in the case of the sandfly Phlebotomus papatasi, usig cap3 in a two step assembly, where the first step deals only with mate-pairs resulted in an improved assembly.

## Dependencies:
<a href="http://www.phrap.org/phredphrapconsed.html">Phrap/Phred </a><br>
<a href="http://www.complex.iastate.edu/download/Lucy2/index.html">Lucy</a><br>
<a href="ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/seqclean/">SeqClean </a><br>
<a href="http://seq.cs.iastate.edu/">CAP3</a><br>
Bioperl <br>
Perl 5 <br>

### Runining instructions:
Each script can be run individually. 
1. cleaning_pipeline.pl
2. cap_mate_assembly2.pl
3. cap3_files_assembly.pl
4. ace_parser.pl

### Script description:
1. cleaning_pipeline.pl<br>
The script basically takes in the Lucy and SeqClean parameters and calls on the two external programs for cleaning and trimining of the initial raw sequences, ensuring that the polyA tails are clipped, vector and addapter sequences have been removed as well as regions with low quality score regions from th ends of the reads. 
```
perl cleaning_pipeline.pl input_fastafile input_quality_file vectorfile splicefile contaminantfile{untrimed vector/linker/adapter} temporary_folder path_to_lucy path_to_seqclean vectordb [lucy/seqclean params : -lucy:lucy_param:value / -sqcln:sqcln_param:value]
```
2. cap_mate_assembly2.pl<br>
The script handles the first step of the EST assembly, where the mated pairs are assbled with more loose parameters. The resulting files are also the resulting FASTA and QUAL files by separating and re-naming the resutling contigs (assembled ESTs) and singlets (un-assembled reads). 
```
perl cap3_mate_assembly2.pl input_fasta_file input_qual_file path_to_cap3 temporary_files_path output_filename_templat [cap3_specific_parameters]
```
3. cap3_files_assembly2.pl<br>
This script 

4. ace_parser.pl<br>
The script returns the number of contigs and singlets generated  during the assmebly processs. 
```
perl ace_parser.pl input_ace_file
```

