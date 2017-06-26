# Sanger sequencing mate-pair based assembly pipeline
 A set of programs for mate-based assembly for Sanger type sequences. 
## Summary
This assembly pipeline uses several existing tools for calling ((phrap/phred),  cleaning (Lucy/Seqclean) and assembling (CAP3) of long EST reads with mate pair information. The reason for this pipeline is that, in the case of the sandfly Phlebotomus papatasi, using CAP3 in a two-step assembly, where the first step deals only with mate-pairs resulted in an improved assembly.

## Dependencies:
<a href="http://www.phrap.org/phredphrapconsed.html">Phrap/Phred </a><br>
<a href="http://www.complex.iastate.edu/download/Lucy2/index.html">Lucy</a><br>
<a href="ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/seqclean/">SeqClean </a><br>
<a href="http://seq.cs.iastate.edu/">CAP3</a><br>
Bioperl <br>
Perl 5 <br>

### Running instructions:
Each script can be run individually. 
1. cleaning_pipeline.pl
2. lucy_trim.pl
3. cap_mate_assembly2.pl
4. cap3_files_assembly.pl
5. ace_parser.pl
6. blast_parser_complete.pl

### Script description:
1. cleaning_pipeline.pl<br>
The script basically takes in the Lucy and SeqClean parameters and calls on the two external programs for cleaning and trimming of the initial raw sequences, ensuring that the polyA tails are clipped, vector and adapter sequences have been removed as well as regions with low quality score regions from the ends of the reads. 
```
perl cleaning_pipeline.pl input_fastafile input_quality_file vectorfile splicefile contaminantfile{untrimed vector/linker/adapter} temporary_folder path_to_lucy path_to_seqclean vectordb [lucy/seqclean params : -lucy:lucy_param:value / -sqcln:sqcln_param:value]
```
2. lucy_trim.pl<br>
The Lucy program identifies regions of low quality but does not trim both the quality files. This script takes in the results from Lucy and trims quality files to match the FASTA files prior to the first step of assembly. 
```
perl lucy_trim.pl fastafile  quailty_file out_filename_template
```
3. cap_mate_assembly2.pl<br>
The script handles the first step of the EST assembly, where the mated pairs are assembled with more loose parameters. The resulting files are also the resulting FASTA and QUAL files by separating and re-naming the resulting contigs (assembled ESTs) and singlets (un-assembled reads). 
```
perl cap3_mate_assembly2.pl input_fasta_file input_qual_file path_to_cap3 temporary_files_path output_filename_templat [cap3_specific_parameters]
```
4. cap3_files_assembly2.pl<br>
This script performs the last step of the assembly and assembles all the reads together regardless of their mate pair information. This assembly step is usually more stringent than the first.
```
perl cap3_files_assembly.pl template_name path_to_cap3 cap3_options
```

5. ace_parser.pl<br>
The script returns the number of contigs and singlets generated during the assembly process. 
```
perl ace_parser.pl input_ace_file
```
6. blast_paerser_complete.pl<br>
This script is not directly part of the assembly pipeline but can be used for identifying the assembled sequences by parsing the results file from a BLAST job against NCBI's NR database. The command example bellow exemplifies a couple of the parsing options of the script that can be used for automated identification of the sequences.
```
perl blast_parser_complete.pl blast_result_file\n
                              [-bhit_species_dist:[eval]:outputfile:[species]]
                              [-nhits:kword:e-val]
                              [-bhit_pident:outputfile]
```
-bhit_species_dist: - top hit for a specific species<br>
-nhits:kword:eval - top hits containing keyword in annotation and have a e-value better than eval <br>
-bhit_pident - best hit's percent identity to the sequence <br>
