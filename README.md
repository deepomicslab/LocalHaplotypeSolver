# LocalHaplotypeSolver

Requires Python2.7 and PuLP package for linear programming.

```shell script
python -m pip install pulp
```

## Example

```shell script
python LocalHap.py TestCases/HPV16_LRP1B_Complex.InputForLocalHapAlgorithm.txt -o output.txt
```

## File Format
We use a custom designed text format for input and output files.

#### Input
```text
##VERSION=v1.1
##_please_find_HEADER_in_online_demo_
SAMPLE	HPV16_LRP1B_Complex
PLOIDY	2m0
SOURCE	H1
SINK	H4
FVGM	ID=FVGM1;SEGSTR=V1+,V2+,V3+
SEG	ID=H1;CN_ORIG=1.97;W=1;INTERVAL=chr2:141740000-141759732;PTDP=29.56
SEG	ID=H2;CN_ORIG=9.96;W=1;INTERVAL=chr2:141759733-141790204;PTDP=149.47
SEG	ID=H3;CN_ORIG=12.11;W=1;INTERVAL=chr2:141790205-141802259;PTDP=181.66
SEG	ID=H4;CN_ORIG=1.92;W=1;INTERVAL=chr2:141802260-141806500;PTDP=28.86
SEG	ID=V1;CN_ORIG=9.84;W=1;INTERVAL=HPV16:1-4648;RAWDP=147.54
SEG	ID=V2;CN_ORIG=0.01;W=1;INTERVAL=HPV16:4649-6531;RAWDP=0.12
SEG	ID=V3;CN_ORIG=9.78;W=1;INTERVAL=HPV16:6532-7906;RAWDP=146.65
JUNC	ID=JUNC1;SEGLINK=H1+,H3+;CN_ORIG=1.33;W=1;REFSEG_5P=chr2;BKPOS_5P=141759732;STRD_5P=+;REFSEG_3P=chr2;BKPOS_3P=141790205;STRD_3P=+;NCPR=0;RAWJRC=20;PTJRC=20
JUNC	ID=JUNC2;SEGLINK=H2-,H3+;CN_ORIG=9.27;W=1;REFSEG_5P=chr2;BKPOS_5P=141759733;STRD_5P=-;REFSEG_3P=chr2;BKPOS_3P=141790205;STRD_3P=+;NCPR=0;RAWJRC=139;PTJRC=139
JUNC	ID=JUNC3;SEGLINK=V1+,H3-;CN_ORIG=6.80;W=1;REFSEG_5P=HPV16;BKPOS_5P=4648;STRD_5P=+;REFSEG_3P=chr2;BKPOS_3P=141802259;STRD_3P=-;NCPR=0;RAWJRC=102;PTJRC=102
JUNC	ID=JUNC4;SEGLINK=H2+,V3+;CN_ORIG=9.33;W=1;REFSEG_5P=chr2;BKPOS_5P=141790204;STRD_5P=+;REFSEG_3P=HPV16;BKPOS_3P=6532;STRD_3P=+;NCPR=0;RAWJRC=140;PTJRC=140
```

`#` lines starting with `#` are comments

`SAMPLE` Name of the sample. It's for your own bookkeeping and it won't affect how this tool functions.

`PLOIDY` Ploidy of virus integrated region, formatted as `<TotalNumberOfAlleles>m<NumberOfMinorAlleles>`. 
e.g. `1m0` means only one allele. `2m0` means 2 same alleles (LOH), 
`3m1` means 2 major alleles and 1 minor allele.

`SOURCE` and `SINK` are the leftmost and rightmost segment in virus integrated region. They should be 
human segments `H` and numbered with the lowest and highest index.

`FVGM` (Optional) Free virus genome, used by the tool to detect if a ucyc can be considered as a free virus.

`SEG` Segments, named by unique `ID` tag with value prefixed by `H` and `V` followed by an integer index.
Tag `CN_ORIG` refers to the original copy number of the segment. 
Optional tag `W` is the weight you would like to assign to the segment during linear programming.

`JUNC` Junctions. Need tag `SEGLINK` to specify the source and target vertex of the junction. Notice that 
`H37+,H2+` is the same as `H2-,H37-`, you only need to specify one. 
It also need `CN_ORIG` tag to specify it's copy number.

For segments and junctions, feel free to add your own tags if you find it convenient for your work.

(Read `readGraph2()` method for more details and optional tags)

#### Output

```text
@DOC	VERSION=1.0.2;SAMPLE_ID=HPV16_LRP1B_Complex;PLOIDY=GROUP1:2M:0m;ARGS=['TestCases/HPV16_LRP1B_Complex.InputForLocalHapAlgorithm.txt']
@HSEG	ID=H1;GROUP=1;CN_ORIG=1.97;CN_REF=0.00;CN_ALT=1.97;CN_ALT_LP=2.00;LP_WEIGHT=1.00;LP_OFFSET%=1.52%;NOTES=SOURCE;W=1;INTERVAL=chr2:141740000-141759732;PTDP=29.56
@HSEG	ID=H2;GROUP=1;CN_ORIG=9.96;CN_REF=0.00;CN_ALT=9.96;CN_ALT_LP=10.00;LP_WEIGHT=1.00;LP_OFFSET%=0.40%;W=1;INTERVAL=chr2:141759733-141790204;PTDP=149.47
@HSEG	ID=H3;GROUP=1;CN_ORIG=12.11;CN_REF=0.00;CN_ALT=12.11;CN_ALT_LP=12.00;LP_WEIGHT=1.00;LP_OFFSET%=-0.91%;W=1;INTERVAL=chr2:141790205-141802259;PTDP=181.66
@HSEG	ID=H4;GROUP=1;CN_ORIG=1.92;CN_REF=0.00;CN_ALT=1.92;CN_ALT_LP=2.00;LP_WEIGHT=1.00;LP_OFFSET%=4.17%;NOTES=SINK;W=1;INTERVAL=chr2:141802260-141806500;PTDP=28.86
@VSEG	ID=V1;CN_ORIG=9.84;CN_ALT=9.84;CN_ALT_LP=10.00;LP_WEIGHT=1.00;LP_OFFSET%=1.63%;W=1;INTERVAL=HPV16:1-4648;RAWDP=147.54
@VSEG	ID=V2;CN_ORIG=0.01;CN_ALT=0.00;CN_ALT_LP=N/A;LP_WEIGHT=N/A;LP_OFFSET%=N/A;W=1;INTERVAL=HPV16:4649-6531;RAWDP=0.12
@VSEG	ID=V3;CN_ORIG=9.78;CN_ALT=9.78;CN_ALT_LP=10.00;LP_WEIGHT=1.00;LP_OFFSET%=2.25%;W=1;INTERVAL=HPV16:6532-7906;RAWDP=146.65
@JUNC	ID=JUNC1;SEGLINK=H1+,H3+;CN_ORIG=1.33;CN_ALT=1.33;CN_ALT_LP=2.00;LP_WEIGHT=1.00;LP_OFFSET%=50.38%;W=1;REFSEG_5P=chr2;BKPOS_5P=141759732;STRD_5P=+;REFSEG_3P=chr2;BKPOS_3P=141790205;STRD_3P=+;NCPR=0;RAWJRC=20;PTJRC=20
@JUNC	ID=JUNC2;SEGLINK=H2-,H3+;CN_ORIG=9.27;CN_ALT=9.27;CN_ALT_LP=10.00;LP_WEIGHT=1.00;LP_OFFSET%=7.87%;W=1;REFSEG_5P=chr2;BKPOS_5P=141759733;STRD_5P=-;REFSEG_3P=chr2;BKPOS_3P=141790205;STRD_3P=+;NCPR=0;RAWJRC=139;PTJRC=139
@JUNC	ID=JUNC3;SEGLINK=V1+,H3-;CN_ORIG=6.80;CN_ALT=6.80;CN_ALT_LP=10.00;LP_WEIGHT=1.00;LP_OFFSET%=47.06%;W=1;REFSEG_5P=HPV16;BKPOS_5P=4648;STRD_5P=+;REFSEG_3P=chr2;BKPOS_3P=141802259;STRD_3P=-;NCPR=0;RAWJRC=102;PTJRC=102
@JUNC	ID=JUNC4;SEGLINK=H2+,V3+;CN_ORIG=9.33;CN_ALT=9.33;CN_ALT_LP=10.00;LP_WEIGHT=1.00;LP_OFFSET%=7.18%;W=1;REFSEG_5P=chr2;BKPOS_5P=141790204;STRD_5P=+;REFSEG_3P=HPV16;BKPOS_3P=6532;STRD_3P=+;NCPR=0;RAWJRC=140;PTJRC=140
@UCYC	ID=UCYC1;SEGSTR=H1+,H3+,H4+,H1+
@UCYC	ID=UCYC2;SEGSTR=H3+,V1-,V3-,H2-,H3+
SOLUTION	NO=1;GROUP=1;SOLUTION_UCYC=UCYC1:2,UCYC2:10;ALLELE_NO=1;ALLELE_CN=2;ALLELE_UCYC=UCYC1:1,UCYC2:5;ALLELE_FVGM_UCYC=N/A
```

`@HSEG` and `VSEG` 
* `CN_REF` and `CN_ALT` refer to the copy number in original copy number in reference allele and the mutated allele.
* `CN_ALT_LP` is the copy number mutated allele is assigned after LP
* `LP_WEIGHT` and `LP_OFFSET` are the weight of the segment in LP and the percentage of 
how much CN has changed because of LP.

`SOLUTION`
* `NO` solution number, as there can be many solutions
* `SOLUTION_UCYC` the unit cycles used in the solution, `UCYC2:2` means UCYC2 is used twice in this solution
* `ALLELE_NO` the allele that this line is talking about. If there are more than one alleles in this sample, then
this solution will occupy multiple lines with different ALLELE_NO
* `ALLELE_CN` How many copies this allele has
* `ALLELE_UCYC` the unit cycles used in this allele. Notice this is for this allele specifically, 
adding all `ALLELE_UCYC` together you get `SOLUTION_UCYC`
* `ALLELE_FVGM_UCYC` how many free virus unit cycles are there.

## Glossary

We call `H1` or `V4`, etc. as "Segments", a segment contain two vertices, e.g. `V4+` and `V4-` are the 
two vertices of segment `V4`. Similarly, "Junctions" are made up of two "Edges", e.g. 
`H37+,H2+` and `H2-,H37-` belong to the same Junction. Basically, "segment" and "junction" refer to the 
biological/physical concept of a DNA segment and the connection of two DNA segments. 
"vertex" and "edge" are used as concepts that the algorithm is based on.
 

## How does it work
The general steps that this tool takes include

1. Create the Graph data structure from user input
1. Try to add normal junctions, like `H3+=>H4+` or `V3-=>V2-`, 
since upstream software usually only report abnormal junctions 
(structural variations and virus integrations, etc.)
1. Enforce reachability, i.e. make sure every segment is reachable
from source and can reach sink.
1. Adjust segment and junction weight using Integer Linear Programming
1. Traverse graph and get all unit cycles


## Usage

```
usage: LocalHap.py [-h] [-o OUTPUT_FILE] [--host_seg_weight_factor]
                   [--virus_seg_weight_factor] [--junction_weight_factor]
                   [--cap_sv_weight] [--cap_integ_weight]
                   [--inferred_jun_weight] [--cap_jun_weight] [--auto_weight]
                   [--weight_by_length] [--weight_by_cn] [--no_infer_normal]
                   [--keep_ref_allele] [--lp_with_original_ploidy]
                   [--equally_assign_cycles_to_alleles] [--drop_ref_allele]
                   [-e] [--linear_virus] [-s] [-v] [-m] [-S] [-L] [--write_lp]
                   input_file

Reconstructing the local haplotype surrounding virus integrated regions
on human genome.

Version 1.0.2 (13-SEP-2020)

positional arguments:
  input_file                  Input file. Format requirements can be found
                              on our GitHub page.

optional arguments:
  -h, --help                  show this help message and exit
  -o OUTPUT_FILE              Output result to file. Visit our GitHub page
                              for explanation of output formats.
                              [<input_file>.ucyc.txt]
  --host_seg_weight_factor    Weight of host segments. Linear programming
                              coefficient of all host segments are scaled
                              with this factor. 1 to disable scaling. [1.0]
  --virus_seg_weight_factor   Weight of virus segments. Linear programming
                              coefficient of all virus segments are scaled
                              with this factor. 1 to disable scaling. [1.0]
  --junction_weight_factor    Weight of junctions. Linear programming
                              coefficient of all junctions are scaled
                              with this factor. 1 to disable scaling. [1.0]
  --cap_sv_weight             Depending on the detection rate of upstream
                              software, you can cap the LP weight of
                              junctions of structural variations on same
                              reference. 1 to disable capping. [1.0]
  --cap_integ_weight          Depending on the detection rate of upstream
                              software, you can cap the LP weight of
                              junctions supporting virus integrations.
                              1 to disable capping. [1.0]
  --inferred_jun_weight       Set the LP weight of inferred junctions. [0.01]
  --cap_jun_weight            Cap junctions' weight to reflect the effect
                              of sequencing coverage fluctuations. [1.0]
  --auto_weight               If --auto_weight option is set,
                              --weight_by_length and --weight_by_cn
                              will become effective. All "W" tags in input
                              file will be overridden. At the same time,
                              you may use other parameters to further tune
                              weights, when multiple weighing rules are
                              imposed, the program takes the lowest one.
                              [False]
  --weight_by_length          Assign different LP weight for segments of
                              different length. Longer segments should have
                              larger LP weight since the calculation of
                              their weight is more reliable. This option
                              also affects the weighing of junctions
                              connecting segments of varying lengths.
                              format:
                               --weight_by_length 0-50:0.2,51-100:0.4,
                               101-150:0.6,151-350:0.8,351-:1
  --weight_by_cn              Assign different LP weight for segments or
                              junctions of different copy number. Segments
                              or junctions with high copy number are more
                              likely to be changed by LP.
                              format:
                               --weight_by_cn 0-100:1,101-200:0.9,201-300:
                              0.8,301-400:0.7,401-500:0.6,501-:0.5
  --no_infer_normal           By default, normal junctions are inferred
                              based on copy number of existing junctions
                              and surrounding segments. Use this option
                              to disable the feature. [False]
  --keep_ref_allele           Keep the allele which is considered a reference
                              allele by default. In the mean time, LP will
                              be solved with original copy number. [False]
  --lp_with_original_ploidy   Do linear programming with original ploidy
                              setting instead of lowering it.
                              i.e. do not change 2m0 into 1m0 for better
                              LP performance. [False]
  --equally_assign_cycles_to_alleles
                              (EXPERIMENTAL)
                              When multiple alleles are found for a sample
                              (and the alleles are different), try to assign
                              unit cycles to each allele with equal number
                              of copies. Used with --drop_ref_allele and
                              --lp_with_original_ploidy [False]
  --drop_ref_allele           Drop major ("M") or minor allele ("m") as the
                              normal copy.
                              format:
                                > If more than one group is involved.
                                    --drop_normal GROUP1:M,GROUP2:m,...
                                > If only default group is involved.
                                    --drop_normal M/m
  -e                          Upon failure at solving LP problem, do not
                              try to recover by deleting incompatible
                              edges. [False]
  --linear_virus              If the virus concerned is linear in stead of
                              circular. [False]
  -s , --select_edge          Mode of edge selection [4]:
                              1.conservative (prefer normal junction first)
                              2.random (always randomly select next junction)
                              3.random_with_memory (random but stick with 
                                previous choice of next junction)
                              4.try_all (try all possible junction choice)
  -v                          Verbose mode. [False]
  -m                          Maximum number of loops while using "try_all"
                              option to select edges. [100]
  -S                          Only output simple contig composition. [False]
  -L                          Limit number of simple contigs to be output.
                              Input a number less than 1 to disable. [-1]
  --write_lp                  Write the given LP problem to a .lp file.
```