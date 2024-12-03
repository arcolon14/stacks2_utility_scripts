# Stacks 2 utility scripts
Some utility scripts for processing and filtering data generated using the Stacks v2 software

## Sumstats to whitelist

Filter the variant sites from a POPULATIONS `populations.sumstats.tsv` file to generate a whitelist for future runs.

### Usage:

```sh
$ python3 sumstats_to_whitelist.py -h

usage: sumstats_to_whitelist.py [-h] -s SUMSTATS [-o OUTD] [-p MIN_POPS] [-r MIN_PERC_INDV]
                                [-i MIN_NUM_INDV] [-n NUMBER_SITES] [-e] [-a HWE_ALPHA]
                                [-f MIN_MAF] [-c MIN_MAC] [-t MAX_OBS_HET] [-g]

Filter a SUMSTATS file from the STACKS POPULATIONS program (populations.sumstats.tsv) to generate a whitelist
containing the STACKS catalog IDs of a subset of selected SNPs. The sites in the SUMSTATS table can be filtered
to remove sites out of Hardy-Heinberg equilibrium, under a given allele frequency/count cutoff, and under a
given number of specified populations/ samples. The user can specify the number of sites to be retained after
filtering. The final whitelist can export one SNP per locus, equivalent to the 'write-random-snp' option in
POPULATIONS.

options:
  -h, --help            show this help message and exit
  -s SUMSTATS, --sumstats SUMSTATS
                        (str) Path to the populations sumstats TSV file
  -o OUTD, --outd OUTD  (str) Path to output directory [default=./]
  -p MIN_POPS, --min-pops MIN_POPS
                        (int) Minimum number of populations required to retain a site [default=1]
  -r MIN_PERC_INDV, --min-perc-indv MIN_PERC_INDV
                        (float) Minimum percentage of individuals in a population required to retain a site
                        [default=None]
  -i MIN_NUM_INDV, --min-num-indv MIN_NUM_INDV
                        (int) Minimum number of individuals in a population required to retain a site. Cannot
                        be used alongside 'min-perc-indv' [default=None]
  -n NUMBER_SITES, --number-sites NUMBER_SITES
                        (int) Number of sites to export [default=1000]
  -e, --hwe             (bool) Retain only sites in in HWE, applied per population [default=False]
  -a HWE_ALPHA, --hwe-alpha HWE_ALPHA
                        (float) p-value cutoff to clasify a site as out of HWE [default=0.05]
  -f MIN_MAF, --min-maf MIN_MAF
                        (float) Minimum allele frequency cutoff to retain a site, applied per-population.
                        Cannot be used alongside 'min-mac' [default=None]
  -c MIN_MAC, --min-mac MIN_MAC
                        (int) Minimum allele count cutoff to retain a site, applied per-population. Cannot be
                        used alongside 'min-mac'. [default=None]
  -t MAX_OBS_HET, --max-obs-het MAX_OBS_HET
                        (float) Max observed heterozygosity cutoff to retain a site, applied per-population
                        [default=1.0, no cutoff]
  -g, --write-random-snp
                        (bool) Export only one random SNP per locus [default=False]
```

### Example

```sh
$ python3 sumstats_to_whitelist.py \
  --sumstats populations.sumstats.tsv \  # Path to POPULATIONS SUMSTATS file
  --min-pops 5 \                         # Min 5 populations to keep a site
  --min-perc-indv 0.8 \                  # Min 80% samples per pop to keep a site
  --min-maf 0.1 \                        # Min minor allele freq of 10% to keep a site
  --max-obs-het 0.8 \                    # Max observed heterozygosity of 80% to keep a site
  --number-sites 2000 \                  # Export 2000 sites
  --write-random-snp                     # Export only a single site per locus
```

### Output

Output is a STACKS whitelist file, as specified in the [documentation](https://catchenlab.life.illinois.edu/stacks/manual/#wl). File has two columns: `<locus_id><tab><column>`

```sh
8<tab>148
80    14
92    321
113   15
195   89
199   10
200   137
203   81
204   255
270   84
```

## Calculate per-individual heterozygosity

Calculate heterozygosity per-individuals including variant and invariant sites using the `populations.all.vcf` from POPULATIONS.

### Usage

```sh
$ python3 calculate_idnv_obs_het.py -h
usage: calculate_idnv_obs_het.py [-h] [-v VCF] [-o OUTDIR]

Supply a populationns.all.vcf for a set of individuals and calculate the observed
heterozygosity, i.e., the proportion of heterozygous sites across all the
genotyped sites.

options:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     (str) Path to populations.all.vcf
  -o OUTDIR, --outdir OUTDIR
                        (str) Path to output directory
```

### Example

```sh
$ python3 calculate_idnv_obs_het.py \
    --vcf populations.all.vcf \   # Path to all-sites VCF
    --outdir het_output/          # Path to output
```

### Output

```sh
#indv_id   num_total_sites  num_total_var_sites  num_sites_in_indv  num_var_sites_in_indv  num_hets_in_indv  prop_hets_total_indv_sites  prop_hets_indv_var_sites
sample_01  6366             19                   4104               19                     7                 0.00170565                  0.36842105
sample_02  6366             19                   3319               15                     0                 0.00000000                  0.00000000
sample_03  6366             19                   3358               12                     1                 0.00029780                  0.08333333
sample_04  6366             19                   3529               19                     0                 0.00000000                  0.00000000
sample_05  6366             19                   4295               19                     3                 0.00069849                  0.15789474
sample_06  6366             19                   3892               19                     0                 0.00000000                  0.00000000
sample_07  6366             19                   4135               19                     8                 0.00193470                  0.42105263
sample_08  6366             19                   3944               19                     7                 0.00177485                  0.36842105
sample_09  6366             19                   3261               16                     4                 0.00122662                  0.25000000
```

## Extract GSTACKS catalog counts

For a given locus in the GSTACKS `catalog.calls` file, extract counts of read and nucleotide per-sites, as well as genotype model outputs. This script was desiged for troubleshooting and development (e.g., exploring the behavior of the genotyping model). The program likely has limited application for other uses.

### Usage

```sh
$ python3 extract_catalog_ng_cnts.py -h
usage: extract_catalog_ng_cnts.py [-h] --catalog CATALOG --locus-id LOCUS_ID [--outdir OUTDIR]

Extract the read count and genotype information for a specific locus in the gstacks catalog.calls file.

options:
  -h, --help           show this help message and exit
  --catalog CATALOG    Path to the GSTACKS catalog.calls file
  --locus-id LOCUS_ID  ID of target locus to export the site.
  --outdir OUTDIR      Path to output directory
```

### Example

```sh
$ python3 extract_catalog_ng_cnts.py \
    --catalog catalog.calls \           # Path to catalog file
    --outdir . \                        # Path to output directory
    --locus-id 11235                    # ID of locus of interest        
```

### Output

Example of output table:

```sh
snp       ref  alt  af       totCnts     sample     gt   gq  lnl_MM  lnl_Mm  lnl_mm   cntA  cntC  cntG  cntT
11235_42  T    G    0.00495  0,0,6,1205  sample_01  0/0  40  -0.0    -10.62  -133.33  0     0     0     30
11235_42  T    G    0.00495  0,0,6,1205  sample_02  0/0  40  -0.0    -10.02  -124.56  0     0     0     28
11235_42  T    G    0.00495  0,0,6,1205  sample_03  0/0  40  -0.0    -6.71   -76.33   0     0     0     17
11235_42  T    G    0.00495  0,0,6,1205  sample_04  0/0  40  -0.0    -5.81   -63.18   0     0     0     14
11235_42  T    G    0.00495  0,0,6,1205  sample_05  0/0  40  -0.0    -7.31   -85.1    0     0     0     19
11235_42  T    G    0.00495  0,0,6,1205  sample_06  0/0  40  -0.0    -10.62  -133.33  0     0     0     30
11235_42  T    G    0.00495  0,0,6,1205  sample_07  0/0  40  -0.0    -4.6    -45.64   0     0     0     10
11235_42  T    G    0.00495  0,0,6,1205  sample_08  0/0  40  -0.0    -5.81   -63.18   0     0     0     14
11235_42  T    G    0.00495  0,0,6,1205  sample_09  0/0  40  -0.0    -8.21   -98.25   0     0     0     22
```

## Calculate catalog coverage

Calculate the per-site coverage distribution from a GSTACKS catalog.calls file. This script was desiged for troubleshooting and development and may have limited application for other uses.

### Usage

```sh
$ python3 calc_gstacks_coverage.py -h
usage: calc_gstacks_coverage.py [-h] -c CATALOG_CALLS [-o OUT_DIR]

Calculate the average-per site coverage distribution from a raw GSTACKS catalog.calls file

options:
  -h, --help            show this help message and exit
  -c CATALOG_CALLS, --catalog-calls CATALOG_CALLS   GSTACKS `catalog.calls` file.
  -o OUT_DIR, --out-dir OUT_DIR   Output directory.
```

### Example

```sh
$ python3 calc_gstacks_coverage.py \
    --catalog-calls catalog.calls \   # Path to GSTACKS catalog.calls file
    --out-dir .                       # Path to outpit directory
```

### Output

```sh
#pos  avg_cov    med_cov    sd_cov     n_loci  f_var     avg_miss_gen  sd_miss_geno
0     26.421443  25.950000  15.176338  45636   0.000000  0.182154      0.237802
1     26.421975  25.950000  15.176499  45636   0.000007  0.182150      0.237802
2     26.422181  25.950000  15.176585  45636   0.000003  0.182148      0.237797
3     26.420613  25.950000  15.175916  45636   0.000000  0.182152      0.237799
4     26.418829  25.950000  15.176294  45636   0.000003  0.182178      0.237808
...
570   1.049830   0.633333   1.170826   38030   0.000649  0.553804      0.255442
571   1.039204   0.616667   1.160621   37865   0.000604  0.555897      0.255079
572   1.027679   0.616667   1.149575   37729   0.000609  0.558289      0.254542
573   1.017073   0.600000   1.138907   37565   0.000561  0.560004      0.253805
574   1.005943   0.600000   1.128170   37428   0.000537  0.562528      0.253439
```

## Insert size histogram

Generate a histogram describing the insert size distribution in the STACKS catalog. Uses from the `gstacks.details.gz` file generating by using the `--details` option in GSTACKS. File writen by Nicolas Rochette.

### Usage

```sh
$ python3 insert_size_hist.py
Usage:
  {} gstacks.details.gz

Creates a histogram of insert sizes, according to the aln_matrix section. Outputs
to `gstacks.details.ins_z.tsv`.
```

### Example

```sh
$ python3 insert_size_hist.py /path/to/gstacks.details.gz
```

### Output

The output table describes:
1. `ins_z`: insert size in bp
2. `n_pairs`: number of read pairs of that size

```sh
ins_z   n_pairs
1       0
2       0
3       0
4       0
5       0
...
348     58957
349     55692
350     180962
351     45265
352     39688
```

## Calculate Tajima's D

Calculate Tajima's D from a reference-aligned RADseq dataset. It takes as input a VCF (`populations.snp.vcf`) and a FASTA (`populations.loci.fasta`) generated by POPULATIONS. It does the locus-wise Tajima's D calculation over the span of the locus only. This skips un-genotyped regions and avoids associated issues.

The script has two dependencies: 1) `numpy` (<https://numpy.org/>), used for storing certain internal data structures, and 2) `scikit allel` (<https://scikit-allel.readthedocs.io/>), which is used for the calculation of Tajima's D using the `allel.tajima_d()` function (see [documentation](https://scikit-allel.readthedocs.io/en/stable/stats/diversity.html?highlight=tajima#allel.tajima_d)). 

### Usage

```sh
$ python3 calc_tajimad_rad.py -h
usage: . [-h] -v STACKS_SNPS_VCF -l STACKS_LOCI_FASTA -o OUT_DIR

options:
  -h, --help            show this help message and exit
  -v STACKS_SNPS_VCF, --stacks-snps-vcf STACKS_SNPS_VCF
                        Stacks SNPs VCF (compatible with v2.57 or higher).
  -l STACKS_LOCI_FASTA, --stacks-loci-fasta STACKS_LOCI_FASTA
                        Stacks loci FASTA.
  -o OUT_DIR, --out-dir OUT_DIR
                        Output directory.
```

### Example

```sh
$ python3 calc_tajimad_rad.py \
    --stacks-snps.vcf ./populations.snps.vcf \      # Path to VCF
    --stacks-loci-fasta ./populations.loci.fasta \  # Path to loci FASTA
    --out-dir ./tajima_d                            # Output directory
```

### Output

The script produces a table reporting the Tajima's D for each RAD locus, and their genomic coordinates.

```sh
Locus_id  Chrom   Basepair  TajimasD
37        chr1    139113    -1.92150146
44        chr1    150806    -1.75730810
88        chr1    211014    -1.72341806
96        chr1    224072    -0.21944326
97        chr1    224075    -0.84519211
111       chr1    251699    1.24644541
118       chr1    252466    -1.57684205
119       chr1    252469    1.07285659
142       chr1    313116    -0.96078502
```

## Phasing RADseq haplotypes

Phase the haplotypes of the adjacent RADtags using the `phase_rad_loci.py` script.

### Basic usage

The script takes the default haplotype VCF (`populations.haps.vcf`) and loci FASTA (`populations.loci.fa`), phases the haplotypes with `PHASE` and generates a new VCF containin the new phased sequence and a FASTA of the merged loci. These input files *must* be generated from a reference-based STACKS catalog, as the genomic coordinates are used to determine adjacent tags.

The script filters both the input and output sequence. For the input, `phase_rad_loci.py` can remove haploptypes with more than a given number of variant sites (`--max-sites-in-hap`) as well as loci with more than a specified number of alleles (`--max-alleles-in-loc`). For the output, the script can filter individuals at a locus based on their phased probability (`--min-phase-prob`), and can remove loci with a low proportion of phased samples (`--min-number-samples`).

The main dependency of the script is the program `PHASE` ([Stephens et al. 2001](https://doi.org/10.1086/319501)), which is available from: <https://stephenslab.uchicago.edu/phase/download.html>.

### Usage

```
$ phase_rad_loci.py -h

  phase_rad_loci.py

  -h, --help            show this help message and exit
  -v STACKS_HAPS_VCF, --stacks-haps-vcf STACKS_HAPS_VCF
                        Stacks haplotypes VCF (compatible with v2.57 or
                        higher).
  -l STACKS_LOCI_FASTA, --stacks-loci-fasta STACKS_LOCI_FASTA
                        Stacks loci FASTA.
  -o OUT_DIR, --out-dir OUT_DIR
                        Output directory.
  --run-name RUN_NAME   (str) Name of current run. Defaults to datetime.
  -e RES_ENZYME, --res-enzyme RES_ENZYME
                        (str) Restriction enzyme.
  -s MAX_SITES_IN_HAP, --max-sites-in-hap MAX_SITES_IN_HAP
                        (int) Max number of sites in a haplotype.
  -a MAX_ALLELES_IN_LOC, --max-alleles-in-loc MAX_ALLELES_IN_LOC
                        (int) Max number of alleles in a locus.
  -x PHASE_EXE_PATH, --phase-exe-path PHASE_EXE_PATH
                        (str) Path to PHASE executable. Default `PHASE`.
  -r MIN_NUMBER_SAMPLES, --min-number-samples MIN_NUMBER_SAMPLES
                        (float) Minumim percentage of phased samples needed to
                        retain a locus.
  -p MIN_PHASE_PROB, --min-phase-prob MIN_PHASE_PROB
                        (float) Minimum phasing probability required to keep a
                        haplotype.
  --phase-dry-run       Run with existing PHASE output.
  --keep-single-tags    Keep haplotypes/loci for single (unpaired) tags.
  --delete-phase-outs   Delete the output files from the individual PHASE
                        runs.
```

### Example

```sh
phase_rad_loci.py
    --stacks-haps-vcf populations.haps.vcf \   # Path to POPULATIONS haplotype VCF
    --stacks-loci-fasta populations.loci.fa \  # Path to POPULATIONS loci FASTA
    --out-dir output_dir/ \                    # Path to output directory
    --res-enzyme 'sbfI' \                      # Restriction enzyme, SbfI
    --phase-exe-path /path/to/PHASE \          # Path to the PHASE executable
    --max-sites-in-hap 25 \                    # Max of 25 variant sites in the input haplotypes
    --max-alleles-in-loc 40 \                  # Max of 40 alleles seen at an input locus
    --run-name 'example_run' \                 # Name of the current run, basename of outputs    
    --min-number-samples 0.8 \                 # Keep loci with 80% of samples after phasing
    --min-phase-prob 0.9                       # Only keep haplotypes of phased with at least 90% probability                  
```

### Output

After phasing, the new ID of the phased loci is the merge of the two adjacent tags. So if in the original STACKS catalog, loci `123` and loci `124` are adjacent to one another in the genome (i.e., they share the same restriction cutsite), they will be merged into a new phased locus, `123_124`. This new name will be propagated to all the outputs.

## Author

**Angel G. Rivera-Colon**<sup>1,2</sup>  
<sup>1</sup>Institute of Ecology and Evolution, University of Oregon, Eugene, OR, USA <ariverac@uoregon.edu>
<sup>2</sup>Department of Evolution, Ecology, and Behavior, University of Illinois at Urbana-Champaign, Urbana, IL, USA <angelgr2@illinois.edu>
