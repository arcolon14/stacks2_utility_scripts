# Stacks 2 utility scripts
Some utility scripts for processing and filtering data generated using the Stacks v2 software

## Sumstats to whitelist

Filter the variant sites from a `populations` SUMSTATS file to generate a whitelist for future runs.

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
