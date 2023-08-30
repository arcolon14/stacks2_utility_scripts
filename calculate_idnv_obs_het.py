#!/usr/bin/env python3
import sys, os, gzip, argparse, datetime, statistics

# Global variables
PROG = sys.argv[0].split('/')[-1]

def parse_command_line():
    desc = '''Supply a populationns.all.vcf for a set of \
    individuals and calculate the observed heterozygosity, \
    i.e., the proportion of heterozygous sites across all \
    the genotyped sites.'''
    p = argparse.ArgumentParser(description=desc, prog=PROG)
    p.add_argument('-v', '--vcf', type=str,
                   help='(str) Path to populations.all.vcf')
    p.add_argument('-o', '--outdir', type=str,
                   required=False, default='.',
                   help='(str) Path to output directory')
    args = p.parse_args()
    if not os.path.exists(args.vcf):
        sys.exit(f'Error: \'{args.vcf}\' does not exist.')
    if not os.path.exists(args.outdir):
        sys.exit(f'Error: \'{args.outdir}\' does not exist.')
    return args

def now() -> str:
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

def time():
    '''Print today's date'''
    return f'{datetime.datetime.now().strftime("%Y%m%d")}'

def simplify_genotype_str(genotype_field: str):
    '''
    Simplify the genotype field for a sample

    Args:
        genotypes_field (str): genotype string as they appear in the VCF

    Returns:
        simplified_genotype (tuple): tuple of the alleles of the individual
    '''
    indv_genos = [None, None]
    genotype = genotype_field.split(':')[0]
    if '/' in genotype:
        alleles = genotype.split('/')
    elif '|' in genotype:
        alleles = genotype.split('|')
    for i in range(len(alleles)):
        allele = alleles[i]
        if allele.isnumeric():
            indv_genos[i] = int(allele)
    simplified_genotype = list(indv_genos)
    return simplified_genotype

def is_genotype_missing(genotype: tuple):
    '''
    Check if the genotype is missing (None) at the individual

    Args:
        genotype (tuple): length-2, contains the index of each allele or None with missing

    Returns:
        missing (bool)
    '''
    missing = False
    if None in genotype:
        missing = True
    return missing

def is_indv_heterozygous(genotype: tuple):
    '''
    Check if the individual is heterozygoys

    Args:
        genotype (tuple): length-2, contains the index of each allele or None with missing

    Returns:
        het (bool)
    '''
    missing = is_genotype_missing(genotype)
    assert missing is False, f'{genotype} {missing}'
    assert len(genotype) == 2
    if genotype[0] == genotype[1]:
        return False
    else:
        return True

def process_vcf(vcf_f: str, outdir: str):
    '''
    Parse the input VCF and tally the proportion of heterozygous sites in each individual. Save into a file.

    Args:
        vcf_f (str): Path to input VCF
        outdir (str): Path to output directory

    Returns:
        None. Saves tally into output table.
    '''
    print('Parsing VCF...',flush=True)
    # Tallies
    total_sites = 0          # Number of genotyped sites
    variant_sites = 0        # Number of variant sites
    per_indv_sites = list()  # Number of sites seen per individual
    per_indv_vars = list()   # Number of variant sites in the individual
    per_indv_hets = list()   # Number of heterozygous sites per individual
    samples = list()         # Order and ID of indviduals
    # Open VCF file handle
    infh = None
    if vcf_f.endswith('gz'):
        infh = gzip.open(vcf_f, 'rt')
    else:
        infh = open(vcf_f, 'r')
    # Parse VCF
    for line in infh:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            samples = line.strip('\n').split('\t')[9:]
            print(f'    Processing genotypes for {len(samples)} samples.\n', flush=True)
            # Initialize tallies
            per_indv_sites = [ 0 for _ in samples ]
            per_indv_vars  = [ 0 for _ in samples ]
            per_indv_hets  = [ 0 for _ in samples ]
            continue
        # Process the genotype rows
        total_sites += 1
        fields = line.strip('\n').split('\t')
        # Extract the VCF row attributes
        variant = False
        ref   = fields[3]
        alt   = fields[4]
        genos = fields[9:]
        assert len(genos) == len(samples)
        # Check if site is variant
        if '.' not in {ref, alt}:
            variant = True
            variant_sites += 1
        for i in range(len(samples)):
            genotype = genos[i]
            genotype = simplify_genotype_str(genotype)
            missing = is_genotype_missing(genotype)
            if not missing:
                per_indv_sites[i] += 1
                if variant:
                    per_indv_vars[i] += 1
                    heterozygous = is_indv_heterozygous(genotype)
                    if heterozygous:
                        per_indv_hets[i] += 1
    # Close the input file
    infh.close()
    # Prepare the output file
    outfh = open(f'{outdir}/per_individual_heterozygosity.tsv', 'w')
    header = ['#indv_id', 'num_total_sites', 'num_total_var_sites',
              'num_sites_in_indv', 'num_var_sites_in_indv', 'num_hets_in_indv',
              'prop_hets_total_indv_sites', 'prop_hets_indv_var_sites']
    header = '\t'.join(header)
    outfh.write(f'{header}\n')
    # Print the per-individual rows
    for i in range(len(samples)):
        sample = samples[i]
        indv_total_sites = per_indv_sites[i]
        indv_var_sites = per_indv_vars[i]
        indv_hets = per_indv_hets[i]
        prop_hets_total = indv_hets/indv_total_sites
        prop_hets_vars = indv_hets/indv_var_sites
        row = f'{sample}\t{total_sites}\t{variant_sites}\t{indv_total_sites}\t{indv_var_sites}\t{indv_hets}\t{prop_hets_total:.08f}\t{prop_hets_vars:.08f}\n'
        outfh.write(row)
    outfh.close()
    print(f'Read a total of {total_sites:,} sites in the VCF.', flush=True)
    print(f'    A total of {variant_sites:,} ({(variant_sites/total_sites):.03%}) were variant sites.\n')
    mean_hets = statistics.mean(per_indv_hets)
    meadian_hets = statistics.median(per_indv_hets)
    sd_hets = statistics.stdev(per_indv_hets)
    print(f'''Average heterozygous sites per individual:
    Mean:   {mean_hets:,.03f}
    Median: {meadian_hets:,.03f}
    StDev:  {sd_hets:,.03f}''')

def main():
    args = parse_command_line()
    # Start program
    print(f'{PROG} started on {now()}', flush=True)
    print(f'''
    Input VCF:  {args.vcf}
    Output Dir: {args.outdir}\n''', flush=True)
    # Process VCF
    process_vcf(args.vcf, args.outdir)
    # End program
    print(f'\n{PROG} finished on {now()}', flush=True)

if __name__ == '__main__':
    main()
