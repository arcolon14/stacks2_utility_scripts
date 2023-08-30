#!/usr/bin/env python3
import sys, os, gzip, argparse
import statistics as st
# (c) 2023 Angel G. Rivera-Colon

PROG = sys.argv[0].split('/')[-1]

#
# Command line options
#
def parse_args():
    desc = '''Calculate the average-per site coverage distribution from a raw GSTACKS catalog.calls file'''
    p = argparse.ArgumentParser(prog=PROG, description=desc)
    p.add_argument('-c', '--catalog-calls', required=True, help='GSTACKS `catalog.calls` file.')
    p.add_argument('-o', '--out-dir', required=False, default='.', help='Output directory.')
    # Check input arguments
    args = p.parse_args()
    args.out_dir = args.out_dir.rstrip('/')
    if not os.path.exists(args.out_dir):
        sys.exit(f"Error: '{args.out_dir}': output directory does not exist.")
    if not os.path.exists(args.catalog_calls):
        sys.exit(f"Error: '{args.catalog_calls}' not found.")
    return args

#
# Calculate proportion of missing data at a site
#
def calc_proportion_missing_genos(genotypes_list, n_samples):
    # Check inputs
    assert type(genotypes_list) is list
    assert type(n_samples) is int
    assert n_samples > 0
    assert len(genotypes_list) == n_samples
    # Parse genotypes and tally
    n_missing = 0
    for geno in genotypes_list:
        if geno[0] == '.':
            n_missing += 1
    # Convert to proportion
    p_missing = n_missing/n_samples
    return p_missing

#
# Parse Gstacks catalog.calls
#
def parse_catalog_calls(calls_f):
    assert os.path.exists(calls_f)
    fh = gzip.open(calls_f, 'rt')
    per_site_cov_dict = dict()
    per_site_var_dict = dict()
    per_site_gen_dict = dict()
    n_samples = None
    for i, line in enumerate(fh):
        if line[0:2] == '##':
            continue
        fields = line.strip('\n').split('\t')
        if fields[0] == '#CHROM':
            n_samples = len(fields[9:])
            continue
        loc_id = int(fields[0])
        col = int(fields[1])-1
        info = fields[7].split(';')
        total_depth = int(info[0].split('=')[1])
        avg_depth = total_depth/n_samples
        variant = False
        genotypes_list = fields[9:]
        if len(info) == 3:
            # If the info column has depth for two alleles, site is variant
            variant = True
        # Populate dictionaries
        per_site_cov_dict.setdefault(col, [])
        per_site_cov_dict[col].append(avg_depth)
        per_site_var_dict.setdefault(col, 0)
        if variant is True:
            per_site_var_dict[col] += 1
        # Calculate the proportion of missing genotypes at a given site
        prop_missing_gen = calc_proportion_missing_genos(genotypes_list, n_samples)
        per_site_gen_dict.setdefault(col, [])
        per_site_gen_dict[col].append(prop_missing_gen)
    return per_site_cov_dict, per_site_var_dict, per_site_gen_dict

#
# Summarize dictionaries and save output
#
def summarize_sites(per_site_cov_dict, per_site_var_dict, per_site_gen_dict, outdir='.'):
    assert type(per_site_cov_dict) is dict
    assert type(per_site_var_dict) is dict
    assert type(per_site_gen_dict) is dict
    assert os.path.exists(outdir)
    out = open(f'{outdir}/per_site_cov_stats.tsv','w')
    out.write('#pos\tavg_cov\tmed_cov\tsd_cov\tn_loci\tf_var\tavg_miss_gen\tsd_miss_geno\n')
    for site in sorted(per_site_cov_dict):
        # Calculate coverage statistics
        cov_list = per_site_cov_dict[site]
        avg_cov = st.mean(cov_list)
        med_cov = st.median(cov_list)
        n_loci = len(cov_list)
        sd_cov = 0
        if n_loci > 1:
            sd_cov = st.stdev(cov_list)
        # Calculate the frequency of variants at a given site
        total_var = sum(per_site_var_dict.values())
        f_var = per_site_var_dict[site]/total_var
        # Calculate genotye statistics
        genos = per_site_gen_dict.get(site, [0])
        avg_gen = st.mean(genos)
        sd_gen = 0
        if len(genos) > 1:
            sd_gen = st.stdev(genos)
        row = f'{site}\t{avg_cov:.06f}\t{med_cov:.06f}\t{sd_cov:.06f}\t{n_loci}\t{f_var:.06f}\t{avg_gen:.06f}\t{sd_gen:.06f}\n'
        out.write(row)
    out.close()

#
# Main
#
def main():
    args = parse_args()
    calls = args.catalog_calls
    outdir = args.out_dir
    per_site_cov_dict, per_site_var_dict, per_site_gen_dict = parse_catalog_calls(calls)
    summarize_sites(per_site_cov_dict, per_site_var_dict, per_site_gen_dict, outdir)

# Run Code
if __name__ == '__main__':
    main()
