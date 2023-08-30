#!/usr/bin/env python3
import sys, os, gzip, argparse
# (c) 2023 Angel G. Rivera-Colon

PROG = sys.argv[0].split('/')[-1]

#
# Command line options
#
def parse_args():
    desc = '''Extract the read count and genotype information for a specific locus in the gstacks catalog.calls file.'''
    p = argparse.ArgumentParser(prog=PROG, description=desc)
    p.add_argument('--catalog', required=True, help='Path to the GSTACKS catalog.calls file')
    p.add_argument('--locus-id', required=True, type=int, help='ID of target locus to export the site.')
    p.add_argument('--outdir', required=False, default='.', help='Path to output directory')
    args = p.parse_args()
    args.outdir = args.outdir.rstrip('/')
    assert os.path.exists(args.outdir), f'Error: {args.outdir} does not exist.'
    assert os.path.exists(args.catalog), f'Error: {args.catalog} does not exist.'
    assert args.locus_id > 0, f'Error: locus-id ({args.locus_id}) must be > 0.'
    return args

def main():
    args = parse_args()
    # Prepare output
    outf = open(f'{args.outdir}/gstacks_bgc_stats.tsv', 'w')
    outf.write('snp\tref\talt\taf\ttotCnts\tsample\tgt\tgq\tlnl_MM\tlnl_Mm\tlnl_mm\tcntA\tcntC\tcntG\tcntT\n')
    samples = list()
    site_cnts = dict()
    loc_found = False
    # Parse input catalog VCF
    with gzip.open(args.catalog, 'rt') as fh:
        for i, line in enumerate(fh):
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                fields = line.strip('\n').split('\t')
                samples = fields[9:]
                continue
            fields = line.strip('\n').split('\t')
            locus = int(fields[0])
            # Skip all the unwanted loci
            if locus != args.locus_id:
                continue
            else:
                loc_found = True
            col = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            snp = f'{locus}_{col}'
            # Skip invariant sites
            if alt == '.':
                continue
            # Get the total counts and frequencies
            total_cnts = fields[7].split(';')[3].split('=')[1]
            allele_freq = float(fields[7].split(';')[2].split('=')[1])
            # Process the individual counts
            for s, sam in enumerate(fields[9:]):
                geno_fields = sam.split(':')
                # GT:PS:FT:GQ:DP:AD:GL
                genotype = geno_fields[0]
                genoqual = 0
                genoliks = [-0.00,-0.00,-0.00]
                ntcounts = [0,0,0,0]
                if genotype != './.':
                    genoqual = geno_fields[3]
                    genoliks = geno_fields[6].split(',')
                    ntcounts = geno_fields[7].split(',')
                # Adjust types
                genoliks = [ float(gt) for gt in genoliks ]
                ntcounts = [ int(nt) for nt in ntcounts ]
                row = f'{snp}\t{ref}\t{alt}\t{allele_freq}\t{total_cnts}\t{samples[s]}\t{genotype}\t{genoqual}\t{genoliks[0]}\t{genoliks[1]}\t{genoliks[2]}\t{ntcounts[0]}\t{ntcounts[1]}\t{ntcounts[2]}\t{ntcounts[3]}\n'
                outf.write(row)
                # Add to the total counts per site
                nts = ['A','C','G','T']
                site_cnts.setdefault(snp, { nt : 0 for nt in nts } )
                for j, nt in enumerate(nts):
                    nt_cnt = ntcounts[j]
                    site_cnts[snp][nt] += nt_cnt
    outf.close()
    if not loc_found:
        sys.exit(f'Error: specified locus ({args.locus_id}) was never found in the input catalog.calls file.')

# Run Code
if __name__ == '__main__':
    main()
