#!/usr/bin/env python3
import sys, os, gzip

invcf = sys.argv[1]
outdir = sys.argv[2]

def main():
    # Prepare output
    outf = open(f'{outdir}/gstacks_bgc_stats.tsv', 'w')
    outf.write('snp\tref\talt\taf\ttotCnts\tsample\tgt\tgq\tlnl_MM\tlnl_Mm\tlnl_mm\tcntA\tcntC\tcntG\tcntT\n')
    samples = list()
    site_cnts = dict()
    # Parse input catalog VCF
    with gzip.open(invcf, 'rt') as fh:
        for i, line in enumerate(fh):
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                fields = line.strip('\n').split('\t')
                samples = fields[9:]
                continue
            fields = line.strip('\n').split('\t')
            locus = int(fields[0])
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
    for snp in site_cnts:
        cnts = site_cnts[snp]
        row = f'{snp}'
        for nt in sorted(cnts):
            cnt = cnts[nt]
            row += f'\t{cnt}'
        print(row)


# Run Code
if __name__ == '__main__':
    main()
