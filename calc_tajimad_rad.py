#!/usr/bin/env python3
import sys, os, gzip, argparse
import numpy as np
import allel as al

# Calculate Tajima's D based on Sckit-Allel: https://scikit-allel.readthedocs.io/en/stable/stats/diversity.html?highlight=tajima#allel.tajima_d

PROG = sys.argv[0].split('/')[0]

# ==========
# Input Data
# ==========

def parse_args():
    p = argparse.ArgumentParser(prog=PROG)
    p.add_argument('-v', '--stacks-snps-vcf', required=True, help='Stacks SNPs VCF (compatible with v2.57 or higher).')
    p.add_argument('-l', '--stacks-loci-fasta', required=True, help='Stacks loci FASTA.')
    p.add_argument('-o', '--out-dir', required=True, help='Output directory.')
    # Check input arguments
    args = p.parse_args()
    args.out_dir = args.out_dir.rstrip('/')
    assert os.path.exists(args.stacks_snps_vcf)
    assert os.path.exists(args.stacks_loci_fasta)
    if not os.path.exists(args.out_dir):
        sys.exit(f"Error: '{args.out_dir}': output directory does not exist.")
    return args

# =======
# Classes
# =======

#
# Stacks RAD locus FASTA
#
class StacksRadLocus:
    def __init__(self, locus_id, chromosome, basepair, direction, sequence):
        assert type(locus_id) is int
        assert type(basepair) is int
        assert direction in ['-', '+']
        nuc_set = {'A', 'C', 'G', 'T', 'N'}
        assert set(sequence).issubset(nuc_set)
        self.locid = locus_id
        self.chrom = chromosome
        self.seq   = sequence
        self.dir   = direction
        self.len   = len(self.seq)
        # This is the position of the start of the locus bases on stacks (always 5')
        self.bp    = basepair
    def __str__(self):
        return f'{self.locid} {self.dir} {self.len:,} {self.chrom} {self.bp:,}'

#
# Locus Genotypes List
#
class LocusGenotypes:
    def __init__(self, locus_id, genotype_array, variant_positions, direction, number):
        assert genotype_array.ndim == 3
        assert len(variant_positions) == len(genotype_array)
        assert direction in ['-', '+']
        self.id    = locus_id
        self.dir   = direction
        self.num   = number
        self.genos = genotype_array
        self.pos   = variant_positions
        # Check the orientation of the locus, reverse sites if in negative strand
        if self.dir == '-':
            self.pos   = variant_positions[::-1]
            self.genos = genotype_array[::-1]
    def __str__(self):
        pos_str_list = [ str(p) for p in self.pos ]
        pos_str = ','.join(pos_str_list)
        gen_str = ''
        for site in self.genos:
            site_str = '['
            for indvs in site:
                for allele in indvs:
                    if allele == -1:
                        site_str += '.'
                    else:
                        site_str += str(allele)
            site_str += ']'
            gen_str += site_str
        return f'{self.id}:{self.dir} {pos_str} {gen_str}'

#
# Parse `samples.fa` file
#
def parse_loci_fasta(loci_fasta):
    # Check inputs
    assert os.path.exists(loci_fasta)
    # Dictionary containing loci
    loci_dict = dict()
    header   = None
    sequence = None
    dropped = 0
    # Open fasta
    for line in open(loci_fasta):
        if line[0] == '#':
            continue
        line = line.strip('\n')
        # Check header
        if line[0] == '>':
            # Populate the locus object
            if header is not None and sequence is not None:
                # Structure of the FASTA header
                # >CLocus_(locus id) [(chromosome id), (basepair), (+/-)]
                # >CLocus_306493 [HiC_scaffold_21, 15104, +]
                header = header[1:]
                fields = header.split(' ')
                loc_id = int(fields[0].split('_')[1])
                chrom  = fields[1][1:].rstrip(',')
                start  = int(fields[2].rstrip(','))
                direc  = fields[3][0]
                stacks_locus = StacksRadLocus(loc_id, chrom, start, direc, sequence)
                loci_dict[loc_id] = stacks_locus
            header = line
        else:
            sequence = line
    # Process the last sequence
    header = header[1:]
    fields = header.split(' ')
    loc_id = int(fields[0].split('_')[1])
    chrom  = fields[1][1:].rstrip(',')
    start  = int(fields[2].rstrip(','))
    direc  = fields[3][0]
    stacks_locus = StacksRadLocus(loc_id, chrom, start, direc, sequence)
    loci_dict[loc_id] = stacks_locus
    return loci_dict


#
# Parse Stacks SNP VCF
#
def parse_stacks_snp_vcf(vcf_f):
    # Check inputs
    assert os.path.exists(vcf_f)
    # Main Output
    locus_genotypes_dict = dict()
    # Intermediary variables
    locus_id  = None
    locus_dir = None
    loc_num   = 0
    locus_variant_pos = list()
    locus_genotypes = list()
    # Open VCF
    with gzip.open(vcf_f, 'rt') if vcf_f.endswith('.gz') else open(vcf_f) as fh:
        for line in fh:
            if line[0] == '#':
                continue
            fields = line.strip('\n').split('\t')
            # VCF Columns
            #   0    CHROM
            #   1    POS
            #   2    ID
            #   3    REF
            #   4    ALT
            #   5    QUAL
            #   6    FILTER
            #   7    INFO
            #   8    FORMAT
            #   9+   SAMPLES
            curr_locus = int(fields[2].split(':')[0])
            snp_column = int(fields[2].split(':')[1])+1 # Scikit Allel wants 1-based coordinates
            direction  = fields[2].split(':')[2]
            site_genotypes = list()
            # Process individual genotypes
            for geno in fields[9:]:
                sample_genotypes = list()
                geno = geno.split(':')[0]
                for g in geno.split('/'):
                    if g.isnumeric() is True:
                        sample_genotypes.append(int(g))
                    else:
                        sample_genotypes.append(-1)
                site_genotypes.append(sample_genotypes)
            # If processing the very first locus, initialize
            if locus_id is None:
                locus_variant_pos.append(snp_column)
                locus_genotypes.append(site_genotypes)
                locus_id  = curr_locus
                locus_dir = direction
            # If the current site is part of the current locus
            elif curr_locus == locus_id:
                locus_variant_pos.append(snp_column)
                locus_genotypes.append(site_genotypes)
            # If encountering a new locus
            else:
                # Complete processing of the old locus
                geno_arr = np.array(locus_genotypes)
                locus_genotypes_dict[locus_id] = LocusGenotypes(locus_id, geno_arr, locus_variant_pos, locus_dir, loc_num)
                # Reset Variables
                locus_variant_pos = list()
                locus_genotypes = list()
                # Move to the new locus
                locus_variant_pos.append(snp_column)
                locus_genotypes.append(site_genotypes)
                locus_id  = curr_locus
                locus_dir = direction
                loc_num += 1
    # Process the last locus
    geno_arr = np.array(locus_genotypes)
    locus_genotypes_dict[locus_id] = LocusGenotypes(locus_id, geno_arr, locus_variant_pos, locus_dir, loc_num)
    # Return output
    return locus_genotypes_dict

#
# Calculate Tajima's D
#
def calc_tajimas_d(locus, genotype, min_sites=1):
    # Check inputs
    assert isinstance(locus, StacksRadLocus)
    assert isinstance(genotype, LocusGenotypes)
    # Process genotypes into Scikit Allel objects
    geno_array = al.GenotypeArray(genotype.genos)
    al_cnts = geno_array.count_alleles()
    # Calculate
    D = al.tajima_d(al_cnts, pos=genotype.pos, start=1, stop=locus.len, min_sites=min_sites)
    if D == np.nan:
        D = 0
    return D

#
# Loop over all loci and process
#
def process_loci(locus_dict, genotypes_dict, outdir='./'):
    # Check inputs
    assert os.path.exists(outdir)
    outdir = outdir.rstrip('/')
    # Make output
    out = open(f'{outdir}/locus_tajimas_D.tsv', 'w')
    out.write('Locus_id\tChrom\tBasepair\tTajimasD\n')

    for loc in sorted(genotypes_dict, key=lambda l: genotypes_dict[l].id):
        genotype = genotypes_dict[loc]
        assert isinstance(genotype, LocusGenotypes)
        locus = locus_dict.get(loc, None)
        if locus is None:
            continue
        assert isinstance(locus, StacksRadLocus)
        # Calculate D
        D = calc_tajimas_d(locus, genotype)
        # Write output file
        row = f'{loc}\t{locus.chrom}\t{locus.bp}\t{D:.08f}\n'
        out.write(row)
    out.close()

def main():
    args = parse_args()
    loci_fa = args.stacks_loci_fasta
    vcf     = args.stacks_snps_vcf
    outdir  = args.out_dir
    # Parse Locus Fasta, get locus ID, coordinates, and length
    loci = parse_loci_fasta(loci_fa)
    # Parse Stacks VCF and extract the genotypes
    genotypes = parse_stacks_snp_vcf(vcf)
    # Process all loci
    process_loci(loci, genotypes, outdir)

# Run Code
if __name__ == '__main__':
    main()
