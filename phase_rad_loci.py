#!/usr/bin/env python3
import sys, os, re, argparse, subprocess, datetime, statistics, gzip

#
# Globals
#
DATE = datetime.datetime.now().strftime("%Y%m%d")
PROG = sys.argv[0].split('/')[-1]
# Compatible enzymes
KNOWN_RENZ = { e[0].lower() : e[1] for e in [
    ('BamHI',   'G/GATCC'),
    ('EcoRI',   'G/AATTC'),
    ('HinDIII', 'A/AGCTT'),
    ('PstI',    'C/TGCAG'),
    ('SbfI',    'CC/TGCAGG')
    ]}

#
# Stacks RAD locus FASTA
#
class StacksRadLocus:
    def __init__(self, counter, locus_id, chromosome, basepair, direction, sequence):
        assert type(counter)  is int
        assert type(locus_id) is int
        assert type(basepair) is int
        assert direction in ['-', '+']
        nuc_set = {'A', 'C', 'G', 'T', 'N'}
        assert set(sequence).issubset(nuc_set)
        self.cnt   = counter
        self.locid = locus_id
        self.chrom = chromosome
        self.seq   = sequence
        self.dir   = direction
        # This is the position of the start of the locus bases on stacks (always 5')
        self.stacks_start_bp = basepair
        # These are the adjusted coordinates on the chromosme
        self.min_bp = self.stacks_start_bp
        self.max_bp = self.min_bp + (len(self.seq)-1)
        if self.dir == '-':
            self.max_bp = self.stacks_start_bp 
            self.min_bp = self.max_bp - (len(self.seq)-1)
        self.diff   = self.max_bp - self.min_bp
        assert len(self.seq)-1 == (self.diff), f'{len(self.seq)} {self.diff} {self.dir}'
    def __str__(self):
        return f'{self.cnt} CLocus_{self.locid}({self.dir}) Len: {len(self.seq):,}bp Chr: {self.chrom} {self.stacks_start_bp:,} ({self.min_bp:,} - {self.max_bp:,})'

#
# Restriction Enzyme
# Example: 'SbfI' 'CC/TGCAGG'
class RestrictionEnzyme:
    def __init__(self, enzyme_name):
        enzyme_name = enzyme_name.lower()
        assert enzyme_name in KNOWN_RENZ
        self.name = enzyme_name
        self.cut_pattern = KNOWN_RENZ[enzyme_name]
        self.cutsite = self.cut_pattern.replace('/', '')
        remainder_length = max([ len(part) for part in self.cut_pattern.split('/') ])
        self.remainder = self.cutsite[len(self.cutsite)-remainder_length:]
        self.removed = self.cutsite[:len(self.cutsite)-remainder_length]
        self.olap_len = 2 * len(self.remainder) - len(self.cutsite)
    def __str__(self):
        return f'{self.name} ({self.cut_pattern})'

#
# RAD Haplotype
class RadHaplotype:
    def __init__(self, locus_id, chromosome, basepair, direction, hap_list, col_list, genotype_list, ploidy=2):
        # Check input
        assert type(locus_id) is int
        assert type(basepair) is int
        assert direction in ['-', '+']
        assert type(hap_list) is list
        assert len(hap_list) > 0
        assert type(col_list) is list
        assert len(col_list) > 0
        assert len(genotype_list)%ploidy == 0
        # Values
        self.locid = locus_id
        self.chrom = chromosome
        self.bp    = basepair
        self.dir   = direction
        self.haps  = hap_list
        self.cols  = col_list
        self.genos = genotype_list
        self.n_sam = len(self.genos)//ploidy
        self.ploid = ploidy
    def __str__(self):
        hap_str = ','.join(self.haps)
        col_str = [ str(c) for c in self.cols]
        col_str = ','.join(col_str)
        gen_lst = list()
        for g in self.genos:
            if g is None:
                gen_lst.append('.')
            else:
                gen_lst.append(str(g))
        gen_str = ','.join(gen_lst)
        return f'{self.locid}:{self.dir} {self.chrom}:{self.bp} {gen_str} {hap_str} {col_str}'
    # Recontruct the VCF line
    def write_vcf_line(self, ploidy=2):
        ref_hap = self.haps[0]
        alt_hap = ','.join(self.haps[1:])
        col_str = [ str(c) for c in self.cols] # NOTE: All SNP columns are 0-based
        col_str = ','.join(col_str)
        gen_str = ''
        for i in range(len(self.genos)):
            geno = self.genos[i]
            if geno is None:
                gen_str += '.'
            else:
                gen_str += str(geno)
            if i%ploidy == 0:
                gen_str += '/'
            else:
                if i < len(self.genos)-1:
                    gen_str += '\t'
        vcf_row = f'{self.chrom}\t{self.bp}\t{self.locid}:1:{self.dir}\t{ref_hap}\t{alt_hap}\t.\tPASS\tsnp_columns={col_str}\tGT\t{gen_str}'
        return vcf_row
    # Function to extract specific genotype from 'geno' string
    def get_genotype(self, sample_i, allele_n, ploidy=2):
        assert allele_n in list(range(ploidy))
        genotype = self.genos[ (sample_i * ploidy) + allele_n ]
        return genotype

#
# Merged RAD Haplotype
class MergedRadHaplotype:
    def __init__(self, for_locus_id, rev_locus_id, chromosome, basepair, hap_list, col_list, genotype_list, ploidy=2):
        # Check input
        assert type(basepair) is int
        assert type(hap_list) is list
        assert len(hap_list) > 0
        assert type(col_list) is list
        assert len(col_list) == len(hap_list[0])
        assert len(genotype_list)%ploidy == 0
        # Set the locus ID - Gio
        self.locid = None
        locus_ids = list()
        if for_locus_id is not None:
            locus_ids.append(str(for_locus_id))
        if rev_locus_id is not None:
            locus_ids.append(str(rev_locus_id))
        if len(locus_ids) > 0:
            self.locid = '_'.join(locus_ids)
        else:
            sys.exit(f"Error: For and Rev Locus cannot both be None ({chromosome} {basepair})")
        # Set the other variables
        self.chrom = chromosome
        self.bp    = basepair
        self.haps  = hap_list
        self.cols  = col_list
        self.genos = genotype_list
        self.n_sam = len(self.genos)//ploidy
        self.ploid = ploidy
    def __str__(self):
        hap_str = ','.join(self.haps)
        col_str = [ str(c) for c in self.cols]
        col_str = ','.join(col_str)
        gen_lst = list()
        for g in self.genos:
            if g is None:
                gen_lst.append('.')
            else:
                gen_lst.append(str(g))
        gen_str = ','.join(gen_lst)
        return f'{self.locid} {self.chrom}:{self.bp} {gen_str} {hap_str} {col_str}'
    # Recontruct the VCf line
    def write_vcf_line(self, ploidy=2):
        ref_hap = self.haps[0]
        alt_hap = ','.join(self.haps[1:])
        col_str = [ str(c) for c in self.cols] # NOTE: All SNP columns are 0-based
        col_str = ','.join(col_str)
        gen_str = ''
        for i in range(len(self.genos)):
            geno = self.genos[i]
            if geno is None:
                gen_str += '.'
            else:
                gen_str += str(geno)
            if i%ploidy == 0:
                gen_str += '/'
            else:
                if i < len(self.genos)-1:
                    gen_str += '\t'
        vcf_row = f'{self.chrom}\t{self.bp}\t{self.locid}:1:+\t{ref_hap}\t{alt_hap}\t.\tPASS\tsnp_columns={col_str}\tGT\t{gen_str}'
        return vcf_row
    # Function to extract specific genotype from 'geno' string
    def get_genotype(self, sample_i, allele_n, ploidy=2):
        assert allele_n in list(range(ploidy))
        genotype = self.genos[ (sample_i * ploidy) + allele_n ]
        return genotype

#
# Loci Pair
# Stores a pair of loci as their status (overlapping/nonoverlapping)
class LociPair:
    def __init__(self, for_locus_id, rev_locus_id, paired_status=False):
        assert type(rev_locus_id) is int or rev_locus_id is None
        assert type(for_locus_id) is int or for_locus_id is None
        assert type(paired_status) is bool
        self.rev_loc = rev_locus_id
        self.for_loc = for_locus_id
        self.paired  = paired_status
    def __str__(self):
        return f'{self.for_loc}(+) {self.rev_loc}(-) {self.paired}'

# Stores the output of the PHASE `pairs` file
# Contains the phased alleles and their probability of phasing
# If thresholds are not met, it stores the values but flags it as a "bad" phase
class PhasedHaplotype:
    def __init__(self, for_locus_id, rev_locus_id, sample_id, for_allele_1, rev_allele_1, for_allele_2, rev_allele_2, phase_probability, min_prob=0.5):
        # Check inputs
        assert type(for_locus_id) is int
        assert type(rev_locus_id) is int
        assert type(for_allele_1) is int
        assert type(rev_allele_1) is int
        assert type(for_allele_2) is int
        assert type(rev_allele_2) is int
        assert type(phase_probability) is float
        # Class contents
        self.for_loc_id  = for_locus_id
        self.rev_loc_id  = rev_locus_id
        self.sample_id   = sample_id
        self.f_allele_1  = for_allele_1
        self.r_allele_1  = rev_allele_1
        self.f_allele_2  = for_allele_2
        self.r_allele_2  = rev_allele_2
        self.phase_prob  = phase_probability
        self.phase_found = True
        # If the phasing probability is low phase has not been found
        if self.phase_prob < min_prob:
            self.phase_found = False
        elif -1 in [self.f_allele_1, self.r_allele_1, self.f_allele_2, self.r_allele_2]:
            self.phase_found = False
        # Create a list of the phase allele pairs so they can be accessed together later
        self.allele_list = [(self.f_allele_1, self.r_allele_1), (self.f_allele_2, self.r_allele_2)]
    def __str__(self):
        return f'MLocus_{self.for_loc_id}_{self.rev_loc_id} {self.sample_id} {self.f_allele_1} {self.r_allele_1} {self.f_allele_2} {self.r_allele_2} {self.phase_prob:.3f} {self.phase_found}'

#
# Command Line Options
#
def parse_args():
    p = argparse.ArgumentParser(prog=PROG)
    p.add_argument('-v', '--stacks-haps-vcf',    required=True,       help='Stacks haplotypes VCF (compatible with v2.57 or higher).')
    p.add_argument('-l', '--stacks-loci-fasta',  required=True,       help='Stacks loci FASTA.')
    p.add_argument('-o', '--out-dir',            required=True,       help='Output directory.')
    p.add_argument(      '--run-name',           required=False,      default=None,    help='(str) Name of current run. Defaults to datetime.')
    p.add_argument('-e', '--res-enzyme',         required=False,      default='sbfI',  help='(str) Restriction enzyme.')
    p.add_argument('-s', '--max-sites-in-hap',   required=False,      default=25,      help='(int) Max number of sites in a haplotype.', type=int)
    p.add_argument('-a', '--max-alleles-in-loc', required=False,      default=40,      help='(int) Max number of alleles in a locus.',   type=int)
    p.add_argument('-x', '--phase-exe-path',     required=True,       default='PHASE', help='(str) Path to PHASE executable. Default `PHASE`.')
    p.add_argument('-r', '--min-number-samples', required=False,      default=0.8,     help='(float) Minumim percentage of phased samples needed to retain a locus.', type=float)
    p.add_argument('-p', '--min-phase-prob',     required=False,      default=0.75,    help='(float) Minimum phasing probability required to keep a haplotype.',      type=float)
    p.add_argument(      '--phase-dry-run',      action='store_true', default=False,   help='Run with existing PHASE output.')
    p.add_argument(      '--keep-single-tags',   action='store_true', default=False,   help='Keep haplotypes/loci for single (unpaired) tags.')
    p.add_argument(      '--delete-phase-outs',  action='store_true', default=False,   help='Delete the output files from the individual PHASE runs.')

    # Check input arguments
    args = p.parse_args()
    args.out_dir = args.out_dir.rstrip('/')
    assert os.path.exists(args.stacks_haps_vcf)
    assert os.path.exists(args.stacks_loci_fasta)
    if not os.path.exists(args.out_dir):
        sys.exit(f"Error: '{args.out_dir}': output directory does not exist.")
    if args.run_name is None:
        args.run_name = datetime.datetime.now().strftime('%Y%m%d')
    if not os.path.exists(args.phase_exe_path):
        sys.exit(f"Error: '{args.out_dir}': PHASE executable could not be found.\nTry running `which PHASE`.")
    renz_str = ','.join(sorted(KNOWN_RENZ.keys()))
    if args.res_enzyme.lower() not in KNOWN_RENZ:
        sys.exit(f"Error: '{args.res_enzyme}' not among the known enzymes.\nAvailable enzymes are: {renz_str}.")
    if args.max_sites_in_hap < 1:
        sys.exit(f"Error: Max sites in haplotype should be non-zero positive integer.")
    if args.max_alleles_in_loc < 1:
        sys.exit(f"Error: Max alleles in locus should be non-zero positive integer.")
    if args.max_alleles_in_loc >= 50:
        sys.exit(f"Error: Max alleles in locus should be < 50 for compatibility with `PHASE`.")
    if not 0.0 < args.min_number_samples <= 1.0:
        sys.exit(f"Error: Min number of samples should be 0 < num <= 1.")
    if not 0.0 < args.min_phase_prob <= 1.0:
        sys.exit(f"Error: Min phasing probability should be 0 < prob <= 1.")
    return args

#
# Function to reverse complete sequence.
# It finds complement and inverts the order
def rev_comp(sequence):
    rev = []
    for nt in sequence.upper():
        if nt == 'A':
            rev.append('T')
        elif nt == 'C':
            rev.append('G')
        elif nt == 'G':
            rev.append('C')
        elif nt == 'T':
            rev.append('A')
        elif nt not in ['A', 'C', 'G', 'T']:
            rev.append('N')
    return ''.join(rev[::-1])

#
# Parse `samples.fa` file
def parse_loci_fasta(loci_fasta, enzyme, log=False):
    # Check inputs
    assert os.path.exists(loci_fasta)
    assert isinstance(enzyme, RestrictionEnzyme)
    # Dictionary containing loci
    loci_dict = dict()
    header   = None
    sequence = None
    dropped = 0
    # Open fasta
    cnt = 0
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
                # Check for the cutsite
                if enzyme.remainder != sequence[:len(enzyme.remainder)]:
                    dropped += 1
                else:
                    stacks_locus = StacksRadLocus(cnt, loc_id, chrom, start, direc, sequence)
                    loci_dict[loc_id] = stacks_locus
                    cnt += 1
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
    # Check for the cutsite
    if enzyme.remainder != sequence[:len(enzyme.remainder)]:
        dropped += 1
    else:
        stacks_locus = StacksRadLocus(cnt, loc_id, chrom, start, direc, sequence)
        loci_dict[loc_id] = stacks_locus
        cnt += 1


    if log is True:
        print(f'Loaded {cnt:,} sequences from the Stacks loci FASTA.', flush=True)
        if dropped > 0:
            print(f'    Dropped {dropped:,} sequences because of altered cutsite sequences.', flush=True)
        print('', flush=True)
    return loci_dict

#
# Look for pairs of loci that overlap at the cutsite
# TODO: Look for all overlapping loci (???)
def find_overlapping_loci(loci_dict, renz, log=False):
    # Check inputs
    assert isinstance(renz, RestrictionEnzyme)
    assert type(loci_dict) is dict
    # The pair of loci should be composed of a negative strand and a positive
    # strand locus, with overlapped position at the cutsite.
    #                            | forward strand start
    #                          CCTGCAGGNNNNNNNNNNNNN...
    #          ...NNNNNNNNNNNNNGGACGTCC
    #          reverse strand start |
    loci_pairs_list = list()
    # Loci pair in LociPair format
    n_paired = 0
    pair = LociPair(None, None, False)
    prev_locus = None
    # Loop over loci dictionary
    # In each iteration you process the loci that came before
    for locid in sorted(loci_dict):
        curr_locus = loci_dict[locid]
        assert isinstance(curr_locus, StacksRadLocus)
        prev_stat  = 'paired'
        # Skip if the previous loci hasn't been seen yet
        if prev_locus is None:
            prev_locus = curr_locus
            continue
        # The following conditions mark the previous locus as unpaired
        #
        # If the previous loci belongs to a different chromosome
        # These shouldn't overlap anyway, but good to check
        if prev_locus.chrom != curr_locus.chrom:
            prev_stat = 'unpaired'
        # If the loci are not sequential
        # Adjacent loci have sequential IDs
        elif abs(curr_locus.locid - prev_locus.locid) != 1:
            prev_stat = 'unpaired'
        # If coming from different restriction sites
        # TODO: Overlap these???
        elif abs(prev_locus.stacks_start_bp - curr_locus.stacks_start_bp) != renz.olap_len-1:
            prev_stat = 'unpaired'
        # If coming from different restriction sites that are adjacent
        elif prev_locus.dir == '-' and curr_locus.dir == '+':
            if prev_locus.stacks_start_bp < curr_locus.stacks_start_bp:
                prev_stat = 'unpaired'
        # If it is in the same strand as the current locus
        elif prev_locus.dir == curr_locus.dir:
            prev_stat = 'unpaired'
        # Process a pair of loci
        if prev_stat == 'paired':
            pair = LociPair(prev_locus.locid, curr_locus.locid, True)
            loci_pairs_list.append(pair)
            prev_locus = None
            n_paired += 1
        else:
            # These examples are unpaired
            pair = LociPair(prev_locus.locid, None, False)
            if prev_locus.dir == '-':
                pair = LociPair(None, prev_locus.locid, False)
            loci_pairs_list.append(pair)
            # Reset to see the next locus
            prev_locus = curr_locus
    if log is True:
        print(f'Found {n_paired:,} pairs of loci originating from the same restriction enzyme cutsite.\n', flush=True)
    return loci_pairs_list

#
# Reformat the loci FASTA
# Merge loci if necessary
def reformat_loci_fasta(loci_dict, loci_pairs_list, renz, basename, outdir='.', log=False):
    assert os.path.exists(outdir)
    outdir = outdir.rstrip('/')
    assert isinstance(renz, RestrictionEnzyme)
    # Generate the output
    outfa = open(f'{outdir}/{basename}.merged_loci.fa', 'w')
    nseqs = 0
    # Loop over the loci pairs
    for pair in loci_pairs_list:
        assert isinstance(pair, LociPair)
        # If the loci are unpaired
        if pair.paired is False:
            for p in [pair.rev_loc, pair.for_loc]:
                if p is None:
                    continue
                loc = loci_dict[p]
                assert isinstance(loc, StacksRadLocus)
                # Restore the original header
                # >CLocus_(locus id) [(chromosome id), (basepair), (+/-)]
                header = f'>CLocus_{loc.locid} [{loc.chrom}, {loc.stacks_start_bp}, {loc.dir}]'
                outfa.write(f'{header}\n')
                outfa.write(f'{loc.seq}\n')
                nseqs+=1
        # Process paired loci
        else:
            rev_loc = loci_dict[pair.rev_loc]
            assert isinstance(rev_loc, StacksRadLocus)
            for_loc = loci_dict[pair.for_loc]
            assert isinstance(for_loc, StacksRadLocus)
            # Create the merged sequence
            for_seq = for_loc.seq
            rev_seq = rev_comp(rev_loc.seq[renz.olap_len:])
            mer_seq = rev_seq+for_seq
            cut_sta = len(rev_seq)-len(renz.removed)
            # Create a new FASTA header
            # >M(erged)Locus_(for_loc_id)_(rev_loc_id) [(chrom), (cutsite position), (min bp - max bp), (index of cutsite first)]
            header = f'>MLocus_{for_loc.locid}_{rev_loc.locid} [{for_loc.chrom}, {for_loc.min_bp-(len(renz.cutsite) - len(renz.remainder))}, {rev_loc.min_bp}-{for_loc.max_bp}, {cut_sta}]'
            # Check output integrity
            assert len(mer_seq)-1 == (for_loc.max_bp-rev_loc.min_bp), f'>MLocus_{for_loc.locid}_{rev_loc.locid} {len(mer_seq)} {for_loc.max_bp-rev_loc.min_bp}'
            cut_end = cut_sta+len(renz.cutsite)
            assert mer_seq[cut_sta:cut_end] == renz.cutsite, f'>MLocus_{for_loc.locid}_{rev_loc.locid} - RENZ do not match'
            # Write to file
            outfa.write(f'{header}\n')
            outfa.write(f'{mer_seq}\n')
            nseqs+=1
    outfa.close()
    if log is True:
        print(f'Generating merged RAD loci FASTA.\n    Wrote sequences for {nseqs:,} loci.\n', flush=True)

#
# Parse Stacks' haplotype VCF
#
def parse_stacks_haps_vcf(vcf, log=False):
    # Check input
    assert os.path.exists(vcf) is True
    # Outputs
    individuals = None
    locus_haplotypes_dict = dict()
    # Temp Data
    individuals = None
    hap_lens    = list()
    missing_dat = list()
    # Open VCF
    with gzip.open(vcf, 'rt') if vcf.endswith('.gz') else open(vcf) as fh:
        for line in fh:
            if line[0:2] == '##':
                continue
            fields = line.strip('\n').split('\t')
            # Process header
            if fields[0] == '#CHROM':
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
                individuals = fields[9:]
                continue
            # Process the locus-wide information
            chrom = fields[0]
            bp    = int(fields[1])
            locus = int(fields[2].split(':')[0])
            direc = fields[2].split(':')[2]
            haps  = list()
            haps.append(fields[3])
            for h in fields[4].split(','):
                haps.append(h)
            cols = list()
            for c in fields[7].split('=')[1].split(','):
                col = int(c) # NOTE: As of Stacks 2.57 SNP columns in VCF outputs are 0-based.
                cols.append(col)
            # Check the orientation of the SNP columns:
            # The columns are in decreasing orders for negative strand loci
            # This for stacks v2.57+
            if direc == '-':
                if cols != sorted(cols, reverse=True):
                    sys.exit('\nERROR: SNP columns for negative strand loci not sorted by genomic coordinates. Haplotype VCF must have been generated by Stacks v2.57 or higher.')
            # Process all the Individuals genotypes into a genotype string
            genotypes = []
            for geno in fields[9:]:
                for g in geno.split('/'):
                    if g.isnumeric() is True:
                        g = int(g)
                    else:
                        g = None
                    genotypes.append(g)
            rad_haplotype = RadHaplotype(locus, chrom, bp, direc, haps, cols, genotypes)
            locus_haplotypes_dict[locus] = rad_haplotype
    if log is True:
        print(f'Loaded {len(locus_haplotypes_dict):,} loci from the Stacks Haplotype VCF.\n', flush=True)
    return locus_haplotypes_dict, individuals

#
# Function to read haplotype, filter, and calculate statistics
#
def filter_haps_statistics(raw_locus_haplotypes_dict, max_hap_len, max_alleles, outdir='.', log=False):
    assert os.path.exists(outdir)
    outdir = outdir.rstrip('/')
    assert type(raw_locus_haplotypes_dict) is dict
    assert max_hap_len > 0, f'{max_hap_len}'
    assert 0 < max_alleles < 50, f'{max_alleles}'
    if log is True:
        print('Filtering haplotypes and printing statistics...\n', flush=True)
    # Outputs
    filt_locus_haplotypes_dict = dict()
    tsv = open(f'{outdir}/locus_haplotype_stats.tsv', 'w')
    tsv.write('#locus_id\tchrom\tbp\thap_len\tn_alleles\tstatus\n')
    # Intermediary values
    haplo_len = list()
    n_alleles = list()
    # Loop over the loaded haplotypes
    for locus in sorted(raw_locus_haplotypes_dict):
        keep = True
        haplotype = raw_locus_haplotypes_dict[locus]
        assert isinstance(haplotype, RadHaplotype)
        alles = len(haplotype.haps)
        sites = len(haplotype.cols)
        if alles > max_alleles:
            keep = False
        elif sites > max_hap_len:
            keep = False
        k = 'kept'
        if keep is False:
            k = 'discarded'
        tsv.write(f'{locus}\t{haplotype.chrom}\t{haplotype.bp}\t{sites}\t{alles}\t{k}\n')
        # Only store and process kept alleles
        if keep is True:
            n_alleles.append(alles)
            haplo_len.append(sites)
            filt_locus_haplotypes_dict[locus] = haplotype
    # Print some haplotype statistics
    if log is True:
        print(f'''Kept haplotypes for a total of {len(filt_locus_haplotypes_dict):,} loci.

    Haplotype length statistics for the kept loci:
        Mean Hap Len:       {statistics.mean(haplo_len):.03f}
        StDev Hap Len:      {statistics.stdev(haplo_len):.03f}
        Median Hap Len:     {statistics.median(haplo_len):.03f}
        Min Hap Len:        {min(haplo_len)}
        Max Hap Len:        {max(haplo_len)}

    Number of alleles statistics for the kept loci:
        Mean Num Alleles:   {statistics.mean(n_alleles):.03f}
        StDev Num Alleles:  {statistics.stdev(n_alleles):.03f}
        Median Num Alleles: {statistics.median(n_alleles):.03f}
        Min Num Alleles:    {min(n_alleles)}
        Max Num Alleles:    {max(n_alleles)}\n''', flush=True)

    return filt_locus_haplotypes_dict

#
# Take genotypes of pair of loci and return the PHASE format string
# Phase Format:
#    3                               3 samples
#    5                               5 sites
#    P 300 1313 1500 2023 5635       positions of sites (bp)
#    MSSSM                           Multi or Single alleles
#    #1                              Sample 1
#    12 1 0 1 3                      Sample 1 alleles 1
#    11 0 1 0 3                      Sample 1 alleles 2
#    #2                              Sample 2
#    12 1 1 1 2                      Sample 2 alleles 1
#    12 0 0 0 3                      Sample 2 alleles 1
#    #3                              Sample 3
#    -1 ? 0 0 2                      Sample 3 alleles 1
#    -1 ? 1 1 13                     Sample 3 alleles 1
def genotypes_to_phase(for_haplotype, rev_haplotype):
    # Check inputs
    assert isinstance(rev_haplotype, RadHaplotype)
    assert isinstance(for_haplotype, RadHaplotype)
    # PHASE format components
    n_sam = rev_haplotype.n_sam
    sites = 2
    allele_type = 'M'*sites
    pos = f'P {for_haplotype.bp} {rev_haplotype.bp}'
    phase_str = f'{n_sam}\n{sites}\n{pos}\n{allele_type}\n'
    # Loop over samples and add genotype lines
    for sam_i in range(n_sam):
        phase_str += f'#{sam_i+1}\n'
        # Loop over the alleles
        for a in range(for_haplotype.ploid):
            for_geno = for_haplotype.get_genotype(sam_i, a)
            if for_geno is None:
                for_geno = -1
            rev_geno = rev_haplotype.get_genotype(sam_i, a)
            if rev_geno is None:
                rev_geno = -1
            phase_str += f'{for_geno} {rev_geno}\n'
    return phase_str

#
# Function to run a command and save its output.
# Function by Niraj Rayamajhi
def run_command(command, log_file):
    log = open(log_file, 'w')
    assert type(command) == list
    for word in command:
        assert type(word) == str
    log.write('This is log file content of PHASE\n')
    command_str = ' '.join(command)
    log.write(f'{command_str}\n\n')
    # Run the command.
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        encoding='utf8')
    for line in process.stdout:
        log.write(line)
        if line.startswith('Error'):
            log.close()
            sys.exit(1)
    # Check that the command succeeded.
    process.wait()
    if process.returncode != 0:
        for f in [log, sys.stderr]:
            print(f'Command failed, aborting.\n{command_str}', file=f)
        log.close()
        sys.exit(1)
    log.close()

#
# Generate a PHASE command
# Using the PHASE command as specified by Thom Nelson, and further edited by
# Niraj Rayamajhi.
# PHASE -MR locus41_42.phase.in locus41_42.phase.out 1000 1000 100 -d1 -x5 -l2
# -MR use PHASE recombination model
# Use 1000 iterations
# Thinning interval of 1000
# Burn in of 100
# -d1 use a multiallelic model without stepwise mutation (not microsats)
# -x5 Run algorithm 5 times
# -l2 limit to max 2 loci
# -F0.01 Minumim input hap frequency
# -O0.01 Minumim output hap frequency
def phase_command(phase_exe_path, phase_input_file, phase_output_file, recomb_model='MR', n_iters=1000, n_thinning=1000, n_burning=100, mutation_model='d1', n_runs=5, n_loci=2, hap_freq=0.01):
    # Check inputs
    assert os.path.exists(phase_exe_path)
    assert os.path.exists(phase_input_file)
    # Construst phase command
    phase_cmd = [phase_exe_path,
                 f'-{recomb_model}',
                 phase_input_file,
                 phase_output_file,
                 f'{n_iters}',
                 f'{n_thinning}',
                 f'{n_burning}',
                 f'-{mutation_model}',
                 f'-x{n_runs}',
                 f'-l{n_loci}',
                 f'-F{hap_freq}',
                 f'-O{hap_freq}']
    # Run command
    logf = f'{phase_output_file}.log'
    run_command(phase_cmd, logf)

#
# Determine the most likely haplotype
# TODO: What should min_prob be? PHASE docs are unclear on exact meaning of value.
def determine_best_phase_haplotype(phased_haps_list, min_prob=0.5):
    # Check input
    assert type(phased_haps_list) is list
    assert len(phased_haps_list) > 0
    assert 1.0 > min_prob > 0.0
    # Temp variables
    prev_prob = 0.0
    best_phase = None
    # Each element on the list correspond to the parsed line of the out_paired file from PHASE
    for hap in phased_haps_list:
        assert isinstance(hap, PhasedHaplotype)
        if hap.phase_prob >= prev_prob:
            prev_prob = hap.phase_prob
            best_phase = hap
    # After the loop, check that output is present
    assert best_phase is not None
    assert isinstance(best_phase, PhasedHaplotype)
    return best_phase

#
# Parse the PHASE out_pairs file
# Structure
#    IND: #1
#    0 1  , 0 3  , 1.000
#    IND: #2
#    0 1  , 0 2  , 1.000
#    IND: #3
#    1 0  , 2 0  , 1.000
#    IND: #4
#    1 0  , 2 0  , 1.000
def parse_phase_pairs(phase_pairs_file, individuals, for_haplotype, rev_haplotype, min_prob=0.5):
    # Check inputs
    assert os.path.exists(phase_pairs_file)
    assert type(individuals) is list
    assert len(individuals) > 0
    assert isinstance(rev_haplotype, RadHaplotype)
    assert isinstance(for_haplotype, RadHaplotype)
    # Output
    best_phase_list = list()
    # Temp variables
    indv_idx = None
    curr_hap = None
    # Open input file
    for line in open(phase_pairs_file, 'r'):
        line = line.strip('\n')
        if line[0] == '#':
            continue
        # Extract the individual index
        if line[0] == 'I':
            if curr_hap is None:
                indv_idx = int(line[6:])-1
                curr_hap = list()
            else:
                if len(curr_hap) == 0:
                    best_phase = PhasedHaplotype(for_haplotype.locid, rev_haplotype.locid, individuals[indv_idx], -1, -1, -1, -1, 0.0)
                    best_phase_list.append(best_phase)
                    indv_idx = int(line[6:])-1
                    curr_hap = list()
                else:
                    # Select the most likely phase
                    best_phase = determine_best_phase_haplotype(curr_hap, min_prob)
                    best_phase_list.append(best_phase)
                    indv_idx = int(line[6:])-1
                    curr_hap = list()
        else:
            # Structure of the line
            #    1 0  , 2 0  , 1.000
            fields = re.findall(r'(-?\d+)\s+(-?\d+)\s+,\s+(-?\d+)\s+(-?\d+)\s+,\s+(\d+\.\d+)', line)
            fields = fields[0]
            fa1 = int(fields[0]) # Forward locus allele 1
            ra1 = int(fields[1]) # Reverse locus allele 1
            fa2 = int(fields[2]) # Forward locus allele 2
            ra2 = int(fields[3]) # Reverse locus allele 2
            freq = float(fields[4])
            hap = PhasedHaplotype(for_haplotype.locid, rev_haplotype.locid, individuals[indv_idx], fa1, ra1, fa2, ra2, freq)
            curr_hap.append(hap)
    # For the last sample
    # Select the most likely phase
    if len(curr_hap) == 0:
        best_phase = PhasedHaplotype(for_haplotype.locid, rev_haplotype.locid, individuals[indv_idx], -1, -1, -1, -1, 0.0)
        best_phase_list.append(best_phase)
    else:
        best_phase = determine_best_phase_haplotype(curr_hap, min_prob)
        best_phase_list.append(best_phase)
    # Check output and return
    assert len(best_phase_list) == len(individuals)
    return best_phase_list

#
# Generate new, merged haplotype sequence
def reconfigure_haplotype_seq(for_allele, rev_allele, for_haplotype, rev_haplotype):
    # Check inputs
    assert isinstance(rev_haplotype, RadHaplotype)
    assert isinstance(for_haplotype, RadHaplotype)
    assert for_allele in range(len(for_haplotype.haps))
    assert rev_allele in range(len(rev_haplotype.haps))
    # Haplotype sequences
    for_seq = for_haplotype.haps[for_allele]
    # NOTE: In Stacks 2.56 the Haplotype VCF outputs has all sequences written in the orientation of the reference genome. The sequences for negative strand loci DO NOT need to be reverse complimented.
    rev_seq = rev_haplotype.haps[rev_allele]
    # Return the merged sequence
    return rev_seq+for_seq

#
# Merge the SNP columns of the two phased haplotypes
def transpose_phased_columns(for_locus, rev_locus, for_haplotype, rev_haplotype):
    # Check inputs
    assert isinstance(for_locus, StacksRadLocus)
    assert isinstance(rev_locus, StacksRadLocus)
    assert isinstance(rev_haplotype, RadHaplotype)
    assert isinstance(for_haplotype, RadHaplotype)
    # Output
    transposed_columns_list = list()
    # Process the reverse locus
    # NOTE: In Stacks 2.56 the Haplotype VCF outputs has the sequences and columns written in opposite orders for negative strand loci.
    ## The sequences are writen 5'->3' according to the reference (reverse complimented from how they appear in the catalog), but the columns are still based on the locus, thus 3'->5' when compared to the reference.
    ## Solution here is that the column list is accessed in the reverse order to match the position of the SNPs in the sequence.
    for column in rev_haplotype.cols:
        new_column = (len(rev_locus.seq)-1) - column
        transposed_columns_list.append(new_column)
    # Process the forward locus
    for column in for_haplotype.cols:
        new_column = (for_locus.min_bp - rev_locus.min_bp) + column
        transposed_columns_list.append(new_column)
    return transposed_columns_list

#
# Process the haplotypes for a single tag
def proc_single_tag_haps(locus, haplotype):
    assert isinstance(locus, StacksRadLocus)
    assert isinstance(haplotype, RadHaplotype)
    # Components of the merged haplotype object
    for_locus_id = None
    rev_locus_id = None
    chromosome   = locus.chrom
    basepair     = None
    genotypes    = list()
    haplotypes   = list()
    columns      = list()
    # Process forward and reverse loci differently
    if locus.dir == '+':
        for_locus_id = locus.locid
        basepair     = locus.min_bp
        genotypes    = haplotype.genos
        haplotypes   = haplotype.haps
        columns      = haplotype.cols
    else:
        rev_locus_id = locus.locid
        basepair     = locus.min_bp
        genotypes    = haplotype.genos
        haplotypes   = haplotype.haps # As of Stacks 2.56, this are ordered based on the reference
        for col in haplotype.cols:
            new_col = (len(locus.seq) - 1) - col
            columns.append(new_col)
    merged_rad_haplotype = MergedRadHaplotype(for_locus_id, rev_locus_id, chromosome, basepair, haplotypes, columns, genotypes)
    # Return the merged haplotype
    return merged_rad_haplotype

#
# Merge a pair of tags that only contain haplotypes on one end
# Transpose the SNP columns to the new locus lengths
def merge_single_hap_pair(for_locus, rev_locus, for_haplotype, rev_haplotype, individuals):
    assert isinstance(for_locus, StacksRadLocus)
    assert isinstance(rev_locus, StacksRadLocus)
    # Generate the components of the merged haplotype object
    for_locus_id = for_locus.locid
    rev_locus_id = rev_locus.locid
    chromosome   = for_locus.chrom
    basepair     = rev_locus.min_bp # This is the 5'-most position in the locus
    genotypes    = list()
    haplotypes   = list()
    columns      = list()
    # Ignore cases when both tags have no haps
    if rev_haplotype is None and for_haplotype is None:
        return None
    # If the reverse tag has no SNP
    elif rev_haplotype is None and for_haplotype is not None:
        assert isinstance(for_haplotype, RadHaplotype)
        haplotypes = for_haplotype.haps
        genotypes  = for_haplotype.genos
        for col in for_haplotype.cols:
            new_col = (for_locus.min_bp - rev_locus.min_bp) + col
            columns.append(new_col)
    # If the forward tag has no SNP
    elif rev_haplotype is not None and for_haplotype is None:
        assert isinstance(rev_haplotype, RadHaplotype)
        haplotypes = rev_haplotype.haps
        genotypes  = rev_haplotype.genos
        for col in rev_haplotype.cols:
            new_col = (len(rev_locus.seq)-1) - col
            columns.append(new_col)
    merged_rad_haplotype = MergedRadHaplotype(for_locus_id, rev_locus_id, chromosome, basepair, haplotypes, columns, genotypes)
    # Return the merged haplotype
    return merged_rad_haplotype


#
# Merge a pair of phased haplotypes into a single haplotype object
# Create all new possible sequence combinations
# Reconfigure the individual genotypes
# Transpose the SNP columns to the new locus lengths
def merge_phased_haplotypes(best_phase_list, for_locus, rev_locus, for_haplotype, rev_haplotype, individuals):
    # Check inputs
    assert len(best_phase_list) == len(individuals)
    assert isinstance(for_locus, StacksRadLocus)
    assert isinstance(rev_locus, StacksRadLocus)
    assert isinstance(rev_haplotype, RadHaplotype)
    assert isinstance(for_haplotype, RadHaplotype)
    # Output
    merged_rad_haplotype = None
    # Temp variables
    genos_dict = dict()
    haplotype_list = list()
    genotype_list  = list()
    g_idx = 0 # Index of the genotypes as they are seen
    # Loop over the individual phases, reconfigure each
    for i in range(len(best_phase_list)):
        indv = individuals[i]
        phase = best_phase_list[i]
        assert isinstance(phase, PhasedHaplotype)
        assert phase.sample_id == indv
        # Skip if a phase was not found
        # TODO: Still add the genotype even if Phase is poor (???)
        if phase.phase_found is False:
            for a in range(for_haplotype.ploid):
                genotype_list.append(None)
            continue
        # Loop over the alleles in the phase
        for allele in phase.allele_list:
            fal = allele[0] # For foward locus
            ral = allele[1] # For reverse locus
            al_str = f'{fal}-{ral}'
            if al_str not in genos_dict:
                # Add the allele if seen for the first time
                genos_dict[al_str] = g_idx
                # Generate a new sequence for that new allele
                new_hap_seq = reconfigure_haplotype_seq(fal, ral, for_haplotype, rev_haplotype)
                haplotype_list.append(new_hap_seq)
                g_idx += 1
            # Append the genotypes
            genotype_list.append(genos_dict[al_str])
    # Generate the components of the merged haplotype object
    for_locus_id = for_locus.locid
    rev_locus_id = rev_locus.locid
    chromosome   = for_locus.chrom
    basepair     = rev_locus.min_bp # This is the 5'-most position in the locus
    col_list = transpose_phased_columns(for_locus, rev_locus, for_haplotype, rev_haplotype)
    # MergedRadHaplotype object
    # TODO: Include information on phasing probability on the individual genotypes (0/1:0.95), as well as average phasing probability for the whole haplotype.
    merged_rad_haplotype = MergedRadHaplotype(for_locus_id, rev_locus_id, chromosome, basepair, haplotype_list, col_list, genotype_list)
    # Return the merged haplotype
    return merged_rad_haplotype
#
# Tally the proportion of samples that could be phased
def tally_phase_list(best_phase_list, min_phased_samples):
    # Check input
    assert type(best_phase_list) is list
    assert type(min_phased_samples) is float
    assert 1.0 >= min_phased_samples > 0.0
    phase_tally = 0
    for phase in best_phase_list:
        if phase.phase_found is True:
            phase_tally += 1
    phase_perc = phase_tally/len(best_phase_list)
    if phase_perc >= min_phased_samples:
        return True
    else:
        return False

#
# Remove the phase output files
def rm_phase_outs(phase_out_basename):
    cmd = f'rm -rf {phase_out_basename}.*'
    os.system(cmd)

#
# Process data and run PHASE
def phase_haplotypes(run_name, for_haplotype, rev_haplotype, individuals, phase_exe_path, outdir='.', min_phase_prob=0.5, min_phased_samples=0.8, dry_phase_run=False, delete_outs=False):
    # Check inputs
    assert isinstance(rev_haplotype, RadHaplotype)
    assert isinstance(for_haplotype, RadHaplotype)
    assert os.path.exists(outdir)
    outdir = outdir.rstrip('/')
    # assert os.path.exists(phase_exe_path)
    # Skip running Phase if performing a dry run
    if dry_phase_run is False:
        # 1. Convert the loci to a PHASE input file
        phase_str  = genotypes_to_phase(for_haplotype, rev_haplotype)
        phase_in_f = f'{outdir}/{run_name}.inp'
        out = open(phase_in_f, 'w')
        out.write(phase_str)
        out.close()
        # 2. Run PHASE
        phase_out_f = f'{outdir}/{run_name}.out'
        phase_command(phase_exe_path, phase_in_f, phase_out_f, 'MR', 1000, 1000, 100, 'd1', 5, 2, 0.05)
    # 3. Parse PHASE out_pairs file
    phase_pairs_f = f'{outdir}/{run_name}.out_pairs'
    best_phase_list = parse_phase_pairs(phase_pairs_f, individuals, for_haplotype, rev_haplotype, min_phase_prob)
    # 4. Check the PHASE output for missing data
    phased_haplotype = tally_phase_list(best_phase_list, min_phased_samples)
    # 5. Remove the output files (optional)
    if delete_outs:
        rm_phase_outs(f'{outdir}/{run_name}')
    # If the haplotype can be phased, return a list of the best phased
    if phased_haplotype is True:
        return best_phase_list
    else:
        return None

#
# Process overlapping genotypes
def process_overlapping_genotypes(loci_pairs_list, loci_dict, locus_haplotypes_dict, individuals, phase_exe_path, min_samples, min_prob, outdir='.', dry_phase_run=False, log=False, keep_single_tags=False, delete_outs=False):
    assert os.path.exists(outdir)
    outdir = outdir.rstrip('/')
    assert os.path.exists(phase_exe_path)
    # Output
    merged_loci_list = list()
    merged_haps_list = list()
    # Tallies for the log file.
    unpaired_loci   = 0
    single_tag_snps = 0
    phase_not_found = 0
    fully_phased    = 0
    phase_calls     = 0
    if log is True:
        print('Processing RAD loci pairs...', flush=True)
    # Loop over the loci pairs and find the overlapping ones:
    for i in range(len(loci_pairs_list)):
        pair = loci_pairs_list[i]
        assert isinstance(pair, LociPair)
        #
        # 1. Loci that that are not paired at all
        # If kept, the coordindates of the SNPs need to be adjusted
        #
        if pair.paired is False or None in [pair.for_loc, pair.rev_loc]:
            # Process the single tags when selected
            if keep_single_tags:
                # The loci "pair" can be added and processed later
                merged_loci_list.append(pair)
                # When the reverse is missing process the forward locus
                if pair.rev_loc is None:
                    # Get and process the haplotypes
                    for_loc = loci_dict[pair.for_loc]
                    for_hap = locus_haplotypes_dict.get(for_loc.locid, None)
                    if for_hap is not None:
                        merged_rad_haplotype = proc_single_tag_haps(for_loc, for_hap)
                        merged_haps_list.append(merged_rad_haplotype)
                elif pair.for_loc is None:
                    # Get and process the haplotypes
                    rev_loc = loci_dict[pair.rev_loc]
                    rev_hap = locus_haplotypes_dict.get(rev_loc.locid, None)
                    if rev_hap is not None:
                        merged_rad_haplotype = proc_single_tag_haps(rev_loc, rev_hap)
                        merged_haps_list.append(merged_rad_haplotype)
            unpaired_loci += 1
            continue
        # Loci that are paired should have both loci and haplotypes
        # Get Loci
        rev_loc = loci_dict[pair.rev_loc]
        assert isinstance(rev_loc, StacksRadLocus)
        for_loc = loci_dict[pair.for_loc]
        assert isinstance(for_loc, StacksRadLocus)
        # Get the haplotypes
        rev_hap = locus_haplotypes_dict.get(rev_loc.locid, None)
        for_hap = locus_haplotypes_dict.get(for_loc.locid, None)
        #
        # 2. Paired that are paired by position, but only have SNPs in one of the tags
        # If kept, they don't have to be phased, but the coordindates need to be
        # adjusted to the new longer locus
        #
        if rev_hap is None or for_hap is None:
            merged_loci_list.append(pair)
            merged_rad_haplotype = merge_single_hap_pair(for_loc, rev_loc, for_hap, rev_hap, individuals)
            if merged_rad_haplotype is not None:
                merged_haps_list.append(merged_rad_haplotype)
            single_tag_snps += 1
            continue
        # Prepare data and run PHASE when data is present
        name = f'MLocus_{for_loc.locid}_{rev_loc.locid}'
        phase_calls += 1
        best_phase_list = phase_haplotypes(name, for_hap, rev_hap, individuals, phase_exe_path, outdir, min_prob, min_samples, dry_phase_run, delete_outs)
        if log is True:
            if phase_calls%1000 == 0:
                print(f'    Running PHASE for the {phase_calls:,}-th haplotype pair.', flush=True)
        #
        # 3. These loci have data, but PHASE could not be found.
        # Treat as separate loci (???)
        #
        if best_phase_list is None:
            pair = LociPair(for_loc.locid, rev_loc.locid, False)
            # TODO: Add and process these loci (split again and adjust coordinates)
            # merged_loci_list.append(pair)
            phase_not_found += 1
            continue
        #
        # 4. Set of fully paired loci
        # Adjust SNP coordinates and genotypes into a single haplotype
        #
        else:
            # Create a single merge haplotype as a MergedRadHaplotype object
            merged_rad_haplotype = merge_phased_haplotypes(best_phase_list, for_loc, rev_loc, for_hap, rev_hap, individuals)
            merged_loci_list.append(pair)
            merged_haps_list.append(merged_rad_haplotype)
            fully_phased += 1
    if log is True:
        print('\nSummary of the Haplotype phasing...')
        print(f'    {unpaired_loci:,} loci had no adjacent tags.', flush=True)
        print(f'    {single_tag_snps:,} pairs of loci only had SNPs in a single tag.', flush=True)
        print(f'    {phase_not_found:,} loci pairs could not be phased.', flush=True)
        print(f'    {fully_phased:,} loci pairs could be phased.\n', flush=True)
    # Check output and return
    #assert len(loci_pairs_list) == len(merged_loci_list)
    return merged_loci_list, merged_haps_list

#
# Generate a new VCF containing the merged Haplotypes
def reformat_haplotype_vcf(merged_haps_list, individuals, basename, outdir, log=False):
    assert os.path.exists(outdir)
    outdir = outdir.rstrip('/')
    nlocs=0
    # Generate the new output VCF
    outvcf = open(f'{outdir}/{basename}.merged_haps.vcf', 'w')
    # Format and write VCF Header
    indv_str = '\t'.join(individuals)
    header_str = f'##fileformat=VCFv4.2\n##fileDate={DATE}\n##source={PROG}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{indv_str}\n'
    outvcf.write(header_str)
    # Loop over all haplotypes
    for haplotype in merged_haps_list:
        assert isinstance(haplotype, MergedRadHaplotype)
        vcf_row = haplotype.write_vcf_line()
        outvcf.write(f'{vcf_row}\n')
        nlocs+=1
    outvcf.close()
    if log is True:
        print(f'Generating merged loci VCF.\n    Wrote haplotypes for {nlocs:,} loci.\n', flush=True)

#
# Report all options
#
def print_options(args):
    print(f'''Selected options:
    Run Name:                  {args.run_name}
    Restriction Enzyme:        {args.res_enzyme}
    Max Haplotype Length:      {args.max_sites_in_hap}
    Max Number of Alleles:     {args.max_alleles_in_loc}
    Min Phasing Probability:   {args.min_phase_prob}
    Min Proportion of Samples: {args.min_number_samples}''')
    if args.delete_phase_outs:
        print('    Deleting PHASE output files.')
    if args.keep_single_tags:
        print('    Keeping Loci/Haplotypes for Single Tags.')
    if args.phase_dry_run:
        print('    Do Not Rerun PHASE (Dry Run).')
    print('', flush=True)

#
# Main Function
#
def main():
    # Parse arguments and generate files
    args        = parse_args()
    renz        = RestrictionEnzyme(args.res_enzyme)
    loci_fa     = args.stacks_loci_fasta
    vcf         = args.stacks_haps_vcf
    outdir      = args.out_dir
    runame      = args.run_name
    max_hap_len = args.max_sites_in_hap
    max_alleles = args.max_alleles_in_loc
    dry         = args.phase_dry_run
    phase_exe   = args.phase_exe_path
    min_phase_p = args.min_phase_prob
    min_n_sams  = args.min_number_samples
    keep_single = args.keep_single_tags
    delete_outs = args.delete_phase_outs
    log_file    = f'{outdir}/phase_rad_loci.log'

    sys_stdout_bak = sys.stdout
    try:
        sys.stdout = open(log_file, 'w')
        log = True

        # Start program
        print(f'phase_rad_loci started on {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n', flush=True)

        # Print command line arguments
        print('Command used:', flush=True)
        arg_str = ' '.join(sys.argv)
        print(f'{arg_str}\n', flush=True)
        # Print the general options
        print_options(args)

        # Parse the Stacks Locus Fasta
        loci_dict = parse_loci_fasta(loci_fa, renz, log)
        # Find loci that come from the same restriction site
        loci_pairs_list = find_overlapping_loci(loci_dict, renz, log)
        # Load variants
        locus_haplotypes_dict, individuals = parse_stacks_haps_vcf(vcf, log)
        # Filter variants and create some haplotype stats
        locus_haplotypes_dict = filter_haps_statistics(locus_haplotypes_dict, max_hap_len, max_alleles, outdir, log)
        # Process the loci that are overlapping
        merged_loci_list, merged_haps_list = process_overlapping_genotypes(loci_pairs_list, loci_dict, locus_haplotypes_dict, individuals, phase_exe, min_n_sams, min_phase_p, outdir, dry, log, keep_single, delete_outs)
        # Generate the merged loci fasta
        reformat_loci_fasta(loci_dict, merged_loci_list, renz, runame, outdir, log)
        # Generate the merged haplotype VCF
        reformat_haplotype_vcf(merged_haps_list, individuals, runame, outdir, log)

        # Finish program
        print(f'phase_rad_loci completed on {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
    finally:
        sys.stdout = sys_stdout_bak

# Run Code
if __name__ == '__main__':
    main()
