#!/usr/bin/env python3
# 2018 Nicolas Rochette
import sys, os, re, gzip

usage='''\
Usage:
  {} gstacks.details.gz

Creates a histogram of insert sizes, according to the aln_matrix section. Outputs
to `gstacks.details.ins_z.tsv`.
'''

cig_ops = '=DHIMNPSX'
op_re = re.compile('([0-9]+)({})'.format('|'.join(cig_ops)))
cigar_re = re.compile('^([0-9]+({}))+$'.format('|'.join(cig_ops)))
def insert_size(cig):
    assert cigar_re.match(cig)
    cig = op_re.findall(cig)
    z = 0
    for n, op in cig:
        if op in 'MD=XN':
            z += int(n)
    return z

def main():
    if len(sys.argv) != 2:
        print(usage, end='', file=sys.stderr)
        sys.exit(1)
    gstacks_details = sys.argv[1]
    hist = {}
    details = gzip.open(gstacks_details, 'rt')
    # We process the lines of the the aln_matrix section, that are preceeded by
    # a 'overlap [>0]' line, and where the read name includes '/m'.
    ref_based = bool(re.search(r'\.ref/', os.getcwd()))
    while True:
        line = details.readline()
        if not line:
            break
        if ref_based:
            if not line.startswith('BEGIN aln_matrix'):
                continue
        else:
            if not (line.startswith('overlap\t')
                    or line.startswith('one_debrjuin_subgraph\t')):
                continue
            if line.startswith('overlap\t0\t'):
                continue
            while True:
                line = details.readline()
                assert line
                if line.startswith('BEGIN aln_matrix'):
                    break
        # Read the contents of the matrix.
        while True:
            line = details.readline()
            assert line
            if line.startswith('END aln_matrix'):
                break
            fields = line.strip().split()
            assert len(fields) == 3
            if fields[0][-2:] != '/m':
                continue
            cig = fields[2]
            if cig[-1] == 'D':
                # Discard that operation.
                i = 2
                while cig[-i] in '0123456789':
                    i += 1
                cig = cig[:-i+1]
            z = insert_size(cig)
            hist.setdefault(z, 0)
            hist[z] += 1
    assert hist
    tsv = open(gstacks_details[:-3] + '.ins_z.tsv', 'w')
    tsv.write('ins_z\tn_pairs\n')
    for z in range(1, max(hist.keys())):
        hist.setdefault(z, 0)
        tsv.write('{}\t{}\n'.format(z, hist[z]))

if __name__ == '__main__':
    main()
