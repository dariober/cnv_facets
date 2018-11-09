#!/usr/bin/env python3

import argparse
import pysam
import os

VERSION= '0.1.0'

parser = argparse.ArgumentParser(description= """
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--vcf', '-V',
    default= '-',
    help='''Input vcf file. Default: %(default)s''') 

parser.add_argument('--tumour-bam', '-t',
    required= True,
    help='''Tumour bam file''')

parser.add_argument('--normal-bam', '-n',
    required= True,
    help='''Normal bam file''')

parser.add_argument('--min-read-count-t', '-rt',
    type= int,
    default= 0,
    help='''\
Output SNPs with at least this many reads in tumour sample. 
Default= %(default)s''')

parser.add_argument('--min-read-count-n', '-rn',
    type= int,
    default= 0,
    help='''\
Output SNPs with at least this many reads in normal sample.  Default= %(default)s''')

parser.add_argument('--pseudo-snps', '-P',
    type= int,
    default= 0,
    help='''\
Add a pseudo count every so many bases if no SNP is found. Disabled if 0.
Default= %(default)s''')

parser.add_argument('--min-base-quality', '-Q',
    type= int,
    default= 0,
    help='''\
Do not count bases with quality below this threshold. Default= %(default)s''')

parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + VERSION)

args= parser.parse_args()

samfile= {'tumour': pysam.AlignmentFile(args.tumour_bam, "rb"),
          'normal': pysam.AlignmentFile(args.normal_bam, "rb")}
vcf= pysam.VariantFile(args.vcf)

# --------------------------

def is_valid_variant(variant):
    """True if variant can be used for genotyping. Make this consistent with
    the original snp-pileup
    """
    if len(variant.ref) > 1 or len(variant.alts) > 1 or len(variant.alts[0]) > 1:
        return False
    return True

def get_pileup(variant, pseudo_snps, pileupcolumn, min_read_count):
    """
    """
    pileup_pos= pileupcolumn.pos+1

    if pileup_pos == variant.pos:
        pseudo= False
        ref= variant.ref
        alt= variant.alts[0]
    elif pseudo_snps > 0 and pileup_pos % pseudo_snps == 0:
        ref= '.'
        alt= '.'
        pseudo= True
    else:
        # Pileup position does not match the variant or the pseudo SNP
        return None

    nref= 0
    nalt= 0
    nother= 0
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            if pseudo or pileupread.alignment.query_sequence[pileupread.query_position] == ref:
                nref += 1
            elif pileupread.alignment.query_sequence[pileupread.query_position] == alt:
                nalt += 1
            else:
                nother += 1
    xpass= True
    if (nref + nalt) < min_read_count:
        xpass= False
    return {'chrom': pileupcolumn.reference_name, 'pos': pileup_pos, 'nref': nref, 'nalt': nalt, 'nother': nother, 'ref': ref, 'alt': alt, 'pass': xpass}

def printer(snp_1, snp_2):
    """Return a printable string from the two dicts of pileup data
    """
    if (snp_1 is None or not snp_1['pass']) or (snp_2 is None or not snp_2['pass']):
        return None

    if snp_2 is None:
        return ','.join([str(snp_1['chrom']), str(snp_1['pos']), snp_1['ref'], snp_1['alt'], 
            str(snp_1['nref']), str(snp_1['nalt']), str(snp_1['nother']), 
            '0', '0', '0'])
    
    if snp_1 is None:
        return ','.join([str(snp_2['chrom']), str(snp_2['pos']), snp_2['ref'], snp_2['alt'], 
            '0', '0', '0',
            str(snp_2['nref']), str(snp_2['nalt']), str(snp_2['nother'])])

    assert snp_1['chrom'] == snp_2['chrom']
    assert snp_1['pos'] == snp_2['pos']

    return ','.join([str(snp_1['chrom']), str(snp_1['pos']), snp_1['ref'], snp_1['alt'], 
        str(snp_1['nref']), str(snp_1['nalt']), str(snp_1['nother']), 
        str(snp_2['nref']), str(snp_2['nalt']), str(snp_2['nother'])])

# --------------------------

prev_chrom= None
prev_varpos= None
print(','.join(['Chromosome', 'Position', 'NormalR', 'NormalA', 'NormalE', 'TumourR', 'TumourA', 'TumourE']))
while True:
    try:
        variant= next(vcf)
    except StopIteration:
        break

    if not is_valid_variant(variant):
        continue

    if prev_chrom is None or prev_chrom != variant.chrom:
        prev_chrom= variant.chrom
        pileups= {'tumour': samfile['tumour'].pileup(variant.chrom, 1, min_base_quality= args.min_base_quality), 
                  'normal': samfile['normal'].pileup(variant.chrom, 1, min_base_quality= args.min_base_quality)}
        column_t= next(pileups['tumour'])
        column_n= next(pileups['normal'])

    if prev_varpos is not None and (variant.pos == prev_varpos and prev_chrom == variant.chrom):
        # If multiple variants occur at the same position, process the first and skip the duplicates
        continue
  
    while column_t.pos+1 <= variant.pos and column_n.pos+1 <= variant.pos:
        # Walk along the pileup files until you reach the variant position
        snp_t= None
        snp_n= None
        if column_t.pos < column_n.pos:
            snp_t= get_pileup(variant, args.pseudo_snps, column_t, args.min_read_count_t)
            column_t= next(pileups['tumour'])
        elif column_t.pos == column_n.pos:
            snp_t= get_pileup(variant, args.pseudo_snps, column_t, args.min_read_count_t)
            snp_n= get_pileup(variant, args.pseudo_snps, column_n, args.min_read_count_n)
            column_t= next(pileups['tumour'])
            column_n= next(pileups['normal'])
        elif column_t.pos > column_n.pos:
            snp_n= get_pileup(variant, args.pseudo_snps, column_n, args.min_read_count_n)
            column_n= next(pileups['normal'])
        else:
            raise Exception()
        line= printer(snp_n, snp_t)
        if line is not None:
            print(line)

    prev_varpos= variant.pos
