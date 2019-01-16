#!/usr/bin/env python

import argparse
import re


def convert_gff_to_gtf(gff_path, gtf_path):
    t = re.compile(r'ID=[0-9]+_([0-9]+);')
    t2 = re.compile(r'ID=cds([0-9]+);')
    with open(gff_path, "r") as gff:
        with open(gtf_path, "w") as gtf:
            for line in gff:
                if line.startswith("#"):
                    continue
                toks = line.strip().split("\t")
                try:
                    orf = t.findall(toks[-1])[0]
                except:
                    try:
                        orf = t2.findall(toks[-1])[0]
                    except:
                        continue
                gene_id = toks[0] + "_" + orf
                toks[-1] = 'gene_id "%s"; %s' % (gene_id, toks[-1])
                gtf.write("\t".join(toks) + "\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts GFF to GTF from Prodigal output.')

    # required
    parser.add_argument('gff', help='Path to input file.')
    parser.add_argument('gtf', help='Path to output file.')

    args = parser.parse_args()

    convert_gff_to_gtf(args.gff, args.gtf)

