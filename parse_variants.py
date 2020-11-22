#!/usr/bin/env python

import re
import sys
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def read_gene_variants(filename):
    """
    Returns a dictionary of all variant types and corresponding genes from snpEff summary
    geneId file.
    """
    variant_map = {}  # column_index -> variant_type
    variants = {}  #  type -> [locus_tag]
    var_re = "variants_effect_([a-z_]+)"
    with open(filename) as vars:
        vars.readline()  # Skip first line about formatting
        header = vars.readline().strip().split("\t")
        for i, col in enumerate(header[8:]):  # Get mapping info for variants
            var_type = re.search(var_re, col).group(1)
            variant_map[i] = var_type
            variants[var_type] = []
        for line in vars:
            line_content = line.split("\t")
            var_info = line_content[8:]
            gene_id = line_content[1]
            for i, var in enumerate(var_info):
                if int(var) > 0:
                    variants[variant_map[i]].append(gene_id)
    return variants


def read_qs_genes(go_gene_file):
    """
    Returns a dictionary of Pseudomonas genes with a certain GO accession from a tsv file
    and associated information.
    """
    gene_info = {}  # locus_tag -> (gene_name, product_description)
    with open(go_gene_file) as genes:
        genes.readline()
        for line in genes:
            line_content = line.strip().split("\t")
            gene_info[line_content[0]] = (line_content[1], line_content[2])
    return gene_info


def print_qs_variants(variants, qs_genes, outfile):
    """
    Prints variants in query_genes that are in qs_genes in tsv format to outfile.
    Format: locus_id    gene_name   mutation_type   product_description
    """
    if outfile is None:
        o = sys.stdout
    else:
        o = open(outfile, "w+")
    print("locus_id\tgene_name\tmutation_type\tproduct_description", file=o)
    for var_type in variants:
        for locus in variants[var_type]:
            if locus in qs_genes:
                print(locus, qs_genes[locus][0], var_type, qs_genes[locus][1], sep="\t", file=o)


def plot_variants(variants):



def get_args():
    """Get command line arguments."""
    parser = argparse.ArgumentParser(description="Parse variants in genes with a given Gene Ontology ID")
    parser.add_argument("vars",
                        type=str,
                        help="snpEff annotated variants tsv file with geneId")
    parser.add_argument("qs",
                        type=str,
                        help="Tsv file containing PAO1 genes with the 'quorum sensing' GO accession")
    parser.add_argument("-o", "--outfile",  # Change to prefix for all files
                        type=str,
                        default=None,
                        help="Output file for quorum sensing variants to be printed to [stdout]")
    return parser.parse_args()


def main():
    args = get_args()
    pao_qs_genes = read_qs_genes(args.qs)
    genes_with_variants = read_gene_variants(args.vars)
    print_qs_variants(genes_with_variants, pao_qs_genes, args.outfile)


if __name__ == "__main__":
    main()

# Add plotting
# Print all variant types to file
# Make stacked bar chart