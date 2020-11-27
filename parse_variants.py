#!/usr/bin/env python

import re
import sys
import argparse
import pandas as pd
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


def count_variant_types(variants):
    """Returns a dataframe with columns [type of variant, count]."""
    variant_types = []
    variant_counts = []
    for v in variants:
        variant_types.append(v)
        variant_counts.append(len(variants[v]))
    var_df = pd.DataFrame(zip(variant_types, variant_counts), columns=["type", "count"]).sort_values(by="count", ascending=False)
    return var_df


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


def get_qs_variants(variants, qs_genes, prefix):
    """
    Prints variants in query_genes that are in qs_genes in tsv format to outfile.
    Also returns a list of mutation types corresponding to the qs gene variants.
    Format: locus_id    gene_name   mutation_type   product_description
    """
    outfile = prefix + ".qs_variants.tsv"
    mutation_type = []
    # Check if file is writable
    with open(outfile, "w+") as o:
        print("locus_id\tgene_name\tmutation_type\tproduct_description", file=o)
        for var_type in variants:
            for locus in variants[var_type]:
                if locus in qs_genes:
                    mutation_type.append(var_type)
                    print(locus, qs_genes[locus][0], var_type, qs_genes[locus][1], sep="\t", file=o)
    return mutation_type


def plot_var_counts(variant_counts, prefix):
    """Prints a bar plot of variant counts to a file."""
    outfile = prefix + ".variant_types.png"
    # Check if file is writable
    bar = sns.barplot(y=variant_counts["type"], x=variant_counts["count"],
        color="grey", label="big")
    bar.set(xscale="log")
    bar.set_xlabel(xlabel="Count", fontsize=15)
    bar.set_ylabel(ylabel="Variant Type", fontsize=15)
    bar.set_title("Variant Types and Frequencies", fontsize=20)
    bar.tick_params(labelsize=15)
    bar_plot = bar.get_figure()
    bar_plot.savefig(outfile, bbox_inches="tight")
    bar_plot.close()


def plot_qs_variants(qs_variants, prefix):
    """Prints a histogram of variant types in quorum sensing genes."""
    outfile = prefix + ".qs_variants.png"
    # Check if file is writable
    qs_var_df = pd.DataFrame(qs_variants, columns=["mutation_type"])
    qs = sns.histplot(qs_var_df, y="mutation_type", color="grey", alpha=1)
    qs.set_title("Types of Quorum Sensing Gene Variants", fontsize=20)
    qs.set_xlabel("Count", fontsize=15)
    qs.set_ylabel("Variant Type", fontsize=15)
    qs.tick_params(labelsize=15)
    qs_plot = qs.get_figure()
    qs_plot.savefig(outfile, bbox_inches="tight")
    qs_plot.close()


def get_args():
    """Get command line arguments."""
    parser = argparse.ArgumentParser(description="Parse variants in genes with a given Gene Ontology ID")
    parser.add_argument("vars",
                        type=str,
                        help="snpEff annotated variants tsv file with geneId")
    parser.add_argument("qs",
                        type=str,
                        help="Tsv file containing PAO1 genes with the 'quorum sensing' GO accession")
    parser.add_argument("-p", "--prefix",
                        type=str,
                        default="PAO1",
                        help="Prefix for output files to be printed to [PAO1]")
    return parser.parse_args()


def main():
    args = get_args()
    # Parse all variants
    genes_with_variants = read_gene_variants(args.vars)
    # Get variant counts
    variant_counts = count_variant_types(genes_with_variants)
    sns.set_theme(style="whitegrid")
    # Plot variant counts
    # plot_var_counts(variant_counts, args.prefix)
    # Find QS genes with variants
    qs_gene_variants = read_qs_genes(args.qs)
    # Print to tsv file
    qs_variant_counts = get_qs_variants(genes_with_variants, qs_gene_variants, args.prefix)
    # Plot QS variants
    plot_qs_variants(qs_variant_counts, args.prefix)


if __name__ == "__main__":
    main()
