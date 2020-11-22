#!/usr/bin/env python

import sys
import argparse
import requests

ANNOTATIONS = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch"
ACCEPT = {"Accept": "text/tsv"}

def get_missense_genes(filename):
    """
    Returns a set of all genes with missense variants from snpEff summary
    geneId file.
    """
    missense_genes = set()
    with open(filename) as vars:
        for line in vars:
            if line[0] == "#":
                continue
            else:
                var_info = line.split("\t")
                if int(var_info[13]) > 0:
                    missense_genes.add(var_info[1])
    return missense_genes


def read_gene_variants(filename):
     """
     Returns a set of all genes with variants from snpEff summary
     geneId file.
     """
     genes = set()
     with open(filename) as vars:
        for line in vars:
            if line[0] == "#":
                continue
            else:
                var_info = line.split("\t")
                genes.add(var_info[1])
     return genes


def get_all_go_genes(go_id):
    """Returns a set of all genes annotated with go_id."""
    params = {"selectedFields": "symbol", "goId": go_id}
    r = requests.get(ANNOTATIONS, headers=ACCEPT, params=params)
    if not r.ok:
        r.raise_for_status()
        sys.exit(1)
    associated_genes = [i.upper() for i in r.text.strip("SYMBOL\n").split("\n")]
    return associated_genes


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


def find_matching_genes(query_genes, go_genes):
    """
    Returns a list of all genes in query_genes that are present in go_genes."""
    matching_genes = []
    for gene in query_genes:
        if gene in go_genes:
            matching_genes.append(gene)
    return matching_genes


def get_args():
    """Get command line arguments."""
    parser = argparse.ArgumentParser(description="Parse variants in genes with a given Gene Ontology ID")
    parser.add_argument("vars",
                        type=str,
                        help="snpEff annotated variants tsv file with geneId")
    parser.add_argument("qs",
                        type=str,
                        help="Tsv file containing PAO1 genes with the 'quorum sensing' GO accession")
    return parser.parse_args()


def main():
    args = get_args()
    # missense_genes = get_missense_genes(args.vars)
    pao_qs_genes = read_qs_genes(args.qs)
    genes_with_variants = read_gene_variants(args.vars)
    for gene in genes_with_variants:
        if gene in pao_qs_genes:
            print(gene)
    # print(genes_with_variants)
    # go_genes = get_all_go_genes(args.go)
    # matching_genes = find_matching_genes(genes_with_variants, go_genes)
    # print("go_id\tnum_genes")
    # print("{0}\t{1}".format(args.go, len(matching_genes)))


if __name__ == "__main__":
    main()
