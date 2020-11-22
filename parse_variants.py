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


def get_all_go_genes(go_id):
    """Returns a set of all genes annotated with go_id."""
    params = {"selectedFields": "symbol", "goId": go_id}
    r = requests.get(ANNOTATIONS, headers=ACCEPT, params=params)
    if not r.ok:
        r.raise_for_status()
        sys.exit(1)
    associated_genes = set(r.text.strip("SYMBOL\n").split("\n"))
    return associated_genes


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
    parser.add_argument("go",
                        type=str,
                        help="Gene Ontology ID query")
    return parser.parse_args()


def main():
    args = get_args()
    missense_genes = get_missense_genes(args.vars)
    go_genes = get_all_go_genes(args.go)
    matching_genes = find_matching_genes(missense_genes, go_genes)
    print(matching_genes)


if __name__ == "__main__":
    main()
