import argparse
import os
import sys
from subprocess import Popen, PIPE

import pandas as pd

from Bio import SeqIO


DEFAULT_THRES = 5
COLUMNS = ["target.name", "target.accession",
           "query.name", "query.accession",
           "full.E-value", "full.score", "full.bias",
           "best.E-value", "best.score", "best.bias",
           "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc",
           "target.description"]


def process_query(query_fn):
    result_seqs = SeqIO.to_dict(SeqIO.parse(query_fn, "fasta"))
    print(f'    There are {len(result_seqs.values())} protein(s) before filtering.', file=sys.stderr)
    for seqobj in result_seqs.values():
        seqobj.description = seqobj.description.replace(" ", "_")
    with open("###query.fasta", "w") as tempfile:
        SeqIO.write(result_seqs.values(), tempfile, "fasta")
    return pd.Series(result_seqs)


def run_hmm(hmm):
    proc = Popen(["hmmsearch", "--tblout", "###table.tab", hmm, "###query.fasta"], stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate()
    table = pd.read_csv("###table.tab", sep="\s+", header=None, names=COLUMNS, comment="#")
    if table.empty:
        sys.stderr.write("No viable hits found!\n")
        sys.exit(1)
    return table


def remove_temp():
    if os.path.isfile("###table.tab"):
        os.remove("###table.tab")
    if os.path.isfile("###query.fasta"):
        os.remove("###query.fasta")


def choose_at_threshold(seqs, table, thres=DEFAULT_THRES):
    seqs_at_thres = seqs[table["target.name"][table["best.E-value"] <= thres].to_list()].to_list()
    for seqobj in seqs_at_thres:
        seqobj.id = seqobj.id.partition("/")[0]
        seqobj.description = ""
    return seqs_at_thres


def score_hmm(query_fn, out_fn, hmm, thres=DEFAULT_THRES):
    result_seqs = process_query(query_fn)
    result_table = run_hmm(hmm)
    remove_temp()
    seqs_at_thres = choose_at_threshold(result_seqs, result_table, thres=thres)
    lng = SeqIO.write(seqs_at_thres, out_fn, "fasta")
    print(f"    {lng} protein(s) are filtered.", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filters sequences with an HMM",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--query",
                        metavar="fasta",
                        required=True,
                        help="input fasta file")
    parser.add_argument("--out",
                        metavar="fasta",
                        required=True,
                        help="output fasta file")
    parser.add_argument("--hmm",
                        metavar="hmm",
                        required=True,
                        help="HMM file")
    parser.add_argument("--thres",
                        default=DEFAULT_THRES,
                        metavar="x",
                        type=float,
                        help=f"E-value threshold")

    args = parser.parse_args()

    score_hmm(query_fn=args.query,
              out_fn=args.out,
              hmm=args.hmm,
              thres=args.thres)
