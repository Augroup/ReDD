"""Compare two ReDD molecule-level prediction files read-by-read.

Both files have the ReDD layout: <label> <read_id> <contig> <pos> <strand> <prob>.
Records are matched on (read_id, contig, pos); --offset / --rename map the coordinates of
file A onto file B (e.g. a sub-reference "chr11_sub" with positions shifted by -700000
is compared with hg38 "chr11" via --rename chr11_sub=chr11 --offset 700000).

Usage: python compare_predictions.py A.txt B.txt [--offset N] [--rename A_contig=B_contig ...] [--tol 1e-3]
"""
import argparse
import gzip
import sys


def load(path, offset=0, rename=None, contigs=None):
    rename = rename or {}
    opener = gzip.open if path.endswith(".gz") else open
    recs = {}
    with opener(path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            read_id, contig, pos, prob = fields[1], fields[2], int(fields[3]), float(fields[5])
            contig = rename.get(contig, contig)
            if contigs is not None and contig not in contigs:
                continue
            recs[(read_id, contig, pos + offset)] = prob
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--offset", type=int, default=0, help="added to positions of A")
    ap.add_argument("--rename", nargs="*", default=[], help="A_contig=B_contig")
    ap.add_argument("--tol", type=float, default=1e-3, help="max |prob_A - prob_B| to count as identical")
    ap.add_argument("--reads", default=None, help="optional file of read ids to restrict B to")
    args = ap.parse_args()

    rename = dict(x.split("=", 1) for x in args.rename)
    a = load(args.a, args.offset, rename)
    contigs = {k[1] for k in a}
    b = load(args.b, contigs=contigs)
    if args.reads:
        keep = {l.split()[0] for l in open(args.reads) if l.strip()}
        b = {k: v for k, v in b.items() if k[0] in keep}

    common = sorted(set(a) & set(b))
    only_a, only_b = len(set(a) - set(b)), len(set(b) - set(a))
    diffs = [abs(a[k] - b[k]) for k in common]
    n_ok = sum(d <= args.tol for d in diffs)
    print(f"A records: {len(a)}   B records (restricted): {len(b)}   matched: {len(common)}   only A: {only_a}   only B: {only_b}")
    if common:
        diffs_sorted = sorted(diffs)
        print(f"|dprob| max: {diffs_sorted[-1]:.6f}   median: {diffs_sorted[len(diffs)//2]:.6f}   "
              f"within tol {args.tol}: {n_ok}/{len(common)} ({n_ok/len(common):.2%})")
        worst = sorted(common, key=lambda k: -abs(a[k] - b[k]))[:5]
        for k in worst:
            print(f"  worst: {k}  A={a[k]:.6f}  B={b[k]:.6f}")
    return 0 if common and n_ok == len(common) else 1


if __name__ == "__main__":
    sys.exit(main())
