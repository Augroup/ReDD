"""Convert aligned reads back into an unaligned Dorado-style BAM keeping the basecaller tags.

Used by make_test_data.sh. Sequences are restored to basecall orientation; Dorado tags
(mv, ts, pi, sp, ns, ...) are kept, minimap2 tags (including its conflicting ts:A) are dropped.
"""
import sys

import pysam

KEEP = {"mv", "pi", "sp", "ns", "MN", "qs", "du", "sm", "sd", "sv", "RG"}


def main(in_bam, out_bam):
    inp = pysam.AlignmentFile(in_bam, "rb")
    header = {"HD": {"VN": "1.6", "SO": "unsorted"},
              "PG": [{"ID": "ReDD_test_data", "PN": "make_test_data",
                      "DS": "K562 RNA004 reads, chr11:700001-860000 subset, tags from Dorado basecalling"}]}
    out = pysam.AlignmentFile(out_bam, "wb", header=header)
    comp = str.maketrans("ACGTN", "TGCAN")
    n = 0
    for r in inp:
        a = pysam.AlignedSegment(out.header)
        a.query_name = r.query_name
        a.flag = 4
        seq = r.query_sequence
        qual = list(r.query_qualities) if r.query_qualities is not None else None
        if r.is_reverse:
            seq = seq.translate(comp)[::-1]
            if qual is not None:
                qual = qual[::-1]
        a.query_sequence = seq
        a.query_qualities = qual
        tags = []
        for tag, val, vtype in r.get_tags(with_value_type=True):
            if tag == "ts":
                if isinstance(val, int):  # Dorado ts:i (trimmed samples); minimap2's ts:A is dropped
                    tags.append(("ts", val, "i"))
            elif tag in KEEP:
                tags.append((tag, val) if vtype == "B" else (tag, val, vtype))
        a.set_tags(tags)
        out.write(a)
        n += 1
    out.close()
    inp.close()
    print(f"wrote {n} reads to {out_bam}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
