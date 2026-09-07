"""Build a read_id -> signal file index for `uncalled4 align -x`.

uncalled4 auto-generates this index inside the signal directory when it is missing, which
(a) requires write access to the raw-data directory and (b) races when several jobs start at
once. The pipeline builds it once here instead.

Usage: python build_read_index.py <pod5|slow5|blow5|fast5 file or directory> [more paths...] -o read_index.txt
"""
import argparse
import os
import sys

EXTS = (".pod5", ".slow5", ".blow5", ".fast5")


def iter_files(paths, recursive):
    for path in paths:
        if os.path.isdir(path):
            if recursive:
                for root, _, files in os.walk(path):
                    for f in sorted(files):
                        if f.endswith(EXTS):
                            yield os.path.join(root, f)
            else:
                for f in sorted(os.listdir(path)):
                    if f.endswith(EXTS):
                        yield os.path.join(path, f)
        elif path.endswith(EXTS):
            yield path
        else:
            raise ValueError(f"not a signal file or directory: {path}")


def read_ids(path):
    if path.endswith(".pod5"):
        import pod5
        with pod5.Reader(path) as reader:
            return [str(r) for r in reader.read_ids]
    if path.endswith((".slow5", ".blow5")):
        import pyslow5
        f = pyslow5.Open(path, "r")
        try:
            return list(f.get_read_ids()[0])
        finally:
            f.close()
    if path.endswith(".fast5"):
        from ont_fast5_api.fast5_interface import get_fast5_file
        with get_fast5_file(path, mode="r") as f:
            return list(f.get_read_ids())
    raise ValueError(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-r", "--recursive", action="store_true")
    args = parser.parse_args()

    n_files = n_reads = 0
    with open(args.output, "w") as out:
        out.write("read_id\tfilename\n")
        for path in iter_files(args.paths, args.recursive):
            ids = read_ids(path)
            for rid in ids:
                out.write(f"{rid}\t{path}\n")
            n_files += 1
            n_reads += len(ids)
            sys.stderr.write(f"{path}: {len(ids)} reads\n")
    if n_files == 0:
        raise SystemExit("no signal files found")
    sys.stderr.write(f"indexed {n_reads} reads in {n_files} files -> {args.output}\n")


if __name__ == "__main__":
    main()
