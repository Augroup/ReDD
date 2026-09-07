"""Molecule-level A-to-I prediction for RNA004 ReDD feature files.

Input : one or more HDF5 files produced by `uncalled4 align --redd-out` (datasets X, y_ref, y_call, info).
Output: tab-separated text, one line per read x A-site:
        <label>\t<read_id>\t<contig>\t<pos_1based>\t<strand>\t<probability>
        (the first five columns are the `info` string written by uncalled4; the label is
        "I:<ratio>" for candidate sites when --redd-candidate was given, otherwise "I").

The checkpoint is a dict with keys `model_state_dict`, `model_config` (kwargs for ReDDModel),
and `output_head` (which logit column of the modification head to report).
"""
import argparse
import gzip
import os
import sys
import time

import h5py
import hdf5plugin  # noqa: F401  (registers the LZ4 filter used by uncalled4 output)
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import ReDDModel  # noqa: E402


def gpu_supported_by_torch(device_index=0):
    """True if this torch build has kernels for the GPU's compute capability.

    A cubin built for sm_XY runs on devices with the same major version and minor >= Y; a device
    with a newer major version (e.g. sm_120 Blackwell with a cu124 build for sm_50..sm_90) has no
    usable kernel image and fails with 'no kernel image is available for execution on the device'.
    """
    major, minor = torch.cuda.get_device_capability(device_index)
    for arch in torch.cuda.get_arch_list():
        if arch.startswith("sm_"):
            a_major, a_minor = int(arch[3:-1]), int(arch[-1])
            if a_major == major and a_minor <= minor:
                return True
    return False


def resolve_device(name, force_gpu=False):
    name = name.lower()
    if name in ("gpu", "cuda"):
        if not torch.cuda.is_available():
            sys.stderr.write("WARNING: GPU requested but CUDA is not available, falling back to CPU\n")
            return torch.device("cpu")
        if not force_gpu and not gpu_supported_by_torch():
            major, minor = torch.cuda.get_device_capability()
            sys.stderr.write(
                f"WARNING: {torch.get_device_name(0)} (sm_{major}{minor}) is not supported by this torch build "
                f"({torch.__version__}, kernels for {', '.join(torch.cuda.get_arch_list())}); falling back to CPU. "
                "For Blackwell (sm_120) and newer GPUs create the environment from environment_rna004_cu128.yaml "
                "(see README), or pass --force_gpu to try anyway.\n")
            return torch.device("cpu")
        return torch.device("cuda")
    return torch.device("cpu")


def load_model(checkpoint_path, device):
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    cfg = ckpt["model_config"]
    model = ReDDModel(**cfg).float().to(device)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    return model, cfg, int(ckpt.get("output_head", 0))


def predict_file(h5path, model, device, out, window_size, featuredim, head, batch_size, verbose):
    with h5py.File(h5path, "r") as h5f:
        total = h5f["X"].shape[0]
        if total == 0:
            if verbose:
                print(f"{h5path}: no sites, skipping")
            return 0
        extracted_window = h5f["X"].shape[1]
        if extracted_window < window_size:
            raise ValueError(f"{h5path}: extracted window {extracted_window} < model window {window_size}")
        center = extracted_window // 2
        nt = window_size // 2
        n_written = 0
        for start in range(0, total, batch_size):
            stop = min(start + batch_size, total)
            x = h5f["X"][start:stop, center - nt: center + nt + 1, :featuredim]
            x = np.nan_to_num(x, copy=False)
            info = h5f["info"][start:stop]
            with torch.no_grad():
                _, _, pred_mod = model(torch.from_numpy(x).to(device))
                prob = torch.sigmoid(pred_mod[:, head]).cpu().numpy()
            lines = []
            for info_str, p in zip(info, prob):
                if isinstance(info_str, bytes):
                    info_str = info_str.decode("utf-8")
                lines.append(f"{info_str}\t{p:.6f}\n")
            out.write("".join(lines))
            n_written += stop - start
            if verbose:
                print(f"{h5path}: {stop}/{total} ({stop / total:.1%})", flush=True)
        return n_written


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--checkpoint", required=True, help="model checkpoint (.pt)")
    parser.add_argument("--input", nargs="+", required=True, help="one or more ReDD HDF5 feature files")
    parser.add_argument("--output", required=True, help="output text file (.gz for compressed)")
    parser.add_argument("--device", default="CPU", help="CPU or GPU [CPU]")
    parser.add_argument("--batch_size", type=int, default=4096)
    parser.add_argument("--threads", type=int, default=1, help="torch CPU threads")
    parser.add_argument("--allow_missing", action="store_true", help="skip input files that do not exist")
    parser.add_argument("--force_gpu", action="store_true", help="use the GPU even if this torch build has no kernels for its architecture")
    parser.add_argument("--verbose", type=int, default=1)
    args = parser.parse_args()

    torch.set_num_threads(max(1, args.threads))
    device = resolve_device(args.device, args.force_gpu)
    t0 = time.time()
    model, cfg, head = load_model(args.checkpoint, device)
    window_size, featuredim = cfg["window_size"], cfg["featuredim"]
    if args.verbose:
        print(f"model loaded on {device} (window={window_size}, features={featuredim}, head={head}) in {time.time() - t0:.1f}s")

    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)
    opener = gzip.open if args.output.endswith(".gz") else open
    n_total = 0
    with opener(args.output, "wt") as out:
        for h5path in args.input:
            if not os.path.exists(h5path):
                if args.allow_missing:
                    sys.stderr.write(f"WARNING: {h5path} not found, skipping\n")
                    continue
                raise FileNotFoundError(h5path)
            n_total += predict_file(h5path, model, device, out, window_size, featuredim, head,
                                    args.batch_size, args.verbose)
    if args.verbose:
        print(f"wrote {n_total} predictions to {args.output} in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
