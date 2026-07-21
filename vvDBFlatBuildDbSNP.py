#!/usr/bin/env python3
"""
vvDBFlatBuildDbSNP.py

Build bgzipped, indexed dbSNP flat files directly from dbSNP BCF files.

rsID index stores BGZF *block offsets* (not decompressed byte offsets),
using a binary, delta-encoded, block-compressed format to reduce size
and allow direct random-access lookup.
"""

import os
import subprocess
import pysam
import struct
import zstandard as zstd


RSID_BLOCK_SIZE = 4096   # entries per index block
ZSTD_LEVEL = 1           # <= ~3% performance impact


class FlatDbSNPBuilder:
    def __init__(self, out_dir="dbsnp_flat"):
        self.out_dir = os.path.abspath(out_dir)
        os.makedirs(self.out_dir, exist_ok=True)

    def _extract_version_tag(self, path):
        """
        vvdbSNP_GRCh37_b157_20241205.bcf → b157_20241205
        """
        name = os.path.basename(path).rsplit(".", 1)[0]
        return "_".join(name.split("_")[-2:])

    def _list_contigs(self, bcf_path):
        """
        Get contigs in file order from BCF index.
        """
        p = subprocess.run(
            ["bcftools", "index", "-s", bcf_path],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        return [line.split("\t")[0] for line in p.stdout.splitlines()]

    def build(self, assembly, bcf_path):
        version = self._extract_version_tag(bcf_path)

        version_dir = os.path.join(self.out_dir, version)
        os.makedirs(version_dir, exist_ok=True)

        out_tsv = os.path.join(
            version_dir,
            f"dbsnp_{assembly}_{version}.variants.tsv"
        )
        out_bgz = out_tsv + ".gz"
        out_tbi = out_bgz + ".tbi"

        rsid_idx_bin = os.path.join(
            version_dir,
            f"dbsnp_{assembly}_{version}.rsid.idx.zst"
        )
        rsid_idx_meta = os.path.join(
            version_dir,
            f"dbsnp_{assembly}_{version}.rsid.idx.meta"
        )

        # ---------------------------------------------------------
        # Skip if already completed
        # ---------------------------------------------------------
        if (
            os.path.exists(out_bgz)
            and os.path.exists(out_tbi)
            and os.path.exists(rsid_idx_bin)
            and os.path.exists(rsid_idx_meta)
        ):
            print(f"✔ Flat dbSNP already exists for {assembly}, skipping")
            return

        print(f"Building flat dbSNP for {assembly}")
        print(f"Version tag: {version}")
        print(f"Output directory: {version_dir}")
        print(f"Output file: {out_bgz}")

        # ---------------------------------------------------------
        # Clean partial outputs
        # ---------------------------------------------------------
        for path in (
            out_tsv,
            out_bgz,
            out_tbi,
            rsid_idx_bin,
            rsid_idx_meta,
        ):
            if os.path.exists(path):
                os.unlink(path)

        contigs = self._list_contigs(bcf_path)

        # ---------------------------------------------------------
        # Process chromosome-by-chromosome
        # ---------------------------------------------------------
        for contig in contigs:
            print(f"  → processing {contig}")

            cmd = (
                f"bcftools view -r {contig} {bcf_path} | "
                f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\n' | "
                f"awk '{{OFS=\"\\t\"; "
                f"split($4,a,\",\"); "
                f"for (i in a) if (a[i] !~ /^</) "
                f"print $1, $2, $2+length($3)-1, substr($5,3), $3, a[i]}}' | "
                f"LC_ALL=C sort -k2,2n >> {out_tsv}"
            )

            subprocess.run(
                ["bash", "-o", "pipefail", "-c", cmd],
                check=True,
            )

        # ---------------------------------------------------------
        # bgzip + tabix
        # ---------------------------------------------------------
        subprocess.run(["bgzip", "-f", out_tsv], check=True)

        subprocess.run(
            ["tabix", "-f", "-s", "1", "-b", "2", "-e", "3", out_bgz],
            check=True,
        )

        # ---------------------------------------------------------
        # rsID → BGZF block-offset index (compressed binary)
        # ---------------------------------------------------------
        print("Building rsID compressed binary index")

        compressor = zstd.ZstdCompressor(level=ZSTD_LEVEL)

        with (
            pysam.BGZFile(out_bgz, "rb") as bgz,
            open(rsid_idx_bin, "wb") as bin_out,
            open(rsid_idx_meta, "w") as meta_out,
        ):
            last_rsid = 0
            buf = []
            out_offset = 0
            first_rsid_in_block = None

            while True:
                vfo = bgz.tell()
                block_offset = vfo >> 16

                line = bgz.readline()
                if not line:
                    break

                fields = line.rstrip(b"\n").split(b"\t")
                rsid = int(fields[3])

                if rsid == last_rsid:
                    continue

                # -------------------------------------------------
                # Block start: reset rsID baseline (OPTION 2)
                # -------------------------------------------------
                if first_rsid_in_block is None:
                    first_rsid_in_block = rsid
                    meta_out.write(f"{rsid}\t{out_offset}\n")
                    rsid_delta = 0
                else:
                    rsid_delta = rsid - last_rsid

                buf.append(struct.pack("<qQ", rsid_delta, block_offset))

                last_rsid = rsid

                if len(buf) == RSID_BLOCK_SIZE:
                    raw = b"".join(buf)
                    comp = compressor.compress(raw)
                    bin_out.write(comp)
                    out_offset += len(comp)
                    buf.clear()
                    first_rsid_in_block = None

            if buf:
                raw = b"".join(buf)
                comp = compressor.compress(raw)
                bin_out.write(comp)

        print("✔ Flat dbSNP build complete")


# ---------------------------------------------------------
# CLI
# ---------------------------------------------------------
if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Build flat dbSNP files from existing dbSNP BCFs"
    )

    parser.add_argument(
        "--out",
        default="dbsnp_flat",
        help="Output directory (default: dbsnp_flat)",
    )

    parser.add_argument(
        "--grch37",
        help="Path to GRCh37 dbSNP BCF",
    )

    parser.add_argument(
        "--grch38",
        help="Path to GRCh38 dbSNP BCF",
    )

    args = parser.parse_args()

    if not args.grch37 and not args.grch38:
        parser.error("At least one of --grch37 or --grch38 must be specified")

    builder = FlatDbSNPBuilder(out_dir=args.out)

    if args.grch37:
        if not os.path.exists(args.grch37):
            sys.exit(f"ERROR: GRCh37 BCF not found: {args.grch37}")
        builder.build("grch37", args.grch37)

    if args.grch38:
        if not os.path.exists(args.grch38):
            sys.exit(f"ERROR: GRCh38 BCF not found: {args.grch38}")
        builder.build("grch38", args.grch38)