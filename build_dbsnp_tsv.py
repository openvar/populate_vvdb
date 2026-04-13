#!/usr/bin/env python3

import os
import subprocess
import mysql.connector
import gzip


# ---------------------------------------------------------
# MySQL helpers
# ---------------------------------------------------------

def _get_mysql_conn(db_config):
    """Standalone binary-safe MySQL connection."""
    return mysql.connector.connect(
        user=db_config["user"],
        password=db_config["password"],
        host=db_config["host"],
        port=int(db_config["port"]),
        database=db_config["database"],
        autocommit=True,
    )


def _load_contig_cache(conn, assembly):
    """Load contig_lookup_<assembly> into an in-memory dict."""
    table = f"contig_lookup_{assembly}"
    cur = conn.cursor()
    cur.execute(f"SELECT contig_name, contig_id FROM {table}")
    cache = dict(cur.fetchall())
    cur.close()
    return cache


def _insert_contig(conn, assembly, contig_name):
    """Insert missing contig and return contig_id."""
    table = f"contig_lookup_{assembly}"
    cur = conn.cursor()
    cur.execute(
        f"INSERT INTO {table} (contig_name) VALUES (%s)",
        (contig_name,),
    )
    contig_id = cur.lastrowid
    cur.close()
    return contig_id


# ---------------------------------------------------------
# Main TSV builder
# ---------------------------------------------------------

def make_tsv_from_bcf(bcf_path, tsv_out, db_config=None, assembly=None):
    """
    Create loadable TSV.GZ from dbSNP BCF.

    OUTPUT FORMAT (tab-separated, gzipped):
        rs_int  contig_id  pos  ref  alt

    DESIGN NOTES:
        - REF / ALT stored as plain text
        - Compression happens in MySQL via COMPRESS()
        - Reference allele (REF == ALT) is implicit at query time
        - Only non-reference ALT alleles are emitted
    """

    if not os.path.exists(bcf_path):
        raise FileNotFoundError(f"BCF file not found: {bcf_path}")

    if db_config is None:
        raise ValueError("db_config must be provided")

    if assembly not in ("grch37", "grch38"):
        raise ValueError("assembly must be 'grch37' or 'grch38'")

    gz_out = tsv_out + ".gz"

    if os.path.exists(gz_out):
        print(f"⏩ Loadable TSV already exists, skipping: {gz_out}")
        return gz_out

    print(f"Creating loadable TSV from {bcf_path}")
    print(f"Output: {gz_out}")

    # rsID, chrom, pos, ref, alt(s)
    fmt = r"%ID\t%CHROM\t%POS\t%REF\t%ALT\n"
    cmd = ["bcftools", "query", "-f", fmt, bcf_path]

    conn = _get_mysql_conn(db_config)
    contig_cache = _load_contig_cache(conn, assembly)

    BUFFER_SIZE = 10000
    buffer = []

    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc, \
         gzip.open(gz_out, "wt", compresslevel=1) as out:

        for line in proc.stdout:
            rsid, contig_name, pos, ref, alt_list = line.rstrip("\n").split("\t")

            # Skip records without rsID
            if rsid == ".":
                continue

            # Resolve contig_id
            contig_id = contig_cache.get(contig_name)
            if contig_id is None:
                contig_id = _insert_contig(conn, assembly, contig_name)
                contig_cache[contig_name] = contig_id

            rs_int = int(rsid[2:])
            alts = alt_list.split(",")

            # Emit one row per non-reference ALT allele
            for alt in alts:
                # Skip symbolic alleles (e.g. <DEL>, <INS>)
                if alt.startswith("<") and alt.endswith(">"):
                    continue

                buffer.append(
                    f"{rs_int}\t{contig_id}\t{pos}\t{ref}\t{alt}\n"
                )

            if len(buffer) >= BUFFER_SIZE:
                out.writelines(buffer)
                buffer.clear()

        if buffer:
            out.writelines(buffer)

    conn.close()

    print(f"✔ Loadable TSV created: {gz_out}")
    return gz_out

# <LICENSE>
# Copyright (C) 2016-2026 VariantValidator Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
