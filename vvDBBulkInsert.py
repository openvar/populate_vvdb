#!/usr/bin/env python3
"""
vvDBBulkInsertDbSNP.py
----------------------

Robust dbSNP bulk loader for VariantValidator.

Guarantees:
  - No FIFO
  - No LOAD DATA LOCAL INFILE
  - Crash-safe, restart-safe (IGNORE semantics)
  - Deterministic chunk files on disk
  - Immediate filesystem proof of correct directory
"""

import os
import subprocess
import mysql.connector
import uuid


print("### VariantValidator dbSNP CHUNKED LOADER ACTIVE ###")


class BulkInsertDbSNP:
    def __init__(self, db_config, mysql_infile_dir="dbsnp_chunks", chunk_mb=1024):
        """
        :param db_config: VariantValidator DB config dict
        :param mysql_infile_dir: directory MySQL can read via LOAD DATA INFILE
        :param chunk_mb: max size of decompressed TSV per chunk
        """

        # ------------------------------------------------------------------
        # Resolve chunk directory next to THIS script (never cwd-dependent)
        # ------------------------------------------------------------------
        base_dir = os.path.dirname(os.path.abspath(__file__))
        self.mysql_infile_dir = os.path.join(base_dir, mysql_infile_dir)
        os.makedirs(self.mysql_infile_dir, exist_ok=True)

        print(f"Using MySQL INFILE directory: {self.mysql_infile_dir}")

        # ------------------------------------------------------------------
        # HARD filesystem probe — must succeed or we abort immediately
        # ------------------------------------------------------------------
        probe_path = os.path.join(self.mysql_infile_dir, ".vv_loader_ok")
        try:
            with open(probe_path, "w") as fh:
                fh.write("ok\n")
            os.sync()
        except Exception as e:
            raise RuntimeError(
                f"INFILE directory is not usable: {self.mysql_infile_dir}"
            ) from e

        print(f"✅ INFILE directory writable, probe file created: {probe_path}")

        # ------------------------------------------------------------------
        # MySQL config
        # ------------------------------------------------------------------
        self.db_config = {
            "user": db_config["user"],
            "password": db_config["password"],
            "host": db_config["host"],
            "port": int(db_config["port"]),
            "database": db_config["database"],
            "autocommit": False,
            "allow_local_infile": False,
        }

        self.chunk_bytes = chunk_mb * 1024 * 1024

    # ---------------------------------------------------------
    # MySQL connection
    # ---------------------------------------------------------
    def get_conn(self):
        return mysql.connector.connect(**self.db_config)

    # ---------------------------------------------------------
    # Index handling
    # ---------------------------------------------------------
    def drop_indexes(self, table):
        conn = self.get_conn()
        cur = conn.cursor()

        print(f"Dropping secondary indexes on {table} ...")

        for idx in ("lookup_rsid", "lookup_contig_pos"):
            try:
                cur.execute(f"ALTER TABLE {table} DROP INDEX {idx}")
                print(f"  Dropped {idx}")
            except Exception:
                print(f"  {idx} not present; skipping")

        conn.commit()
        cur.close()
        conn.close()

    def add_indexes(self, table):
        conn = self.get_conn()
        cur = conn.cursor()

        print(f"Rebuilding secondary indexes on {table} ...")

        cur.execute(
            f"""
            ALTER TABLE {table}
              ADD INDEX lookup_contig_pos (contig_id, pos),
              ADD INDEX lookup_rsid (rs_int)
            """
        )

        conn.commit()
        cur.close()
        conn.close()

    # ---------------------------------------------------------
    # Core LOAD DATA
    # ---------------------------------------------------------
    def _load_chunk(self, conn, table, tsv_path):
        sql = (
            f"LOAD DATA INFILE '{tsv_path}' "
            f"IGNORE INTO TABLE {table} "
            "FIELDS TERMINATED BY '\\t' "
            "LINES TERMINATED BY '\\n' "
            "(@rs, @contig, @pos, @ref, @alt) "
            "SET "
            "  rs_int    = @rs, "
            "  contig_id = @contig, "
            "  pos       = @pos, "
            "  ref_blob  = COMPRESS(@ref), "
            "  alt_blob  = COMPRESS(@alt)"
        )

        cur = conn.cursor()
        cur.execute(sql)
        conn.commit()
        cur.close()

    # ---------------------------------------------------------
    # Chunked gunzip loader — deterministic
    # ---------------------------------------------------------
    def load_tsv_gz(self, table, gz_path):
        if not os.path.exists(gz_path):
            raise FileNotFoundError(gz_path)

        print(f"Loading {gz_path} into {table}")
        print(f"Chunk size: {self.chunk_bytes // (1024 * 1024)} MB")

        conn = self.get_conn()

        gunzip = subprocess.Popen(
            ["gunzip", "-c", gz_path],
            stdout=subprocess.PIPE,
            text=True,
        )

        try:
            eof = False
            chunk_no = 0

            while not eof:
                chunk_no += 1
                chunk_name = f"chunk_{chunk_no}_{uuid.uuid4().hex}.tsv"
                chunk_path = os.path.join(self.mysql_infile_dir, chunk_name)

                print(f"Creating chunk file: {chunk_path}")

                written = 0
                with open(chunk_path, "w") as out:
                    while written < self.chunk_bytes:
                        line = gunzip.stdout.readline()
                        if not line:
                            eof = True
                            break
                        out.write(line)
                        written += len(line)

                if written == 0:
                    os.unlink(chunk_path)
                    break

                print(
                    f"  → loading chunk {chunk_no} "
                    f"({written // (1024 * 1024)} MB)"
                )

                self._load_chunk(conn, table, chunk_path)
                os.unlink(chunk_path)

        finally:
            gunzip.wait()
            conn.close()

        print(f"Completed load for {table}")

    # ---------------------------------------------------------
    # Public entry point
    # ---------------------------------------------------------
    def bulk_insert_table(self, table, tsv_gz_path):
        print(f"\n=== dbSNP load start: {table} ===\n")

        self.drop_indexes(table)

        if not tsv_gz_path.endswith(".gz"):
            raise ValueError("Expected a .tsv.gz input")

        self.load_tsv_gz(table, tsv_gz_path)
        self.add_indexes(table)

        print(f"\n=== dbSNP load complete: {table} ===\n")


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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>
