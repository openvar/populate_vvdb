
#!/usr/bin/env python3
"""
populate.py
-----------

Full dbSNP build → TSV(.gz) generation → MySQL population pipeline
for VariantValidator.

NOTE:
- dbSNP load runs and then EXITS deliberately.
- Transcript / version updates are intentionally not run here.
"""

import os
import sys
import re
import subprocess
import traceback
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import freeze_support

import VariantValidator
from VariantValidator import update_vv_db
import add_version_info

from build_dbsnp_tsv import make_tsv_from_bcf
from vvDBBulkInsert import BulkInsertDbSNP


# ---------------------------------------------------------
# SETUP
# ---------------------------------------------------------
vval = VariantValidator.Validator()
ROOT = os.path.dirname(os.path.abspath(__file__))
LOGFILE = os.path.join(ROOT, "update_log.txt")
fo = open(LOGFILE, "w")


# ---------------------------------------------------------
# UTILITY: Auto-detect TSV / TSV.GZ
# ---------------------------------------------------------
def find_tsv(base, stem):
    gz = os.path.join(base, stem + ".load.tsv.gz")
    raw = os.path.join(base, stem + ".load.tsv")

    if os.path.exists(gz):
        print(f"✔ Using compressed TSV: {gz}")
        return gz

    if os.path.exists(raw):
        print(f"✔ Using uncompressed TSV: {raw}")
        return raw

    raise FileNotFoundError(f"No TSV or TSV.GZ found for {stem}")


# ---------------------------------------------------------
# STEP 1: Run dbSNP pre-build
# ---------------------------------------------------------
def run_dbsnp_prebuild():
    script_path = os.path.join(ROOT, "scripts", "create_VVdbSNP.sh")

    print("\n==== Running dbSNP Pre-Build ====\n")

    if not os.path.exists(script_path):
        print(f"ERROR: Script not found: {script_path}")
        sys.exit(1)

    subprocess.run(["bash", script_path], check=True)
    print("✔ dbSNP pre-build completed.\n")


# ---------------------------------------------------------
# STEP 2: Locate latest BCFs
# ---------------------------------------------------------
def get_latest_bcf(assembly):
    for fname in os.listdir(ROOT):
        if fname.startswith(f"vvdbSNP_{assembly}") and fname.endswith(".bcf"):
            return os.path.join(ROOT, fname)

    print(f"ERROR: No BCF found for {assembly}")
    sys.exit(1)


# ---------------------------------------------------------
# STEP 3: Build TSVs (skip if present)
# ---------------------------------------------------------
def build_tsvs():

    print("==== Building dbSNP TSVs ====\n")

    bcf37 = get_latest_bcf("GRCh37")
    bcf38 = get_latest_bcf("GRCh38")

    tsv37 = os.path.join(ROOT, "dbsnp_grch37.load.tsv")
    tsv38 = os.path.join(ROOT, "dbsnp_grch38.load.tsv")

    make_tsv_from_bcf(
        bcf37,
        tsv37,
        db_config=vval.dbConfig,
        assembly="grch37"
    )

    make_tsv_from_bcf(
        bcf38,
        tsv38,
        db_config=vval.dbConfig,
        assembly="grch38"
    )


# ---------------------------------------------------------
# STEP 4: Load dbSNP tables (FIXED)
# ---------------------------------------------------------
def load_dbsnp():

    print("==== Loading dbSNP tables ====\n")

    tsv37 = find_tsv(ROOT, "dbsnp_grch37")
    tsv38 = find_tsv(ROOT, "dbsnp_grch38")

    loader = BulkInsertDbSNP(vval.dbConfig)

    loader.bulk_insert_table("dbsnp_grch37", tsv37)
    loader.bulk_insert_table("dbsnp_grch38", tsv38)

    print("\n==== dbSNP build complete ====\n")


# ---------------------------------------------------------
# MAIN WORKFLOW
# ---------------------------------------------------------
def main():

    # dbSNP pipeline ONLY
    run_dbsnp_prebuild()
    build_tsvs()
    load_dbsnp()

    # HARD STOP — intentional
    sys.exit(0)

    # -----------------------------------------------------
    # Everything below here is intentionally unreachable
    # during dbSNP builds.
    # -----------------------------------------------------

    testing, transcript_set, infile = check_args()

    with open(infile) as tx_data, ProcessPoolExecutor(max_workers=4) as executor:
        tx_ids = ((line.split("\t")[0], transcript_set) for line in tx_data)
        executor.map(process_transcript, tx_ids)

    fo.close()

    if testing == "False":
        update_vv_db.update()
        add_version_info.update_version()

    print("\n===== UPDATE COMPLETE =====")


# ---------------------------------------------------------
# ENTRY POINT
# ---------------------------------------------------------
if __name__ == "__main__":
    freeze_support()
    main()


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
