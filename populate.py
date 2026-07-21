#!/usr/bin/env python3
"""
populate.py
-----------

Full dbSNP build → flat-file dbSNP build
→ (intentionally STOP) → transcript pipeline (preserved but unreachable)
"""

import os
import sys
import re
import subprocess
import logging
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import freeze_support

import VariantValidator
from VariantValidator import update_vv_db
import add_version_info
import missing_validator_ens_records

from vvDBFlatBuildDbSNP import FlatDbSNPBuilder

import logger as log_set

log_set.configure_logging(
    console_level="INFO",
    file_level="ERROR",
)

logger = logging.getLogger("populate_vdb")

# ---------------------------------------------------------
# SETUP
# ---------------------------------------------------------
vval = VariantValidator.Validator()
ROOT = os.path.dirname(os.path.abspath(__file__))
LOGFILE = os.path.join(ROOT, "update_log.txt")

logger.info(f"Setup with ROOT={ROOT} and LOGFILE={LOGFILE}")

# ---------------------------------------------------------
# BCF finder
# ---------------------------------------------------------
def find_bcf(base, assembly):
    for fname in os.listdir(base):
        if fname.startswith(f"vvdbSNP_{assembly}") and fname.endswith(".bcf"):
            path = os.path.join(base, fname)
            print(f"✔ Using BCF: {path}")
            return path
    raise FileNotFoundError(f"No BCF found for {assembly}")


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
# STEP 2: Build flat dbSNP files
# ---------------------------------------------------------
def build_flat_dbsnp():

    print("==== Building flat dbSNP files ====\n")

    bcf37 = find_bcf(ROOT, "GRCh37")
    bcf38 = find_bcf(ROOT, "GRCh38")

    builder = FlatDbSNPBuilder(out_dir=os.path.join(ROOT, "dbsnp_flat"))

    builder.build("grch37", bcf37)
    builder.build("grch38", bcf38)

    print("\n==== Flat dbSNP build complete ====\n")


# ---------------------------------------------------------
# ARGUMENT HANDLING (VERBATIM)
# ---------------------------------------------------------
def check_args():
    if len(sys.argv) != 3:
        print(
            "Too few arguments. The command required is: python populate.py False all: "
            "Where True/False dictates whether to use the uta_transcripts_testing.txt "
            "file for testing or the default input file and All is one of all, ensembl or refseq"
        )
        sys.exit(1)

    testing = sys.argv[1]
    transcript_set = sys.argv[2]
    print(testing)

    if testing == "False":
        infile = os.path.join(ROOT, "uta_transcripts.txt")
    elif testing == "True":
        infile = os.path.join(ROOT, "uta_transcripts_testing.txt")
    else:
        print(
            "Argument 1 must be True or False where use the uta_transcripts_testing.txt "
            "file is True or False dictating whether to use the uta_transcripts_testing.txt "
            "file for testing or the default input file"
        )
        sys.exit(1)

    if transcript_set not in ["all", "ensembl", "refseq"]:
        print("Argument 2 must be all, ensembl or refseq to select a specific transcript set or not")
        sys.exit(1)

    return testing, transcript_set, infile


# ---------------------------------------------------------
# WORKER INITIALISATION
# ---------------------------------------------------------
def init_worker():
    log_set.configure_logging(
        console_level="INFO",
        file_level="ERROR",
    )

# ---------------------------------------------------------
# TRANSCRIPT PROCESSING (UNCHANGED)
# ---------------------------------------------------------
def process_transcript(args):
    tx_id, transcript_set = args

    tx_id = tx_id.strip()

    if tx_id == "ac":
        return
    if not re.search(r"\.", tx_id):
        return
    if ".." in tx_id:
        return
    if tx_id.startswith("U"):
        return
    if "NG" in tx_id:
        return

    logger.debug("Updating - " + tx_id)

    #  RefSeq processing
    if "ENST" not in tx_id:
        if transcript_set == "ensembl":
            return
        try:
            vval.update_transcript_record(tx_id)
            logger.info(f"Updated transcript record for {tx_id}")
        except BaseException as e:
            logger.exception(f"Error updating transcript record for {tx_id}: {e}")
            with open(LOGFILE, "a") as fo:
                fo.write(f"{tx_id}\t{' '.join(str(e).split())}\n")
                fo.flush()
                os.fsync(fo.fileno())

    #  Ensembl processing
    else:
        if transcript_set == "refseq":
            print("Bypass")
            return
        try:
            vval.update_transcript_record(tx_id, genome_build="GRCh38", test=True)
            logger.info(f"Updated transcript record for {tx_id}")
        except VariantValidator.modules.utils.DatabaseConnectionError:
            try:
                vval.update_transcript_record(tx_id, genome_build="GRCh37", test=True)
                logger.info(f"Updated transcript record for {tx_id}")
            except BaseException as e:

                if "Cannot retrieve data from Ensembl REST for record" in str(e):
                    logger.info(
                        "Unable to populate %s from Ensembl REST - attempting TARK fallback",
                        tx_id
                    )

                    try:
                        if missing_validator_ens_records.add_missing_ensembl_from_tark(
                                tx_id,
                                vval,
                                logger
                        ):
                            logger.info(
                                "Successfully added %s using TARK fallback",
                                tx_id
                            )
                            return

                    except BaseException as tark_error:
                        logger.exception(
                            "TARK fallback failed for %s: %s",
                            tx_id,
                            tark_error
                        )

                logger.exception(
                    f"Error updating transcript record for {tx_id}: {e}"
                )

                with open(LOGFILE, "a") as fo:
                    fo.write(f"{tx_id}\t{' '.join(str(e).split())}\n")
                    fo.flush()
                    os.fsync(fo.fileno())
        except BaseException as e:
            logger.exception(f"Error updating transcript record for {tx_id}: {e}")
            with open(LOGFILE, "a") as fo:
                fo.write(f"{tx_id}\t{' '.join(str(e).split())}\n")
                fo.flush()
                os.fsync(fo.fileno())


# ---------------------------------------------------------
# MAIN WORKFLOW
# ---------------------------------------------------------
def main():

    # Create/empty the retry file
    with open(LOGFILE, "w"):
        pass

    # # ---- dbSNP pipeline ----
    # run_dbsnp_prebuild()
    # build_flat_dbsnp()
    #
    # # ---- deliberate hard stop (INTENTIONAL) ----
    # sys.exit(0)

    # ---- transcript pipeline (UNCHANGED, UNREACHABLE) ----
    testing, transcript_set, infile = check_args()

    with open(infile) as tx_data, ProcessPoolExecutor(
            max_workers=5,
            initializer=init_worker,
    ) as executor:
        tx_ids = ((tx_line.split("\t")[0], transcript_set) for tx_line in tx_data)
        list(executor.map(process_transcript, tx_ids))

    if testing == "False":
        update_vv_db.update()

    if testing == "False":
        # Re index
        db = update_vv_db.connect()
        update_vv_db.drop_core_indexes(db)
        update_vv_db.rebuild_core_indexes(db)

    if testing == "False":
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>