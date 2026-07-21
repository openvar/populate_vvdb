import json
import logging

import VariantValidator
import vvhgvs

import logger as log_set


def add_missing_ensembl_from_tark(accession, validator, logger):
    """
    Add a minimal transcript_info record for an Ensembl transcript that
    cannot be populated through Ensembl REST.

    TARK is used as the source of transcript/version, Ensembl gene and
    MANE data. VVTA is used to establish the transcript gene symbol.

    Information that cannot be established from TARK, such as Ensembl
    canonical status and CCDS, is left unset.
    """

    # ------------------------------------------------------------------
    # Get exact transcript version from TARK
    # ------------------------------------------------------------------
    tark_record = VariantValidator.modules.utils.ensembl_tark(
        id=accession,
        endpoint="/api/transcript/stable_id_with_version/",
        options="expand=transcript_release_set,genes"
    )

    tark_json = tark_record["record"]

    if not tark_json or not tark_json.get("results"):
        logger.error(
            "Cannot create transcript_info record from TARK for %s",
            accession
        )
        return False

    record = tark_json["results"][0]

    # ------------------------------------------------------------------
    # Confirm TARK returned the requested accession.version
    # ------------------------------------------------------------------
    version = "%s.%s" % (
        record["stable_id"],
        record["stable_id_version"]
    )

    if version != accession:
        logger.error(
            "TARK version mismatch for %s: returned %s",
            accession,
            version
        )
        return False

    # ------------------------------------------------------------------
    # Get information from VVTA
    # ------------------------------------------------------------------
    try:
        uta_info = validator.hdp.get_tx_identity_info(version)
    except vvhgvs.exceptions.HGVSDataNotAvailableError:
        logger.error(
            "Cannot create transcript_info record for %s: "
            "transcript not found in VVTA",
            accession
        )
        return False

    uta_symbol = str(uta_info[6])

    if uta_symbol == "":
        uta_symbol = "unassigned"

    # Use the VVTA gene symbol rather than relying on the TARK gene name.
    genbank_symbol = uta_symbol

    # ------------------------------------------------------------------
    # Gene information from TARK
    # ------------------------------------------------------------------
    genes = record.get("genes", [])

    if not genes:
        logger.error(
            "No gene information returned by TARK for %s",
            accession
        )
        return False

    # Multiple entries may represent versions of the same ENSG.
    # Use the most recent gene version.
    gene = max(
        genes,
        key=lambda x: x.get("stable_id_version", 0)
    )

    ensemblgene_id = gene.get("stable_id")

    if not ensemblgene_id:
        logger.error(
            "No Ensembl gene ID returned by TARK for %s",
            accession
        )
        return False

    # TARK does not provide the Ensembl transcript display_name.
    # Use the VVTA gene symbol as a conservative description.
    description = genbank_symbol

    # ------------------------------------------------------------------
    # MANE status
    # ------------------------------------------------------------------
    mane_select = False
    mane_plus_clinical = False
    select_tx = False

    mane_type = record.get("mane_transcript_type")

    if mane_type == "MANE SELECT":
        mane_select = True
        select_tx = "MANE"

    elif mane_type is not None:
        mane_plus_clinical = True
        select_tx = "MANE_Plus_Clinical"

    # ------------------------------------------------------------------
    # TARK location information
    # ------------------------------------------------------------------
    mapped_chr = record.get("loc_region")

    if (
        mapped_chr is not None
        and record.get("loc_start") is not None
        and record.get("loc_end") is not None
    ):
        map_position = "chr%s:%s:%s" % (
            mapped_chr,
            record["loc_start"],
            record["loc_end"]
        )
    else:
        map_position = None

    # ------------------------------------------------------------------
    # Compile minimal metadata dictionary
    # ------------------------------------------------------------------
    variant = {
        "db_xref": {
            "ensemblgene": ensemblgene_id,
            "ncbigene": None,
            "HGNC": None,
            "CCDS": None,
            "select": select_tx
        },
        "chromosome": mapped_chr,
        "map": map_position,
        "note": None,
        "variant": "0",
        "mane_select": mane_select,
        "mane_plus_clinical": mane_plus_clinical,
        "ensembl_select": False,
        "refseq_select": False
    }

    # ------------------------------------------------------------------
    # Update gene information from HGNC
    # ------------------------------------------------------------------
    try:
        hgnc_data = validator.db.update_gene_stable_identifiers(
            genbank_symbol
        )
    except Exception as e:
        logger.warning(
            "Unable to retrieve HGNC information for %s (%s): %s",
            accession,
            genbank_symbol,
            e
        )
        hgnc_data = None

    if hgnc_data is not None:
        variant["map"] = hgnc_data["map_loc"]
        variant["note"] = hgnc_data["gene_name"]
        variant["db_xref"]["HGNC"] = hgnc_data["hgnc_id"]

    # Match normal transcript_info JSON structure
    variant["db_xref"]["hgnc"] = variant["db_xref"].pop("HGNC")

    # ------------------------------------------------------------------
    # Store JSON
    # ------------------------------------------------------------------
    variant = json.dumps(variant)

    query_info = [
        version,
        description,
        variant,
        version,
        genbank_symbol,
        uta_symbol
    ]

    table = "transcript_info"

    returned_data = validator.db.in_entries(
        version,
        table
    )

    if "none" in returned_data:
        validator.db.insert(
            version,
            query_info,
            table
        )

        logger.info(
            "Added missing transcript_info record from TARK for %s",
            version
        )

        return True

    logger.info(
        "Transcript_info record already exists for %s - "
        "skipping TARK fallback",
        version
    )

    return True


def populate_missing_ensembl():
    """
    Read uta_transcripts_testing.txt and add missing ENST transcript_info
    records using TARK.
    """

    log_set.configure_logging(
        console_level="INFO",
        file_level="ERROR",
    )

    logger = logging.getLogger("populate_vdb")
    validator = VariantValidator.Validator()

    with open(
        "uta_transcripts_testing.txt",
        "r"
    ) as transcript_file:

        for line in transcript_file:
            line = line.strip()

            if not line:
                continue

            # First whitespace/tab-delimited element is the accession
            accession = line.split()[0].strip()

            # Only process Ensembl transcripts
            if not accession.startswith("ENST"):
                continue

            # Skip transcripts already present in transcript_info
            returned_data = validator.db.in_entries(
                accession,
                "transcript_info"
            )

            if "none" not in returned_data:
                logger.info(
                    "Transcript_info record already exists for %s - skipping",
                    accession
                )
                continue

            logger.info(
                "No transcript_info record for %s - "
                "attempting TARK fallback",
                accession
            )

            try:
                add_missing_ensembl_from_tark(
                    accession,
                    validator,
                    logger
                )
            except Exception:
                logger.exception(
                    "Unexpected error processing %s",
                    accession
                )


if __name__ == "__main__":
    populate_missing_ensembl()