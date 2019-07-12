# populate_vvdb
Code used to populate the validator database used by VV

# Extract Transcripts from UTA
<Needs to be figured out again!>

## To run
export UTA_DB_URL="postgresql://uta_admin:uta_admin@127.0.0.1/uta/uta_20180821"
export HGVS_SEQREPO_DIR="/Users/pjf9/variant_validator_data/seqrepo/2018-08-21"
python populate.py

***Check update_log.txt and if necessary add records manually***
