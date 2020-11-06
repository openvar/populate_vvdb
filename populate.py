import os
import re
import VariantValidator
vval = VariantValidator.Validator()


"""
Notes for error fixing. 
1. Increase transcriptVariant max cell content to 50

2. Increase hgncSymbol to 50 (??? NR_036569.1) and utaSymbol to 50 (NR_036570.1 ???)

2. local variable 'version' referenced before assignment

3. local variable 'hgnc_symbol' referenced before assignment
"""
ROOT = os.path.dirname(os.path.abspath(__file__))
infile = ROOT

# Use pgAdmin4 to download the transcripts table from UTA then replce commas for tabs and remove " characters
infile = os.path.join(ROOT, 'uta_transcripts.txt')
logfile = os.path.join(ROOT, 'update_log.txt')
fo = open(logfile, "w")

# open file and loop through to populate.py
with open(infile) as tx_data:
    for tx_line in tx_data:
        tx_cells = tx_line.split('\t')
        tx_id = tx_cells[0]

        # populate
        if tx_id == 'ac':
            continue
        if not re.search('.', tx_id):
            continue
        if 'ENST' in tx_id:
            continue

        print('Updating - ' + tx_id)
        accession = tx_id
        # Look for the accession in our database
        # Connect to database and send request
        try:
            vval.update_transcript_record(tx_id)
            print('True')
        except BaseException:
            print('False')
            fo.write(tx_id + '\n')

fo.close()

# Update everything else
vv.update_vv_data()
