import os
import re
import sys
import traceback
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

        print('Updating - ' + tx_id)
        accession = tx_id
        # Look for the accession in our database
        # Connect to database and send request
        if 'ENST' not in tx_id:
            try:
                vval.update_transcript_record(tx_id)
                print('True')
            except BaseException as e:
                fo.write(tx_id + '\n')
        else:
            try:
                vval.update_transcript_record(tx_id, genome_build='GRCh38', test=True)
                print('True')
                exit()
            except VariantValidator.modules.utils.DatabaseConnectionError:
                try:
                    vval.update_transcript_record(tx_id, genome_build='GRCh37', test=True)
                    print('True')
                    exit()
                except BaseException as e:
                    exc_type, exc_value, last_traceback = sys.exc_info()
                    print(str(exc_type) + " " + str(exc_value))
                    traceback.print_tb(last_traceback)
                    print('False')
                    exit()
                    fo.write(tx_id + '\n')
            except BaseException as e:
                exc_type, exc_value, last_traceback = sys.exc_info()
                print(str(exc_type) + " " + str(exc_value))
                traceback.print_tb(last_traceback)
                print('False')
                exit()
                fo.write(tx_id + '\n')


fo.close()

# Update everything else
vval.update_vv_data()
