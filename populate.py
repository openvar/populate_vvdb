import os
import re
import sys
import traceback
import VariantValidator
import VariantValidator.variantanalyser.dbControls as va_dbCrl
import VariantValidator.variantanalyser.batch as va_btch
import VariantValidator.variantValidator as vv
import hgvs
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper

"""
Notes for error fixing. 
1. Increase transcriptVariant max cell content to 50

2. Increase hgncSymbol to 50 (??? NR_036569.1) and utaSymbol to 50 (NR_036570.1 ???)

2. local variable 'version' referenced before assignment

3. local variable 'hgnc_symbol' referenced before assignment
"""

seqrepo_current_version = '2018-08-21'
HGVS_SEQREPO_DIR = '/Users/pjf9/variant_validator_data/seqrepo/' + seqrepo_current_version
os.environ['HGVS_SEQREPO_DIR'] = HGVS_SEQREPO_DIR
uta_current_version = 'uta_20180821'
UTA_DB_URL = 'postgresql://uta_admin:uta_admin@127.0.0.1/uta/' + uta_current_version
os.environ['UTA_DB_URL'] = UTA_DB_URL
os.environ['PYLIFTOVER_DIR'] = '/Users/pjf9/variant_validator_data/pyLiftover/'

hdp = hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()

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
        else:
            alt_aln_method = 'splign'

        evm = hgvs.assemblymapper.AssemblyMapper(hdp,
                                                 assembly_name='GRCh38',
                                                 alt_aln_method=alt_aln_method,
                                                 normalize=True,
                                                 replace_reference=True
                                                 )
        print('Updating - ' + tx_id)
        # try:
        #     va_dbCrl.data.update_transcript_info_record(tx_id, hdp)
        # except Exception as e:
        #     fo.write(tx_id + ' : ' + str(e) + '\n')

        accession = tx_id
        # Look for the accession in our database
        # Connect to database and send request
        entry = va_dbCrl.data.in_entries(accession, 'transcript_info')

        # Analyse the returned data and take the necessary actions
        # If the error key exists
        if 'error' in entry:
            # Open a hgvs exception log file in append mode
            error = entry['description']
            fo.write(tx_id + ' : ' + str(error) + '\n\n')
            print(error)
            continue

        # If the accession key is found
        elif 'accession' in entry:
            description = entry['description']
            # If the current entry is too old
            if entry['expiry'] == 'true':
                dbaction = 'update'
                try:
                    entry = va_btch.data_add(input=input, alt_aln_method=alt_aln_method,
                                             accession=accession, dbaction=dbaction, hp=hp, evm=evm,
                                             hdp=hdp)
                except hgvs.exceptions.HGVSError as e:
                    exc_type, exc_value, last_traceback = sys.exc_info()
                    te = traceback.format_exc()
                    tbk = [str(exc_type), str(exc_value), str(te)]
                    error = str('\n'.join(tbk))
                    fo.write(tx_id + ' : ' + str(error) + '\n\n')
                    print(error)
                    continue
                except Exception as e:
                    exc_type, exc_value, last_traceback = sys.exc_info()
                    te = traceback.format_exc()
                    tbk = [str(exc_type), str(exc_value), str(te)]
                    error = str('\n'.join(tbk))
                    fo.write(tx_id + ' : ' + str(error) + '\n\n')
                    print(error)
                    continue
                else:
                    print('updated')
            else:
                print('up-to-date')
        # If the none key is found add the description to the database
        elif 'none' in entry:
            dbaction = 'insert'
            try:
                entry = va_btch.data_add(input=input, alt_aln_method=alt_aln_method,
                                         accession=accession, dbaction=dbaction, hp=hp, evm=evm,
                                         hdp=hdp)
            except Exception as e:
                exc_type, exc_value, last_traceback = sys.exc_info()
                te = traceback.format_exc()
                tbk = [str(exc_type), str(exc_value), str(te)]
                error = str('\n'.join(tbk))
                fo.write(tx_id + ' : ' + str(error) + '\n\n')
                print(error)
                continue
            else:
                print('inserted')

        # If no correct keys are found
        else:
            # Open a hgvs exception log file in append mode
            error = 'Unknown error type'
            fo.write(tx_id + ' : ' + str(error) + '\n\n')
            continue

fo.close()

# Update everything else
vv.update_vv_data()
