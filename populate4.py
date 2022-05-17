import os
import re
import VariantValidator
from VariantValidator import update_vv_db
import add_version_info
import sys
vval = VariantValidator.Validator()

"""
Notes for error fixing. 
1. Increase transcriptVariant max cell content to 50

2. Increase hgncSymbol to 50 (??? NR_036569.1) and utaSymbol to 50 (NR_036570.1 ???)

2. local variable 'version' referenced before assignment

3. local variable 'hgnc_symbol' referenced before assignment
"""
ROOT = os.path.dirname(os.path.abspath(__file__))

# Check Args
if len(sys.argv) != 2:
    print('Too few arguments. The command required is: python populate.py False : Where True/False dictates whether to '
          'use the uta_transcripts_testing.txt file for testing or the default input file')
    exit()

testing = sys.argv[1]
print("Run test transcript input file")
print(testing)

# Use pgAdmin4 to download the transcripts table from UTA then replce commas for tabs and remove " characters
if str(testing) in "False":
    infile = os.path.join(ROOT, 'uta_transcripts4.txt')
elif str(testing) in "True":
    infile = os.path.join(ROOT, 'uta_transcripts_testing.txt')
else:
    print('Argument 1 must be True or False where use the uta_transcripts_testing.txt file is True or False dictating '
          'whether to use the uta_transcripts_testing.txt file for testing or the default input file')
    exit()

print("Running file")
print(infile)

# Set log and output files
logfile = os.path.join(ROOT, 'update_log4.txt')
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
        if '..' in tx_id:
            continue
        if re.search('^U', tx_id):
            continue
        if 'NG' in tx_id:
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
                print('False')
                fo.write(tx_id + '\t' + str(e) + '\n')
                fo.flush()
                os.fsync(fo.fileno())
                print(tx_id + '\t' + str(e) + '\n')
        else:
            try:
                vval.update_transcript_record(tx_id, genome_build='GRCh38', test=True)
                print('True')
            except VariantValidator.modules.utils.DatabaseConnectionError:
                try:
                    vval.update_transcript_record(tx_id, genome_build='GRCh37', test=True)
                    print('True')
                except VariantValidator.modules.utils.DatabaseConnectionError:
                    pass
                except BaseException as e:
                    print('False')
                    fo.write(tx_id + '\t' + str(e) + '\n')
                    fo.flush()
                    os.fsync(fo.fileno())
                    print(tx_id + '\t' + str(e) + '\n')
            except BaseException as e:
                print('False')
                fo.write(tx_id + '\t' + str(e) + '\n')
                fo.flush()
                os.fsync(fo.fileno())
                print(tx_id + '\t' + str(e) + '\n')

fo.close()
# exit()

if str(testing) in "False":
    # Update everything else
    update_vv_db.update()
    add_version_info.update_version()

# Finish
print("UPDATE COMPLETE: Check update_log for failed transcripts and correct")

# <LICENSE>
# Copyright (C) 2016-2021 VariantValidator Contributors
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
