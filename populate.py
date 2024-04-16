import os
import re
import VariantValidator
from VariantValidator import update_vv_db
import add_version_info
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor

vval = VariantValidator.Validator()
ROOT = os.path.dirname(os.path.abspath(__file__))

def check_args():
    if len(sys.argv) != 3:
        print('Too few arguments. The command required is: python populate.py False all: Where True/False dictates whether to '
              'use the uta_transcripts_testing.txt file for testing or the default input file and All is one of all, ensembl or refseq')
        exit()

    testing = sys.argv[1]
    transcript_set = sys.argv[2]
    print(testing)

    if str(testing) in "False":
        infile = os.path.join(ROOT, 'uta_transcripts.txt')
    elif str(testing) in "True":
        infile = os.path.join(ROOT, 'uta_transcripts_testing.txt')
    else:
        print('Argument 1 must be True or False where use the uta_transcripts_testing.txt file is True or False dictating '
              'whether to use the uta_transcripts_testing.txt file for testing or the default input file')
        exit()

    if str(transcript_set) not in ["all", "ensembl", "refseq"]:
        print('Argument 2 must be all, ensembl or refseq to select a specific transcript set or not')
        exit()

    return testing, transcript_set, infile

def process_transcript(tx_id):
    if tx_id == 'ac':
        return
    if not re.search('.', tx_id):
        return
    if '..' in tx_id:
        return
    if re.search('^U', tx_id):
        return
    if 'NG' in tx_id:
        return

    print('Updating - ' + tx_id)

    if 'ENST' not in tx_id:
        if transcript_set == "ensembl":
            print("Bypass")
            return
        try:
            vval.update_transcript_record(tx_id)
            print('True')
        except BaseException as e:
            print('False')
            fo.write(tx_id + '\t' + str(e) + '\n')
            fo.flush()
            os.fsync(fo.fileno())
            print(tx_id + '\t' + str(e) + '\n')
            traceback.print_exc()
    else:
        if transcript_set == "refseq":
            print("Bypass")
            return
        try:
            vval.update_transcript_record(tx_id, genome_build='GRCh38', test=True)
            print('True')
        except VariantValidator.modules.utils.DatabaseConnectionError:
            try:
                vval.update_transcript_record(tx_id, genome_build='GRCh37', test=True)
                print('True')
            except VariantValidator.modules.utils.DatabaseConnectionError as e:
                print('False')
                fo.write(tx_id + '\t' + str(e) + '\n')
                fo.flush()
                os.fsync(fo.fileno())
                print(tx_id + '\t' + str(e) + '\n')
                traceback.print_exc()
            except BaseException as e:
                print('False')
                fo.write(tx_id + '\t' + str(e) + '\n')
                fo.flush()
                os.fsync(fo.fileno())
                print(tx_id + '\t' + str(e) + '\n')
                traceback.print_exc()
        except BaseException as e:
            print('False')
            fo.write(tx_id + '\t' + str(e) + '\n')
            fo.flush()
            os.fsync(fo.fileno())
            print(tx_id + '\t' + str(e) + '\n')
            traceback.print_exc()

testing, transcript_set, infile = check_args()

logfile = os.path.join(ROOT, 'update_log.txt')
fo = open(logfile, "w")

with open(infile) as tx_data, ProcessPoolExecutor(max_workers=4) as executor:
    tx_ids = (tx_line.split('\t')[0] for tx_line in tx_data)
    executor.map(process_transcript, tx_ids)

fo.close()

if str(testing) in "False":
    update_vv_db.update()
    add_version_info.update_version()

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