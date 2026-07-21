import os

# ----------------------------------------
# DEFAULT LOG FILE
# ----------------------------------------
LOG_FILE = os.path.join(
    os.getcwd(),
    ".vv_populate.log"
)

# ----------------------------------------
# LOG LEVELS
# ----------------------------------------
CONSOLE_LEVEL = "INFO"
FILE_LEVEL = "ERROR"

# ----------------------------------------
# LOGGING CONFIG
# ----------------------------------------
LOGGING_CONFIG = {
    'version': 1,

    'formatters': {
        'simple': {
            'format': (
                '%(asctime)s | %(levelname)-8s | %(name)s | '
                '%(filename)s:%(lineno)d | %(message)s'
            ),
            'datefmt': '%Y-%m-%d %H:%M:%S',
        },
        'detailed': {
            'format': (
                '%(asctime)s | %(levelname)-8s | %(name)s | '
                '%(filename)s:%(lineno)d | %(message)s'
            ),
            'datefmt': '%Y-%m-%d %H:%M:%S',
        }
    },

    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'level': CONSOLE_LEVEL,
            'formatter': 'simple'
        },
        'file': {
            'class': 'logging.FileHandler',
            'level': FILE_LEVEL,
            'filename': LOG_FILE,
            'mode': 'a',
            'formatter': 'detailed',
        },
    },

    'loggers': {
        'populate_vdb': {
            'level': 'DEBUG',
            'handlers': ['console', 'file'],
            'propagate': False,
        },
        'VariantValidator': {
            'level': 'WARNING',
            'handlers': ['console', 'file'],
            'propagate': False,
        },
    }
}


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


