import logging.config
import settings


def configure_logging(
    console_level=None,
    file_level=None,
):
    #
    # STEP 1: Configure RotatingFileHandler
    #
    if file_level is None:
        f_level = settings.FILE_LEVEL
    else:
        f_level = file_level.upper()

    settings.LOGGING_CONFIG["handlers"]["file"] = {
        "class": "logging.handlers.RotatingFileHandler",
        "level": f_level,
        "filename": settings.LOG_FILE,
        "mode": "a",
        "maxBytes": 500000,
        "backupCount": 10,
        "formatter": "detailed",
    }

    #
    # STEP 2: Configure console handler
    #
    if console_level is None:
        c_level = settings.CONSOLE_LEVEL
    else:
        c_level = console_level.upper()

    settings.LOGGING_CONFIG["handlers"]["console"]["level"] = c_level

    #
    # STEP 3: Apply logging configuration
    #
    logging.config.dictConfig(settings.LOGGING_CONFIG)


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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>