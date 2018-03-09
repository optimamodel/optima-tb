# -*- coding: utf-8 -*-
"""
Optima TB can be imported in the following manner:

from optima_tb.project import Project
or
import optima_tb as optb
or
from optima_tb import *

his program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


@author  : sjarvis
@date    : 13jan2017
"""


optimalicense = '''
Optima TB: a TB optimization and analysis tool
Copyright (C) 2017 by the Optima Consortium
'''
print(optimalicense)


from ._version import __version__


## General modules
from uuid import uuid4 as uuid
from datetime import datetime; today = datetime.today
from copy import deepcopy as dcp

#from .tb import *

# NOTE - To configure logging in individual scripts, use the commands below to reset the
# logger settings
import logging.config
logging_conf = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'standard': {
            'format': '%(asctime)-20s %(levelname)-8s %(message)s',
            'datefmt': '%d-%m-%y %H:%M:%S'
        },
    },
    'handlers': {
        'default': {
            'level': 'DEBUG',
            'formatter': 'standard',
            'class': 'logging.StreamHandler',
        },
    },
    'loggers': {
        '': {
            'handlers': ['default'],
            'level': 'DEBUG',
        },
    }
}
logging.config.dictConfig(logging_conf)
import logging
logger = logging.getLogger()
logger.setLevel('INFO')