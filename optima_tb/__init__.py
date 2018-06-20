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

# Configure the logging system
import logging
logger = logging.getLogger() # Get the root logger, keep its level
logger.setLevel('INFO')
h = logging.StreamHandler()
h.setFormatter(logging.Formatter(fmt='%(asctime)-20s %(levelname)-8s %(message)s',datefmt='%d-%m-%y %H:%M:%S'))
logger.addHandler(h)
del h
logger.critical('Optima TB: a TB optimization and analysis tool\nCopyright (C) 2018 by the Optima Consortium')

# Set Numpy to throw errors for numerical warnings (these generally indicate something has gone wrong)
import numpy as np
np.seterr(all='raise')



