# -*- coding: utf-8 -*-
"""
This file performs all necessary imports, so Optima can be used either as

from optima import Project, Parameters [etc.]
or
import optima as op
or
from optima import *


Now, the legal part:

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Version: 2016nov03 by cliffk
"""

optimalicense = '''
Optima: optimization and analysis tool
Copyright (C) 2016 by the Optima Consortium
'''
print(optimalicense)

## General modules
from uuid import uuid4 as uuid
from datetime import datetime; today = datetime.today
from copy import deepcopy as dcp

#from .tb import *
