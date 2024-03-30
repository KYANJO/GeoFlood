"""Additions for ForestClaw support"""

from __future__ import absolute_import
import os
import logging
import logging.config

import sys

sys.path.append('../../../scripts')

# Default logging configuration file
_DEFAULT_LOG_CONFIG_PATH = os.path.join(os.path.dirname(__file__),
                                        '../pyclaw/log.config')
del os

# Setup loggers
logging.config.fileConfig(_DEFAULT_LOG_CONFIG_PATH)

__all__ = []

# Module imports - Note the only difference here is the geometry module
__all__.extend(['Controller', 'Dimension', 'Patch', 'Domain', 'Solution',
                'State', 'CFL'])
from pyclaw.controller import Controller
from pyclaw.solution import Solution
from pyclaw.state import State
from pyclaw.cfl import CFL
from pyclaw.geometry import Dimension
from pyclaw.geometry import Domain
from .geometry import Patch
