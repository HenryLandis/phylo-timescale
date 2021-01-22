#!/usr/bin/env python

"""
Load class objects and set logger to default INFO setting
"""

from .Simulator import Simulator
from .Analysis import Analysis
from .logger import set_loglevel

set_loglevel(loglevel="INFO")
