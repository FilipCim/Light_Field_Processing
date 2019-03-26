#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Lytro Power Tools - cameracontrols script"""

# <copyright>
# Copyright (c) 2011-2015 Lytro, Inc. All rights reserved.
# This software is the confidential and proprietary information of Lytro, Inc.
# You shall not disclose such confidential information and shall use it only in
# accordance with the license granted to you by Lytro, Inc.

# EXCEPT AS EXPRESSLY SET FORTH IN A WRITTEN LICENSE AGREEMENT WITH LICENSEE,
# LYTRO, INC. MAKES NO REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF
# THE SOFTWARE, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR
# NON-INFRINGEMENT. LYTRO, INC. SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY
# LICENSEE AS A RESULT OF USING, COPYING, MODIFYING OR DISTRIBUTING THIS
# SOFTWARE OR ITS DERIVATIVES.
# </copyright>

from __future__ import division

__prog__ = 'cameracontrols'
__version__ = '1.1'

import argparse
import textwrap
import os
import sys
import functools

try:
    import lpt
except ImportError:
    mod_dir = os.path.dirname(os.path.realpath(__file__))
    lpt_dir = os.path.abspath(os.path.join(mod_dir, '../..'))
    sys.path.insert(0, lpt_dir)
    import lpt

from lpt.camera import camerabin
from lpt.camera import builder

bld = builder.Build()
cam = camerabin.CameraControls()

def captureImage(count, sleep):
    bld.check_space_and_framing(count, cam.get_pictures_remaining(True))
    cam.disable_focus_bracketing(True)
    cam.disable_exp_bracketing(True)
    if count > 10:
        sleep_timer = cam.get_timeout_value(True)
        cam.disable_timeout(True)
        cam.capture_wait(count, sleep, True)
        cam.set_timeout(sleep_timer, True)
    else:
        cam.capture_wait(count, sleep, True)



def main():
    ''' CAPTURE IMAGE '''
    count = 1
    sleep = 1

    print('Capturing!')
    captureImage(count, sleep)

    # cam.handle_if_adb_not_running()
    #
    # cam.verbose = True
    # cam.init(False)
    # cam.disable_virtual_cable(True)



if __name__ == "__main__":
    main()
