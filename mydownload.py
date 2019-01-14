#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Lytro Power Tools - cameratool script"""

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

__prog__ = 'cameratool'
__version__ = '1.1'

import os
import sys

try:
    import lpt
except ImportError:
    mod_dir = os.path.dirname(os.path.realpath(__file__))
    lpt_dir = os.path.abspath(os.path.join(mod_dir, '../..'))
    sys.path.insert(0, lpt_dir)
    import lpt

from lpt.camera import builder
import functools
import textwrap

from lpt.camera import camerabin
from lpt.camera import builder

bld = builder.Build()
cam = camerabin.CameraControls()

nImages = 2

def main():

    # bld.cam.download_images("C:\PROJEKT\lytro\data\imagesCamera\Raw", False, nImages)

    subparsers = bld.parser.add_subparsers()


    ''' DOWNLOAD RECENT IMAGES '''
    p_import_pictures = subparsers.add_parser("download-images",
                                              description=textwrap.dedent(''' Downloads images to the local file
                                              system '''),

                                              help=bld.dwl_string)

    p_import_pictures.add_argument("--nImages",
                                   type=int,
                                   default=bld.nImages,
                                   dest="nImages",
                                   metavar=" ",
                                   help="Choose this option to specify a destination path other than the default path: "
                                        "{}".format(bld.nImages))

    p_import_pictures.add_argument("--path",
                                   type=str,
                                   default=bld.currentDir,
                                   dest="path",
                                   metavar=" ",
                                   help="Choose this option to specify a destination path other than the default path: "
                                        "{}".format(bld.currentDir))

    p_import_pictures.add_argument("--all",
                                   action="store_true",
                                   default=False,
                                   dest="dwl_all",
                                   help="import all pictures as opposed to only the pictures captured in latest run.")

    p_import_pictures.set_defaults(func=bld.cam.download_images)
    args = bld.parser.parse_args()

    bld.cam.handle_if_adb_not_running()

    args.func(args)



if __name__ == "__main__":
    main()
