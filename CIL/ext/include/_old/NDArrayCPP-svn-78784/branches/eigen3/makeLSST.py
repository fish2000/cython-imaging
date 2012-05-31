#!/usr/bin/env python

# This is a conversion script to create a customized version of ndarray
# for the Large Synoptic Survey Telescope project, reflecting the fact
# that LSST started using ndarray before it was a stable project in its
# own right.  This script allows selected changes from the main ndarray
# distribution to be merged into the LSST version.  Someday LSST will
# use the main ndarray distribution and this file will go away.

import re
import glob
import os

substitutions = [
    (re.compile(r'NDARRAY_(\w+)_hpp_INCLUDED'), r'LSST_NDARRAY_\1_h_INCLUDED'),
    (re.compile(r'NDARRAY_ASSERT'), r'LSST_NDARRAY_ASSERT'),
    (re.compile(r'ndarray(.*)\.hpp'), r'lsst/ndarray\1.h'),
    (re.compile(r'ndarray(.*)\.cc'), r'lsst/ndarray\1.cc'),
    (re.compile(r'^namespace ndarray {'), r'namespace lsst { namespace ndarray {'),
    (re.compile(r'^} // namespace ndarray'), r'}} // namespace lsst::ndarray'),
    (re.compile(r'^}} // namespace ndarray::(\w+)'), r'}}} // namespace lsst::ndarray::\1'),
    (re.compile(r'([^:])ndarray::'), r'\1lsst::ndarray::'),
    ]

def convertFile(infile, outfile):
    outfile.write("""// -*- lsst-c++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
""")
    for line in infile:
        for regex, replacement in substitutions:
            line = regex.sub(replacement, line)
        outfile.write(line)

def walkTree(inroot, outroot):
    if not os.path.isdir(outroot):
        os.makedirs(outroot)
    for root, dirs, files in os.walk(inroot, topdown=True):
        s = root.split("/")[1:]
        path = os.path.join(*s) if len(s) > 0 else ""
        if ".svn" in dirs:
            dirs.remove(".svn")
        for dirname in dirs:
            outdir = os.path.join(outroot, path, dirname)
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
        for infile in files:
            if infile.endswith(".hpp"):
                outfile = infile[:-4] + ".h"
            elif infile.endswith(".hpp.m4"):
                outfile = infile[:-7] + ".h.m4"
            elif infile.endswith(".cc"):
                outfile = infile
            else:
                continue
            convertFile(open(os.path.join(inroot, path, infile), 'r'), 
                        open(os.path.join(outroot, path, outfile), 'w'))

def main(root):
    walkTree('include', os.path.join(root, "include", "lsst"))
    walkTree('tests', os.path.join(root, "tests"))
    
if __name__ == "__main__":
    import sys
    main(sys.argv[1])
