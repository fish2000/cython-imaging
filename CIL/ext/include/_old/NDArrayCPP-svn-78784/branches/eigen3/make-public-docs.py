#!/usr/bin/env python
import os
import glob
import shutil
import sys

def main(source,output,dry_run=False):
    def run(cmd):
        if dry_run:
            print cmd
        else:
            os.system(cmd)
    old = set(glob.glob(os.path.join(output,"*")))
    types = {".html":"text/html",
             ".css":"text/css",
             ".gif":"image/gif",
             ".png":"image/png",
             }
    new = set()
    for fullname in source:
        path,filename = os.path.split(fullname)
        shutil.copy(fullname,output)
        new.add(os.path.join(output,filename))
    deleted = old - new
    added = new - old
    print "Deleting %i files." % len(deleted)
    for fullname in deleted:
        run("svn delete %s" % fullname)
    print "Adding %i files." % len(added)
    for fullname in added:
        base,ext = os.path.splitext(fullname)
        run("svn add %s" % fullname)
        run("svn propset 'svn:mime-type' '%s' %s" % (types[ext],fullname))

if __name__=="__main__":
    main(glob.glob("doc/html/*"),"../doc/",dry_run=(len(sys.argv)>1 and sys.argv[1]=="dry"))
    
