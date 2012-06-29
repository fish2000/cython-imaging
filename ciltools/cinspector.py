#!/usr/bin/env python
# encoding: utf-8
"""
ciltools/cinspector.py

Created by FI$H 2000 on 2012-06-28.
Copyright (c) 2012 Objects In Space And Time, LLC. All rights reserved.
"""
from __future__ import with_statement

class CInspectorException(Exception):
    pass

from os.path import abspath, dirname, join
#from pybindgen import FileCodeSink
from pybindgen.gccxmlparser import ModuleParser

class CInspector(object):
    
    package_root = dirname(dirname(abspath(__file__)))
    module_names = ('cimg', 'cimg_library')
    
    def __init__(self, cimg_path=None):
        if cimg_path is None:
            cimg_path = join(self.package_root, 'CIL', 'ext', 'include', 'CImg.h')
        self.cimg_path = cimg_path
        self.cimg_methods = {}
        self.modules = dict([(mdname, None) for mdname in self.module_names])
        self.parsers = dict([(str(md), ModuleParser('CImg', md)) for md in self.modules.keys()])
    
    def __getstate__(self):
        out_dict = self.__dict__.copy()
        del out_dict['modules']
        del out_dict['parsers']
        if 'cimg_class_typemap' in out_dict:
            del out_dict['cimg_class_typemap']
        return out_dict
    
    def __setstate__(self, in_dict):
        self.__dict__ = in_dict
        self.modules = dict([(mdname, None) for mdname in self.module_names])
        self.parsers = dict([(str(md), ModuleParser('CImg', md)) for md in self.modules.keys()])
    
    def get_module(self, mdname):
        if not mdname in self.modules.keys():
            raise CInspectorException('Unknown module name: %s' % mdname)
        if self.modules[mdname] is None:
            print "*** Parsing CImg.h for namespace: %s" % mdname
            print "*** This may take a moment..."
            print ""
            self.modules[mdname] = self.parsers[mdname].parse([self.cimg_path])
            print ""
        return self.modules[mdname]
    
    def get_cimg_types(self):
        if not hasattr(self, 'cimg_types'):
            cimg_classes = self.get_module('cimg').classes
            cimg_typeclasses = filter(
                lambda cls: cls.name.startswith('type'),
                cimg_classes)
            cimg_typenames = map(
                lambda tcls: tcls.template_parameters[0],
                cimg_typeclasses)
            self.cimg_types = set(cimg_typenames)
        return self.cimg_types
    
    def get_cimg_class_names(self):
        if not hasattr(self, 'cimg_class_names'):
            cimg_member_names = self.get_module('cimg_library').keys()
            cimg_member_class_names = filter(
                lambda k: k.startswith('cimg_library::CImg<'),
                cimg_member_names)
            self.cimg_class_names = set(cimg_member_class_names)
        return self.cimg_class_names
    
    def get_cimg_class_typemap(self):
        if not hasattr(self, 'cimg_class_typemap'):
            cimg_member_classes = self.get_module('cimg_library').classes
            self.cimg_class_typemap = dict(
                map(
                    lambda cls: (cls.template_parameters[0], cls),
                    filter(
                        lambda cls: cls.name == 'CImg',
                        cimg_member_classes)))
        return self.cimg_class_typemap
    
    def get_cimg_methods_for_type(self, typename='char'):
        cimg_types = self.get_cimg_types()
        if not typename in cimg_types:
            raise CInspectorException("<T> typename %s isn't valid.")
        if not typename in self.cimg_methods.keys():
            pass


def main():
    from pprint import pprint
    import cPickle as pickle
    
    cinspector = CInspector()
    
    print ""
    pprint(cinspector.get_cimg_types())
    
    print ""
    pprint(cinspector.get_cimg_class_names())
    
    print ""
    pprint(cinspector.get_cimg_class_typemap())
    
    print ""
    for T, cls in cinspector.get_cimg_class_typemap().items():
        print "\t %s:" % cls.full_name
        print "\t\t typename %s" % T
        print "\t\t\t %5s Constructors" % len(cls.constructors)
        print "\t\t\t %5s Methods (from cls.methods)" % len(cls.methods)
        print "\t\t\t %5s Methods (from cls.get_all_methods())" % len(list(cls.get_all_methods()))
    
    print ""
    barrel_pth = join(CInspector.package_root, 'cinspect.dump')
    print "Pickling CInspector data to %s" % barrel_pth
    with open(barrel_pth, 'wb') as barrel:
        pickle.dump(cinspector, barrel)
    
    print ""
    print "Done."
    
    
    

if __name__ == '__main__':
    main()