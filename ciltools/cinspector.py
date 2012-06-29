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

import re
from os.path import abspath, dirname, join
from collections import OrderedDict
from pybindgen.gccxmlparser import ModuleParser

class CInspector(object):
    
    package_root = dirname(dirname(abspath(__file__)))
    module_names = ('cimg', 'cimg_library')
    
    def __init__(self, cimg_path=None):
        if cimg_path is None:
            cimg_path = join(self.package_root, 'CIL', 'ext', 'include', 'CImg.h')
        self.cimg_path = cimg_path
        self.cimg_methods = {}
        self.cimg_constructors = {}
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
    
    def get_cimg_methods_for_type(self, typename='char', exclude_pattern=None):
        cimg_types = self.get_cimg_types()
        if not typename in cimg_types:
            raise CInspectorException("get_cimg_methods_for_type(): <T> typename %s isn't valid.")
        if not typename in self.cimg_methods.keys():
            self.cimg_methods[typename] = dict()
            cls = self.get_cimg_class_typemap()[typename]
            for method_name, method in cls.methods.items():
                self.cimg_methods[typename][method.wrapper_name] = list()
                for method_overload in method.wrappers:
                    self.cimg_methods[typename][method.wrapper_name].append(dict(
                        overload_index=method_overload.overload_index,
                        deprecated=method_overload.deprecated,
                        is_const=method_overload.is_const,
                        is_static=method_overload.is_static,
                        is_virtual=method_overload.is_virtual,
                        is_pure_virtual=method_overload.is_pure_virtual,
                        method_name=method_overload.method_name,
                        meth_flags=method_overload.meth_flags,
                        visibility=method_overload.visibility,
                        template_parameters=method_overload.template_parameters,
                        return_value=method_overload.return_value.value,
                        #return_value_class_name=method_overload.return_value.cpp_class.name),
                        #return_value_py_name=method_overload.return_value.py_name,
                        ))
        out = None
        if exclude_pattern is None:
            out = self.cimg_methods[typename]
        else:
            excluder = re.compile(exclude_pattern)
            out = OrderedDict(
                sorted(
                    filter(
                        lambda tup: excluder.match(tup[0]) is None,
                        self.cimg_methods[typename].items()),
                    key=lambda tup: tup[0]))
        return out
    
    def get_cimg_constructors_for_type(self, typename='char'):
        cimg_types = self.get_cimg_types()
        if not typename in cimg_types:
            raise CInspectorException("get_cimg_constructors_for_type(): <T> typename %s isn't valid.")
        if not typename in self.cimg_constructors.keys():
            self.cimg_constructors[typename] = list()
            cls = self.get_cimg_class_typemap()[typename]
            for constructor in cls.constructors:
                self.cimg_constructors[typename].append(dict(
                    overload_index=constructor.overload_index,
                    parameters=[dict(name=p.name, value=p.value, default=p.default_value, ctype=p.ctype) for p in constructor.parameters],
                    meth_flags=constructor.meth_flags,
                    visibility=constructor.visibility))
        return self.cimg_constructors[typename]

def load_pickled_cinspector():
    from ciltools.cinspector import CInspector
    import cPickle as pickle
    out = None
    barrel_pth = join(CInspector.package_root, 'cinspect.dump')
    print "Loading pickled CInspector data from %s" % barrel_pth
    with open(barrel_pth, 'rb') as barrel:
        out = pickle.load(barrel)
    return out

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
        ctors = cinspector.get_cimg_constructors_for_type(T)
        methods = cinspector.get_cimg_methods_for_type(T)
        print "\t %s:" % cls.full_name
        print "\t\t typename %s" % T
        print "\t\t\t %5s Constructors (from cls.constructors)" % len(cls.constructors)
        print "\t\t\t %5s Constructors (from get_cimg_constructors_for_type())" % len(ctors)
        print "\t\t\t %5s Methods (from cls.methods)" % len(cls.methods)
        print "\t\t\t %5s Methods (from cls.get_all_methods())" % len(list(cls.get_all_methods()))
        print "\t\t\t %5s Methods (from get_cimg_methods_for_type())" % len(methods.keys())
    
    print ""
    barrel_pth = join(CInspector.package_root, 'cinspect.dump')
    print "Pickling CInspector data to %s" % barrel_pth
    with open(barrel_pth, 'wb') as barrel:
        pickle.dump(cinspector, barrel)
    
    print ""
    print "Done."


if __name__ == '__main__':
    main()