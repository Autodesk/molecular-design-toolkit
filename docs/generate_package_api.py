#!/usr/bin/env python
""" Auto-generate rst files for the MDT API.

This is more or less our customized replacement for sphinx-apidoc. The goal here is to
organize the auto-generated API documentation in a more user-focused way - to clearly group
the classes and functions exposed by various parts of the MDT *runtime* API.
"""
import os
import types
import collections

SubMod = collections.namedtuple('MDTSubModule', 'modname name toplevel')

SUBMODULES = [SubMod('geom', 'Geometry', True),
              SubMod('widgets', 'Widgets', False),
              SubMod('tools', 'Tools', True),
              SubMod('models', 'Energy models', False),
              SubMod('integrators', 'Integrators', False),
              SubMod('data', 'Data', False)]

DOCPATH = '_mdt_api'

IGNORES = set(('toplevel exports __builtins__ __file__ __name__ __package__ __path__ __module__ '
               '__all__').split())


def main():
    if not os.path.exists(DOCPATH):
        os.mkdir(DOCPATH)

    for m in SUBMODULES:
        rstfile = os.path.join(DOCPATH, 'moldesign.%s.rst' % m.modname)
        with open(rstfile, 'w') as outfile:
            print >> outfile, document_module(m)
            print '   %s' % rstfile


def document_module(submod):
    """ Create a string with appropriate Autodoc directives for a given submodule

    Args:
        submod (SubMod): description of the submodule

    Returns:
        str: RST file content
    """
    module = getattr(moldesign, submod.modname)
    moddoc, mod_members, allnames = get_module_members(module, submod)

    docname = '%s API' % submod.name
    output = ['='*len(docname), docname, '='*len(docname), '']
    sections = collections.OrderedDict([('Classes', []),
                                        ('Functions', []),
                                        ('Data', []),
                                        ('Unexported API', [])
                                        ])

    if moddoc is not None:
        output.extend([moddoc, ''])

    for attrname, attr in mod_members:
        if submod.toplevel and attrname in allnames:
            path = 'moldesign'
        else:
            path = 'moldesign.%s' % submod.modname

        directive = get_autodoc_directive(attr, attrname, path)

        if submod.toplevel and attrname not in allnames:
            sections['Unexported API'].extend(directive)
        elif isinstance(attr, type):
            sections['Classes'].extend(directive)
        elif isinstance(attr, types.FunctionType):
            sections['Functions'].extend(directive)
        else:
            sections['Data'].extend(directive)

    for key, lines in sections.iteritems():
        if not lines: continue
        output.extend([key, '='*len(key), ''])
        output.extend(lines)

    return '\n'.join(output)


def get_module_members(module, submod):
    """ Make a list of the attributes to document for this module

    Returns:
        (str): docstring for the module (or None if it doesn't have one)
        List[(str, object)]: list of tuples of the form (attribute name, attribute reference)
        set(str): list of attribute names in the __all__ variables
    """
    attrnames = dir(module)
    mod_members = []
    for attrname in attrnames:
        if attrname in IGNORES: continue

        attr = getattr(module, attrname)
        if isinstance(attr, types.ModuleType): continue

        if attrname == '__doc__':
            moddoc = attr
            continue
        else:
            mod_members.append((attrname, attr))

    mod_members.sort()

    if hasattr(module, '__all__') and submod.toplevel:
        allnames = set(module.__all__)
    else:
        allnames = set(attrnames)

    return moddoc, mod_members, allnames


def get_autodoc_directive(attr, attrname, path):
    """ Create autodoc directive for a given object in the module

    Returns:
        List[str]: list of lines (to be joined via ``'\n'.join()``)
    """
    attrdocname = '``%s`` %s'%(attrname, _obj_description(attr))
    directive = [attrdocname,
                 '-'*len(attrdocname), '',
                 '.. %s:: %s.%s'%(_autorst(attr), path, attrname)]

    if isinstance(attr, type):
        directive.append('    :members:\n    :undoc-members:\n    :show-inheritance:')

    directive.append('')

    return directive


def _obj_description(o):
    if isinstance(o, type):
        return 'class'
    elif isinstance(o, types.FunctionType):
        return 'function'
    else:
        return ''


def _autorst(o):
    if isinstance(o, type):
        return 'autoclass'
    elif isinstance(o, types.FunctionType):
        return 'autofunction'
    else:
        return 'autoattribute'


if __name__ == '__main__':
    import moldesign
    main()