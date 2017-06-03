"""
    Routines for runtime docstring argument injection

    This file contains HEAVILY modified routines from sphinx.ext.napoleon, from version 1.4.4

    This has been vendored into MDT because the modification makes use of
    private functions which have already changed in the dev branch.

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    :copyright: Copyright 2007-2016 by the Sphinx team, see sphinxlicense/AUTHORS.
    :license: BSD, see sphinxlicense/LICENSE for details.
"""

from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library
standard_library.install_aliases()

from past.builtins import basestring
import collections
import re

import sys

_google_section_regex = re.compile(r'^(\s|\w)+:\s*$')
_google_typed_arg_regex = re.compile(r'\s*(.+?)\s*\(\s*(.+?)\s*\)')
_single_colon_regex = re.compile(r'(?<!:):(?!:)')
_xref_regex = re.compile(r'(:\w+:\S+:`.+?`|:\S+:`.+?`|`.+?`)')
_bullet_list_regex = re.compile(r'^(\*|\+|\-)(\s+\S|\s*$)')
_enumerated_list_regex = re.compile(
    r'^(?P<paren>\()?'
    r'(\d+|#|[ivxlcdm]+|[IVXLCDM]+|[a-zA-Z])'
    r'(?(paren)\)|\.)(\s+\S|\s*$)')


class GoogleDocArgumentInjector(object):

    SECTIONS = set('args arguments parameters'.split())

    def __init__(self, docstring, prepare=True):
        # this routine has been modified - it's been streamlined for the current purpose

        if prepare:
            if docstring is None:
                self.docstring = []
            else:
                self.docstring = prepare_docstring(docstring)

        elif isinstance(docstring, basestring):
            self.docstring = docstring.splitlines()

        else:
            self.docstring = docstring

        self.lines_before_args = []
        self.arg_section = []
        self.lines_after_args = []
        self.args = collections.OrderedDict()
        self.arg_indent = None
        self.arg_section_name = 'Args'  # default, can be overwritten by the actual section name


        self._what = 'function'
        self._lines = list(self.docstring)
        self._line_iter = modify_iter(self.docstring, modifier=lambda s: s.rstrip())
        self._parsed_lines = []
        self._is_in_section = False
        self._section_indent = 0

        self._sections = {
                'args': self._parse_parameters_section,
                'arguments': self._parse_parameters_section,
                'attributes': None,
                'example': None,
                'examples': None,
                'keyword args': None,
                'keyword arguments': None,
                'methods': None,
                'note': None,
                'notes': None,
                'other parameters': None,
                'parameters': self._parse_parameters_section,
                'return': None,
                'returns': None,
                'raises': None,
                'references': None,
                'see also': None,
                'todo': None,
                'warning': None,
                'warnings': None,
                'warns': None,
                'yield': None,
                'yields': None,
            }

        self.parse()

    def new_docstring(self):
        """ Create a new docstring with the current state of the argument list

        Returns:
            str: docstring with modified argument list
        """
        newlines = list(self.lines_before_args)
        if self.args:
            newlines.append(' '*self.arg_indent + self.arg_section_name + ':')
            newlines.extend(self._indent(list(self.args.values()), self.arg_indent+4))
            newlines.append('')
        newlines.extend(self.lines_after_args)

        return '\n'.join(newlines)

    def parse(self):
        """ This method is a modified version of GoogleDocstring._parse
        """
        self._parsed_lines = self._consume_empty()

        found_args = lines_are_args = False

        while self._line_iter.has_next():
            if self._is_section_header():
                try:
                    section = self._consume_section_header()
                    self._is_in_section = True
                    self._section_indent = self._get_current_indent()

                    lines = [section + ':']

                    if section.lower() in self.SECTIONS:
                        lines.extend(self._sections[section.lower()](section))
                        found_args = True
                        lines_are_args = True
                    else:
                        lines.extend(self._consume_to_next_section())

                finally:
                    self._is_in_section = False
                    self._section_indent = 0
            else:
                if not self._parsed_lines:
                    lines = self._consume_contiguous()+self._consume_empty()
                else:
                    lines = self._consume_to_next_section()

            if lines_are_args:
                lines_are_args = False
                self.arg_section.extend(lines)
            elif found_args:
                self.lines_after_args.extend(lines)
            else:
                self.lines_before_args.extend(lines)

            self._parsed_lines.extend(lines)
            if self.arg_indent is None:
                self.arg_indent = self._get_current_indent()


    def _parse_parameters_section(self, section):
        """ This method was heavily modified to store information instead of formatting it for rst
        """
        self.arg_section_name = section
        fields = self._consume_fields()
        num_indent = self._get_current_indent()
        self.arg_indent = num_indent
        lines = []

        for _name, _type, _desc in fields:
            _desc = self._strip_empty(_desc)
            if isinstance(_desc, list):
                _desc = '\n    '.join(_desc)

            if _type:
                line = '%s (%s): %s' % (_name, _type, _desc)
            else:
                line = '%s: %s' % (_name, _desc)
            self.args[_name.lstrip('\*')] = line

            lines.append(line)

        lines = self._indent(lines, num_indent+4)
        if lines[-1].strip():
            lines.append('')
        return lines

    def _indent(self, lines, n=4):
        # MDT: modified to include breaks within lines
        sp = ' ' * n
        return [sp + line.replace('\n', '\n'+sp) for line in lines]


    ######################################################
    ### All routines below are unmodified              ###
    ######################################################
    def lines(self):
        """Return the parsed lines of the docstring in reStructuredText format.

        Returns
        -------
        :obj:`list` of :obj:`str`
            The lines of the docstring in a list.

        """
        return self._parsed_lines

    def _consume_indented_block(self, indent=1):
        lines = []
        line = self._line_iter.peek()
        while(not self._is_section_break() and
              (not line or self._is_indented(line, indent))):
            lines.append(next(self._line_iter))
            line = self._line_iter.peek()
        return lines

    def _consume_contiguous(self):
        lines = []
        while (self._line_iter.has_next() and
               self._line_iter.peek() and
               not self._is_section_header()):
            lines.append(next(self._line_iter))
        return lines

    def _consume_empty(self):
        lines = []
        line = self._line_iter.peek()
        while self._line_iter.has_next() and not line:
            lines.append(next(self._line_iter))
            line = self._line_iter.peek()
        return lines

    def _consume_field(self, parse_type=True, prefer_type=False):
        line = next(self._line_iter)

        before, colon, after = self._partition_field_on_colon(line)
        _name, _type, _desc = before, '', after

        if parse_type:
            match = _google_typed_arg_regex.match(before)
            if match:
                _name = match.group(1)
                _type = match.group(2)

        _name = self._escape_args_and_kwargs(_name)

        if prefer_type and not _type:
            _type, _name = _name, _type
        indent = self._get_indent(line) + 1
        _desc = [_desc] + self._dedent(self._consume_indented_block(indent))
        _desc = self.__class__(_desc, prepare=False).lines()
        return _name, _type, _desc

    def _consume_fields(self, parse_type=True, prefer_type=False):
        self._consume_empty()
        fields = []
        while not self._is_section_break():
            _name, _type, _desc = self._consume_field(parse_type, prefer_type)
            if _name or _type or _desc:
                fields.append((_name, _type, _desc,))
        return fields

    def _consume_section_header(self):
        section = next(self._line_iter)
        stripped_section = section.strip(':')
        if stripped_section.lower() in self._sections:
            section = stripped_section
        return section

    def _consume_to_end(self):
        lines = []
        while self._line_iter.has_next():
            lines.append(next(self._line_iter))
        return lines

    def _consume_to_next_section(self):
        self._consume_empty()
        lines = []
        while not self._is_section_break():
            lines.append(next(self._line_iter))
        return lines + self._consume_empty()

    def _dedent(self, lines, full=False):
        if full:
            return [line.lstrip() for line in lines]
        else:
            min_indent = self._get_min_indent(lines)
            return [line[min_indent:] for line in lines]

    def _escape_args_and_kwargs(self, name):
        if name[:2] == '**':
            return r'\*\*' + name[2:]
        elif name[:1] == '*':
            return r'\*' + name[1:]
        else:
            return name

    def _fix_field_desc(self, desc):
        if self._is_list(desc):
            desc = [''] + desc
        elif desc[0].endswith('::'):
            desc_block = desc[1:]
            indent = self._get_indent(desc[0])
            block_indent = self._get_initial_indent(desc_block)
            if block_indent > indent:
                desc = [''] + desc
            else:
                desc = ['', desc[0]] + self._indent(desc_block, 4)
        return desc

    def _get_current_indent(self, peek_ahead=0):
        line = self._line_iter.peek(peek_ahead + 1)[peek_ahead]
        while line != self._line_iter.sentinel:
            if line:
                return self._get_indent(line)
            peek_ahead += 1
            line = self._line_iter.peek(peek_ahead + 1)[peek_ahead]
        return 0

    def _get_indent(self, line):
        for i, s in enumerate(line):
            if not s.isspace():
                return i
        return len(line)

    def _get_initial_indent(self, lines):
        for line in lines:
            if line:
                return self._get_indent(line)
        return 0

    def _get_min_indent(self, lines):
        min_indent = None
        for line in lines:
            if line:
                indent = self._get_indent(line)
                if min_indent is None:
                    min_indent = indent
                elif indent < min_indent:
                    min_indent = indent
        return min_indent or 0

    def _is_indented(self, line, indent=1):
        for i, s in enumerate(line):
            if i >= indent:
                return True
            elif not s.isspace():
                return False
        return False

    def _is_list(self, lines):
        if not lines:
            return False
        if _bullet_list_regex.match(lines[0]):
            return True
        if _enumerated_list_regex.match(lines[0]):
            return True
        if len(lines) < 2 or lines[0].endswith('::'):
            return False
        indent = self._get_indent(lines[0])
        next_indent = indent
        for line in lines[1:]:
            if line:
                next_indent = self._get_indent(line)
                break
        return next_indent > indent

    def _is_section_header(self):
        section = self._line_iter.peek().lower()
        match = _google_section_regex.match(section)
        if match and section.strip(':') in self._sections:
            header_indent = self._get_indent(section)
            section_indent = self._get_current_indent(peek_ahead=1)
            return section_indent > header_indent

        return False

    def _is_section_break(self):
        line = self._line_iter.peek()
        return (not self._line_iter.has_next() or
                self._is_section_header() or
                (self._is_in_section and
                    line and
                    not self._is_indented(line, self._section_indent)))

    def _partition_field_on_colon(self, line):
        before_colon = []
        after_colon = []
        colon = ''
        found_colon = False
        for i, source in enumerate(_xref_regex.split(line)):
            if found_colon:
                after_colon.append(source)
            else:
                m = _single_colon_regex.search(source)
                if (i % 2) == 0 and m:
                    found_colon = True
                    colon = source[m.start(): m.end()]
                    before_colon.append(source[:m.start()])
                    after_colon.append(source[m.end():])
                else:
                    before_colon.append(source)

        return ("".join(before_colon).strip(),
                colon,
                "".join(after_colon).strip())

    def _strip_empty(self, lines):
        if lines:
            start = -1
            for i, line in enumerate(lines):
                if line:
                    start = i
                    break
            if start == -1:
                lines = []
            end = -1
            for i in reversed(range(len(lines))):
                line = lines[i]
                if line:
                    end = i
                    break
            if start > 0 or end + 1 < len(lines):
                lines = lines[start:end + 1]
        return lines



class peek_iter(object):
    """An iterator object that supports peeking ahead.
    Parameters
    ----------
    o : iterable or callable
        `o` is interpreted very differently depending on the presence of
        `sentinel`.
        If `sentinel` is not given, then `o` must be a collection object
        which supports either the iteration protocol or the sequence protocol.
        If `sentinel` is given, then `o` must be a callable object.
    sentinel : any value, optional
        If given, the iterator will call `o` with no arguments for each
        call to its `next` method; if the value returned is equal to
        `sentinel`, :exc:`StopIteration` will be raised, otherwise the
        value will be returned.
    See Also
    --------
    `peek_iter` can operate as a drop in replacement for the built-in
    `iter <https://docs.python.org/2/library/functions.html#iter>`_ function.
    Attributes
    ----------
    sentinel
        The value used to indicate the iterator is exhausted. If `sentinel`
        was not given when the `peek_iter` was instantiated, then it will
        be set to a new object instance: ``object()``.
    """
    def __init__(self, *args):
        """__init__(o, sentinel=None)"""
        self._iterable = iter(*args)
        self._cache = collections.deque()
        if len(args) == 2:
            self.sentinel = args[1]
        else:
            self.sentinel = object()

    def __iter__(self):
        return self

    def __next__(self, n=None):
        # note: prevent 2to3 to transform self.next() in next(self) which
        # causes an infinite loop !
        return getattr(self, 'next')(n)

    def _fillcache(self, n):
        """Cache `n` items. If `n` is 0 or None, then 1 item is cached."""
        if not n:
            n = 1
        try:
            while len(self._cache) < n:
                self._cache.append(next(self._iterable))
        except StopIteration:
            while len(self._cache) < n:
                self._cache.append(self.sentinel)

    def has_next(self):
        """Determine if iterator is exhausted.
        Returns
        -------
        bool
            True if iterator has more items, False otherwise.
        Note
        ----
        Will never raise :exc:`StopIteration`.
        """
        return self.peek() != self.sentinel

    def next(self, n=None):
        """Get the next item or `n` items of the iterator.
        Parameters
        ----------
        n : int or None
            The number of items to retrieve. Defaults to None.
        Returns
        -------
        item or list of items
            The next item or `n` items of the iterator. If `n` is None, the
            item itself is returned. If `n` is an int, the items will be
            returned in a list. If `n` is 0, an empty list is returned.
        Raises
        ------
        StopIteration
            Raised if the iterator is exhausted, even if `n` is 0.
        """
        self._fillcache(n)
        if not n:
            if self._cache[0] == self.sentinel:
                raise StopIteration
            if n is None:
                result = self._cache.popleft()
            else:
                result = []
        else:
            if self._cache[n - 1] == self.sentinel:
                raise StopIteration
            result = [self._cache.popleft() for i in range(n)]
        return result

    def peek(self, n=None):
        """Preview the next item or `n` items of the iterator.
        The iterator is not advanced when peek is called.
        Returns
        -------
        item or list of items
            The next item or `n` items of the iterator. If `n` is None, the
            item itself is returned. If `n` is an int, the items will be
            returned in a list. If `n` is 0, an empty list is returned.
            If the iterator is exhausted, `peek_iter.sentinel` is returned,
            or placed as the last item in the returned list.
        Note
        ----
        Will never raise :exc:`StopIteration`.
        """
        self._fillcache(n)
        if n is None:
            result = self._cache[0]
        else:
            result = [self._cache[i] for i in range(n)]
        return result


class modify_iter(peek_iter):
    """An iterator object that supports modifying items as they are returned.
    Parameters
    ----------
    o : iterable or callable
        `o` is interpreted very differently depending on the presence of
        `sentinel`.
        If `sentinel` is not given, then `o` must be a collection object
        which supports either the iteration protocol or the sequence protocol.
        If `sentinel` is given, then `o` must be a callable object.
    sentinel : any value, optional
        If given, the iterator will call `o` with no arguments for each
        call to its `next` method; if the value returned is equal to
        `sentinel`, :exc:`StopIteration` will be raised, otherwise the
        value will be returned.
    modifier : callable, optional
        The function that will be used to modify each item returned by the
        iterator. `modifier` should take a single argument and return a
        single value. Defaults to ``lambda x: x``.
        If `sentinel` is not given, `modifier` must be passed as a keyword
        argument.
    Attributes
    ----------
    modifier : callable
        `modifier` is called with each item in `o` as it is iterated. The
        return value of `modifier` is returned in lieu of the item.
        Values returned by `peek` as well as `next` are affected by
        `modifier`. However, `modify_iter.sentinel` is never passed through
        `modifier`; it will always be returned from `peek` unmodified.
    Example
    -------
    >>> a = ["     A list    ",
    ...      "   of strings  ",
    ...      "      with     ",
    ...      "      extra    ",
    ...      "   whitespace. "]
    >>> modifier = lambda s: s.strip().replace('with', 'without')
    >>> for s in modify_iter(a, modifier=modifier):
    ...   print('"%s"' % s)
    "A list"
    "of strings"
    "without"
    "extra"
    "whitespace."
    """
    def __init__(self, *args, **kwargs):
        """__init__(o, sentinel=None, modifier=lambda x: x)"""
        if 'modifier' in kwargs:
            self.modifier = kwargs['modifier']
        elif len(args) > 2:
            self.modifier = args[2]
            args = args[:2]
        else:
            self.modifier = lambda x: x
        if not callable(self.modifier):
            raise TypeError('modify_iter(o, modifier): '
                            'modifier must be callable')
        super().__init__(*args)

    def _fillcache(self, n):
        """Cache `n` modified items. If `n` is 0 or None, 1 item is cached.
        Each item returned by the iterator is passed through the
        `modify_iter.modified` function before being cached.
        """
        if not n:
            n = 1
        try:
            while len(self._cache) < n:
                self._cache.append(self.modifier(next(self._iterable)))
        except StopIteration:
            while len(self._cache) < n:
                self._cache.append(self.sentinel)


def prepare_docstring(s, ignore=1):
    """Convert a docstring into lines of parseable reST.  Remove common leading
    indentation, where the indentation of a given number of lines (usually just
    one) is ignored.
    Return the docstring as a list of lines usable for inserting into a docutils
    ViewList (used as argument of nested_parse().)  An empty line is added to
    act as a separator between this docstring and following content.
    """
    lines = s.expandtabs().splitlines()
    # Find minimum indentation of any non-blank lines after ignored lines.
    margin = sys.maxsize
    for line in lines[ignore:]:
        content = len(line.lstrip())
        if content:
            indent = len(line) - content
            margin = min(margin, indent)
    # Remove indentation from ignored lines.
    for i in range(ignore):
        if i < len(lines):
            lines[i] = lines[i].lstrip()
    if margin < sys.maxsize:
        for i in range(ignore, len(lines)):
            lines[i] = lines[i][margin:]
    # Remove any leading blank lines.
    while lines and not lines[0]:
        lines.pop(0)
    # make sure there is an empty line at the end
    if lines and lines[-1]:
        lines.append('')
    return lines