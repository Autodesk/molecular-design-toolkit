#!/usr/bin/env python2.7
# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Multiple inheritance for your dockerfiles.
Requires: python 2.7, docker-py, pyyaml (RUN: easy_install pip; pip install docker-py pyyaml)
"""
import sys
import os
import textwrap
from collections import OrderedDict
from io import StringIO, BytesIO
import argparse
import pprint

import docker, docker.utils
import yaml


class DockerMaker(object):
    def __init__(self, makefile, repository=None,
                 build_images=True,
                 print_dockerfiles=False,
                 no_cache=False,
                 tag=None,
                 pull=False):

        self._sources = set()
        self.makefile_path = makefile
        self.img_defs = self.parse_yaml(self.makefile_path)
        self.all_targets = self.img_defs.pop('_ALL_', None)

        # Connect to docker daemon if necessary
        if build_images:
            connection = docker.utils.kwargs_from_env()
            connection['tls'].assert_hostname = False
            self.client = docker.Client(**connection)
        else:
            self.client = None

        if repository and repository[-1] not in '/:':
            self.repo = repository + '/'
        elif repository is None:
            self.repo = ''
        else:
            self.repo = repository
        self.tag = tag
        self.build_images = build_images
        self.print_dockerfiles = print_dockerfiles
        self.pull = pull
        self.no_cache = no_cache

    def parse_yaml(self, filename):
        fname = os.path.expanduser(filename)
        print 'READING %s' % os.path.expanduser(fname)
        if fname in self._sources: raise ValueError('Circular _SOURCE_')
        self._sources.add(fname)

        with open(fname, 'r') as yaml_file:
            yamldefs = yaml.load(yaml_file)

        sourcedefs = {}
        for s in yamldefs.get('_SOURCES_', []):
            sourcedefs.update(self.parse_yaml(s))

        sourcedefs.update(yamldefs)
        return sourcedefs

    def build(self, image):
        """
        Drives the build of the final image - get the list of steps and execute them.
        :param image: name of the image from the yaml file to build
        :return: final tagged image name
        """
        print 'docker-make starting build for %s' % image
        build_steps = self.generate_build_order(image)
        for istep, step in enumerate(build_steps):
            print '  **** DockerMake Step %d/%d: %s ***' % (istep + 1, len(build_steps), ','.join(step.images))
            print '     * Build directory: %s' % step.build_dir
            print '     * Target image name: %s' % step.tag
            dockerfile = '\n'.join(step.dockerfile)

            # build the image
            if self.build_images:
                self.build_step(step, dockerfile)

            # Dump the dockerfile to a file
            if self.print_dockerfiles:
                if not os.path.exists('docker_makefiles'):
                    os.makedirs('docker_makefiles')
                if '/' in step.tag:
                    filename = 'docker_makefiles/Dockerfile.%s' % image
                else:
                    filename = 'docker_makefiles/Dockerfile.%s' % step.tag
                with open(filename, 'w') as dfout:
                    print >> dfout, dockerfile

        return step.tag

    def build_step(self, step, dockerfile):
        """
        Drives an individual build step. Build steps are separated by build_directory.
        If a build has zero one or less build_directories, it will be built in a single
        step.
        """
        # set up the build context
        build_args = dict(decode=True, tag=step.tag, pull=self.pull,
                          fileobj=None, path=None, dockerfile=None,
                          nocache=self.no_cache)
        if step.build_dir is not None:
            tempname = '_docker_make_tmp/'
            tempdir = '%s/%s' % (step.build_dir, tempname)
            temp_df = tempdir + 'Dockerfile'
            if not os.path.isdir(tempdir):
                os.makedirs(tempdir)
            with open(temp_df, 'w') as df_out:
                print >> df_out, dockerfile

            build_args['path'] = os.path.abspath(step.build_dir)
            build_args['dockerfile'] = tempname + 'Dockerfile'
        else:
            build_args['fileobj'] = StringIO(unicode(dockerfile))

        # start the build
        stream = self.client.build(**build_args)

        # monitor the output
        for item in stream:
            if item.keys() == ['stream']:
                print item['stream'].strip()
            elif 'errorDetail' in item or 'error' in item:
                raise BuildError(dockerfile, item, build_args)
            else:
                print item

        # remove the temporary dockerfile
        if step.build_dir is not None:
            os.unlink(temp_df)
            os.rmdir(tempdir)

    def generate_build_order(self, image):
        """
        Separate the build into a series of one or more intermediate steps.
        Each specified build directory gets its own step
        """
        repo_name = self.repo + image
        if self.tag:
            if ':' in repo_name:
                repo_name += '-' + self.tag
            else:
                repo_name += ':' + self.tag
        dependencies = self.sort_dependencies(image)
        base = self.get_external_base_image(image, dependencies)

        build_steps = [BuildStep(base)]
        step = build_steps[0]
        for d in dependencies:
            dep_definition = self.img_defs[d]
            mydir = dep_definition.get('build_directory', None)
            if mydir is not None:
                mydir = os.path.expanduser(mydir)  # expands `~` to home directory
                if step.build_dir is not None:
                    # Create a new build step if there's already a build directory
                    step.tag = '%dbuild_%s' % (len(build_steps), image)
                    build_steps.append(BuildStep(step.tag))
                    step = build_steps[-1]
                step.build_dir = mydir

            step.images.append(d)
            if 'build' in dep_definition:
                step.dockerfile.append('\n#Commands for %s' % d)
                step.dockerfile.append(dep_definition['build'])
            else:
                step.dockerfile.append('\n####end of requirements for %s\n' % d)

        # Sets the last step's name to the final build target
        step.tag = repo_name
        for step in build_steps:
            step.dockerfile.insert(0, '#Build directory: %s\n#tag: %s' %
                                   (step.build_dir, step.tag))
        return build_steps

    def sort_dependencies(self, com, dependencies=None):
        """
        Topologically sort the docker commands by their requirements
        TODO: sort using a "maximum common tree"?
        :param com: process this docker image's dependencies
        :param dependencies: running cache of sorted dependencies (ordered dict)
        :return type: OrderedDict
        """
        if dependencies is None: dependencies = OrderedDict()

        if com in dependencies: return
        requires = self.img_defs[com].get('requires', [])
        assert type(requires) == list, 'Requirements for %s are not a list' % com

        for dep in requires:
            self.sort_dependencies(dep, dependencies)
        if com in dependencies:
            raise ValueError('Circular dependency found', dependencies)
        dependencies[com] = None
        return dependencies

    def get_external_base_image(self, image, dependencies):
        """
        Makes sure that this image has exactly one external base image
        """
        base = None
        base_for = None
        for d in dependencies:
            this_base = self.img_defs[d].get('FROM', None)
            if this_base is not None and base is not None and this_base != base:
                error = ('Multiple external dependencies: image %s depends on:\n' % image +
                         '  %s (FROM: %s), and\n' % (base_for, base) +
                         '  %s (FROM: %s).' % (d, this_base))
                raise ValueError(error)
            if this_base is not None:
                base = this_base
                base_for = d
        if base is None:
            raise ValueError("No base image found in %s's dependencies" % image)
        return base


class BuildError(Exception):
    def __init__(self, dockerfile, item, build_args):
        with open('dockerfile.fail', 'w') as dff:
            print>> dff, dockerfile
        with BytesIO() as stream:
            print >> stream, '\n   -------- Docker daemon output --------'
            pprint.pprint(item, stream, indent=4)
            print >> stream, '   -------- Arguments to client.build --------'
            pprint.pprint(build_args, stream, indent=4)
            print >> stream, 'This dockerfile was written to dockerfile.fail'
            stream.seek(0)
            super(BuildError, self).__init__(stream.read())


class BuildStep(object):
    def __init__(self, baseimage):
        self.dockerfile = ['FROM %s\n' % baseimage]
        self.tag = None
        self.build_dir = None
        self.images = []


def main():
    args = make_arg_parser().parse_args()

    # Help and exit
    if args.help_yaml:
        print_yaml_help()
        return

    # Otherwise, parse the yaml file
    maker = DockerMaker(args.makefile, repository=args.repository,
                        build_images=not (args.no_build or args.list),
                        print_dockerfiles=(args.print_dockerfiles or args.no_build),
                        pull=args.pull, no_cache=args.no_cache, tag=args.tag)

    if args.list:
        print 'TARGETS in `%s`' % args.makefile
        for item in maker.img_defs.keys(): print ' *', item
        return

    # Assemble custom requirements target
    if args.requires or args.name:
        assert args.requires and args.name
        assert args.name not in maker.img_defs
        maker.img_defs[args.name] = {'requires': args.requires}
        targets = [args.name]
    elif args.all:
        assert len(args.TARGETS) == 0, "Pass either a list of targets or `--all`, not both"
        if maker.all_targets is not None:
            targets = maker.all_targets
        else:
            targets = maker.img_defs.keys()
    else:
        targets = args.TARGETS

    if not targets:
        print 'No build targets specified!'
        print 'Targets in `%s`:' % args.makefile
        for item in maker.img_defs.keys(): print ' *', item
        return

    # Actually build the images! (or Dockerfiles)
    built, warnings = [], []
    for t in targets:
        name = maker.build(t)
        print '  docker-make built:', name
        built.append(name)
        if args.push_to_registry:
            success, w = push(maker, name)
            warnings.extend(w)
            if not success: built[-1] += ' -- PUSH FAILED'
            else: built[-1] += ' -- pushed to %s' % name.split('/')[0]

    # Summarize the build process
    print '\ndocker-make finished.'
    print 'Built: '
    for item in built: print ' *', item
    if warnings:
        print 'Warnings:'
        for item in warnings: print ' *', item


def push(maker, name):
    success = False
    warnings = []
    if '/' not in name or name.split('/')[0].find('.') < 0:
        warn = 'WARNING: could not push %s - ' \
               'repository name does not contain a registry URL' % name
        warnings.append(warn)
        print warn
    else:
        print '  Pushing %s to %s:' % (name, name.split('/')[0])
        line = {'error': 'no push information received'}
        _lastid = None
        for line in maker.client.push(name, stream=True):
            line = yaml.load(line)
            if 'status' in line:
                if line.get('id', None) == _lastid and line['status'] == 'Pushing':
                    print '\r', line['status'], line['id'], line.get('progress', ''),
                    sys.stdout.flush()
                else:
                    print line['status'], line.get('id', '')
                    _lastid = line.get('id', None)
            else:
                print line
        if 'error' in line:
            warnings.append('WARNING: push failed for %s. Message: %s' % (name, line['error']))
        else:
            success = True
    return success, warnings


def print_yaml_help():
    print "A brief introduction to writing Dockerfile.yml files:\n"

    print 'SYNTAX:'
    print printable_code("""[image_name]:
  build_directory: [relative path where the ADD and COPY commands will look for files]
  requires:
   - [other image name]
   - [yet another image name]
  FROM: [named_base_image]
  build: |
   RUN [something]
   ADD [something else]
   [Dockerfile commands go here]

[other image name]: ...
[yet another image name]: ...""")

    print
    print textwrap.fill("The idea is to write dockerfile commands for each specific "
                        'piece of functionality in the build field, and "inherit" all other'
                        ' functionality from a list of other components that your image requires. '
                        'If you need to add files with the ADD and COPY commands, specify the root'
                        ' directory for those files with build_directory. Your tree of '
                        '"requires" must have exactly one unique named base image '
                        'in the FROM field.')

    print '\n\nAN EXAMPLE:'
    print printable_code("""devbase:
 FROM: phusion/baseimage
 build: |
  RUN apt-get -y update && apt-get -y install build-essential

airline_data:
 requires:
  - devbase
 build_directory: sample_data/airline_data
 build: |
  ADD AirlinePassengers.csv

python_image:
 requires:
  - devbase
 build: |
  RUN apt-get -y update \
  && apt-get install -y python python-pip \
  && pip install pandas

data_science:
 requires:
  - python_image
  - airline_data""")


def printable_code(c):
    output = []
    dedented = textwrap.dedent(c)
    for line in dedented.split('\n'):
        output.append(' >> ' + line)
    return '\n'.join(output)


def make_arg_parser():
    parser = argparse.ArgumentParser(description=
                                     "NOTE: Docker environmental variables must be set.\n"
                                     "For a docker-machine, run "
                                     "`eval $(docker-machine env [machine-name])`")
    bo = parser.add_argument_group('Choosing what to build')
    bo.add_argument('TARGETS', nargs="*",
                    help='Docker images to build as specified in the YAML file')
    bo.add_argument('-f', '--makefile',
                    default='DockerMake.yml',
                    help='YAML file containing build instructions')
    bo.add_argument('-a', '--all', action='store_true',
                    help="Print or build all images (or those specified by _ALL_)")
    bo.add_argument('-l', '--list', action='store_true',
                    help='List all available targets in the file, then exit.')
    bo.add_argument('--requires', nargs="*",
                    help='Build a special image from these requirements. Requires --name')
    bo.add_argument('--name', type=str,
                    help="Name for custom docker images (requires --requires)")

    df = parser.add_argument_group('Dockerfiles')
    df.add_argument('-p', '--print_dockerfiles', action='store_true',
                    help="Print out the generated dockerfiles named `Dockerfile.[image]`")
    df.add_argument('-n', '--no_build', action='store_true',
                    help='Only print Dockerfiles, don\'t build them. Implies --print.')

    ca = parser.add_argument_group('Image caching')
    ca.add_argument('--pull', action='store_true',
                    help='Always try to pull updated FROM images')
    ca.add_argument('--no-cache', action='store_true',
                    help="Rebuild every layer")
    # TODO: add a way to invalidate a specific target

    rt = parser.add_argument_group('Repositories and tags')
    rt.add_argument('--repository', '-r', '-u',
                    help="Prepend this repository to all built images, e.g.\n"
                         "`docker-make hello-world -u quay.io/elvis` will tag the image "
                         "as `quay.io/elvis/hello-world`. You can add a ':' to the end to "
                         "image names into tags:\n `docker-make -u quay.io/elvis/repo: hello-world` "
                         "will create the image in the elvis repository: quay.io/elvis/repo:hello-world")
    rt.add_argument('--tag', '-t', type=str,
                    help='Tag all built images with this tag. If image names are ALREADY tags (i.e.,'
                         ' your repo name ends in a ":"), this will append the tag name with a dash. '
                         'For example: `docker-make hello-world -u elvis/repo: -t 1.0` will create '
                         'the image "elvis/repo:hello-world-1.0')
    rt.add_argument('--push-to-registry', '-P', action='store_true',
                    help='Push all built images to the repository specified '
                         '(only if image repository contains a URL) -- to push to dockerhub.com, '
                         'use index.docker.io as the registry)')

    hh = parser.add_argument_group('Help')
    hh.add_argument('--help-yaml', action='store_true',
                    help="Print summary of YAML file format and exit.")

    return parser


__license__ = """Copyright (c) 2016, Autodesk Research
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."""

if __name__ == '__main__': main()
