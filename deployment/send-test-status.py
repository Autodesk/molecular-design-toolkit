#!/usr/bin/env python
"""
This script accesses the github API to send a custom status message
about test results
"""
import os
import sys

import github

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('exitcode', type=str)
parser.add_argument('msg', type=str)
parser.add_argument('--deployed', action='store_true')
args = parser.parse_args()

status = {'0':'success', 'na':'pending'}.get(args.exitcode, 'failure')

missing_env = []
for key in 'CI_COMMIT_ID TESTENV GITHUB_REPO_TOKEN CI_PROJECT_ID CI_BUILD_ID PYVERSION'.split():
    if not os.environ.get(key, None):
        missing_env.append(key)

# set by codeship CI
sha = os.environ.get('CI_COMMIT_ID', '_no_commitid')
testenv = os.environ.get('TESTENV', '_notestenv')
ghtoken = os.environ.get('GITHUB_REPO_TOKEN', '_notoken')
# projid = os.environ.get('CI_PROJECT_ID', '_no_projid')
projid = '214515'  # hardcoded for now
buildid = os.environ.get('CI_BUILD_ID', '_no_buildid')
pyversion = os.environ.get('PYVERSION', '_no_pyversion')


data = dict(state=status,
            target_url='https://app.codeship.com/projects/%s/builds/%s' %
                       (projid, buildid),
            description=args.msg.replace("=","").strip(),
            context='%s/py%s' % (testenv, pyversion))


if missing_env:
    print("Not sending status update b/c of missing env vars: %s" % ','.join(missing_env))
    print(data)
    sys.exit(0)


g = github.Github(ghtoken)
repo = g.get_repo('autodesk/molecular-design-toolkit')
commit = repo.get_commit(sha)

if args.deployed:
    raise NotImplementedError
else:
    commit.create_status(**data)

