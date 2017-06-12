#!/usr/bin/env bash

# This script has only one purpose - to stop ci from running twice in certain cases. Those cases
# are (currently):
#    1) The --semvertag flag is present in the commit message, but this is not a "tagged" commit


# if --semvertag is present, only run CI if branch is a valid semantic version
if [[ "${CI_COMMIT_MESSAGE}" == *"--semvertag" ]];  then
  if [[ "${CI_BRANCH}" =~ ^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(-(0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(\.(0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*)?(\+[0-9a-zA-Z-]+(\.[0-9a-zA-Z-]+)*)?$ ]]; then
    echo "OK to run CI - '${CI_BRANCH}' is valid semantic version"
    exit 0
  else
    echo "STOPPING CI - '${CI_BRANCH}' is not a valid semantic version"
    exit 1
  fi
fi


echo "OK to run CI!"
