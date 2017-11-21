#!/usr/bin/env bash
set -e

if [ -z "${CI_BRANCH}" ]; then
  echo "Error: env var CI_BRANCH not set"
  exit 1
fi


pyversion=$(python -m moldesign version | tail -n 1)

echo "Expecting moldesign==${pyversion}"
echo "Found moldesign==${pyversion}"


if [ "${pyversion}" == "${CI_BRANCH}" ]; then
  echo "All good"
  exit 0
else
  echo "Error: moldesign package version '${pyversion}' differs from its Git tag '${CI_BRANCH}'"
  exit 1
fi
