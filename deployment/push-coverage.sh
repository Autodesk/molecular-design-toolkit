#!/usr/bin/env bash
# combines coverage from all test environments and pushes it to coveralls.io

set -e

mkdir -p /opt/reports/env-coverage
cp /opt/reports/env-coverage/coverage.* /opt/reports/

if $(coverage combine /opt/reports/coverage.*); then
  coveralls || echo "Failed to upload coverage to coveralls.io"
else
   echo "No coverage files found, skipping coveralls upload"
fi
