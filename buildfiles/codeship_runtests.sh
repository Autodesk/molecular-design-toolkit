#!/usr/bin/env bash

py.test -n 4 --junit-xml=/opt/reports/junit.${TESTENV}.xml | tee /opt/reports/pytest.${TESTENV}.log
exitstat=${PIPESTATUS[0]}

statline="$(tail -n1 /opt/reports/pytest.${TESTENV}.log)"

echo 'Test status:'
echo ${statline}

python ../../buildfiles/send_test_status.py ${exitstat} "${statline}"

exit ${exitstat}
