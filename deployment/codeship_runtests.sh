#!/usr/bin/env bash

# this script expects to run from the moldesign/_tests directory

VERSION="${TESTENV}.py${PYVERSION}"
PYTESTFLAGS="-n 4 --durations=20 --junit-xml=/opt/reports/junit.${VERSION}.xml --timeout=1800"
if [ "${TESTENV}" == "complete" ]; then
       PYTESTFLAGS="--cov .. --cov-config=./.coveragerc ${PYTESTFLAGS}"
fi


function send_status_update(){
     python ../../deployment/send_test_status.py "${1}" "${2}"
}


function check_if_tests_should_run(){
    echo "Should I run the tests in this environment?"

   if [ "${TESTENV}" == "complete" ]; then
       runthem=true
       echo "YES: always run in 'complete' environment"
   fi

   case "${CI_BRANCH}" in
       master|deploy|dev)
       runthem=true
       echo "YES: always run in branch \"${CI_BRANCH}\""
       ;;
   esac

   if [[ "${CI_COMMIT_MESSAGE}" == *"--testall"* ]];  then
       runthem=true
       echo "YES: found \"--testall\" flag in commit message"
   fi

   if [ "${runthem}" != "true" ]; then
     echo "SKIPPING tests in this environment."
     echo "To run these tests, add \"--testall\" to your commit message"
     echo "(or work in the dev or deploy branches)"
     send_status_update 0 "Skipped (dev/deploy branches only)"
     exit 0
   fi
}


function run_tests(){
    send_status_update "na" "Starting tests for ${VERSION}"

    py.test ${PYTESTFLAGS} | tee /opt/reports/pytest.${VERSION}.log
    exitstat=${PIPESTATUS[0]}

    statline="$(tail -n1 /opt/reports/pytest.${VERSION}.log)"

    echo 'Test status:'
    echo ${statline}

    send_status_update "${exitstat}" "${statline}"

    if [ "${TESTENV}" == "complete" ]; then
       if [[ ${exitstat} == 0 && "$CI_BRANCH" != "" ]]; then
          coveralls || echo "Failed to upload code coverage stats"
       else
          echo "Skipping coveralls upload: tests not passing or \$CI_COMMIT not set."
       fi
    fi

    exit ${exitstat}
}


check_if_tests_should_run
run_tests
