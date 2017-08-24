#!/usr/bin/env bash
# this script expects to run from the root of the repository

set -e  # fail immediately if any command fails


VERSION="${TESTENV}.py${PYVERSION}"
PYTESTFLAGS="-n 3 --spec  --durations=20 --junit-xml=/opt/reports/junit.${VERSION}.xml --timeout=3600 --tb=short"
if [ "${VERSION}" == "complete.py3" ]; then
       PYTESTFLAGS="--cov moldesign ${PYTESTFLAGS}"
fi


function send_status_update(){
     python deployment/send_test_status.py "${1}" "${2}"
}


function check_if_tests_should_run(){
   echo "Should I run the tests in this environment?"
   runthem=false

   if [[ "${CI_COMMIT_MESSAGE}" == *"--fast-ci-tests"* && "${VERSION}" != "complete.py3" ]];  then
       echo "NO: found \"--fast-ci-tests\" flag in commit message; run complete.py3 only"
       exit 0
   fi

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
     exit 0
   fi
}


function run_tests(){
    send_status_update "na" "Starting tests for ${VERSION}"

    echo
    echo "Test command running in working dir '$(pwd)':"
    echo "py.test ${PYTESTFLAGS}"
    echo

    py.test ${PYTESTFLAGS} | tee /opt/reports/pytest.${VERSION}.log
    exitstat=${PIPESTATUS[0]}

    statline="$(tail -n1 /opt/reports/pytest.${VERSION}.log)"

    echo 'Test status:'
    echo ${statline}

    send_status_update "${exitstat}" "${statline}"

    if [ "${VERSION}" == "complete.py3" ]; then
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
