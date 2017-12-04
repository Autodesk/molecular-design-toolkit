#!/usr/bin/env bash
# Drives tests for our CI system. This looks for the following environment variables:
# Defined by codeship
# - CI_BRANCH
# - CI_COMMIT_MESSAGE
# - PROJECT_ID
# Defined in ../codeship-services.yml
# - TESTENV
# - PYVERSION

set -e  # fail immediately if any command fails

if [ -z "${CI_BRANCH}" ]; then
   echo "FAILURE: Variable \$CI_BRANCH not defined."
   exit 1
fi

install_location=$(python -c "import moldesign, os; print(moldesign.__path__[0])")
test_location=$(dirname "${install_location}")

VERSION="${TESTENV}.py${PYVERSION}"
PYTESTFLAGS="moldesign/_tests/ -n 2 --spec  --durations=20
        --junit-xml=/opt/reports/junit.${VERSION}.xml --timeout=3600 --tb=short
        --cov moldesign --cov-config /opt/molecular-design-toolkit/.coveragerc"


function send_status_update(){
     python /opt/molecular-design-toolkit/deployment/send-test-status.py "${1}" "${2}"
}


function check_if_tests_should_run(){
   echo "Should I run the tests in this environment?"

   if [[ "${CI_COMMIT_MESSAGE}" == *"--skip-ci-tests"* ]];  then
       echo "NO: found \"--skip-ci-tests\" flag in commit message; will not run any test suites"
       exit 0
   fi

   if [[ "${CI_COMMIT_MESSAGE}" == *"--fast-ci-tests"* && "${VERSION}" != "complete.py3" ]];  then
       echo "NO: found \"--fast-ci-tests\" flag in commit message; run complete.py3 only"
       exit 0
   fi

   if [[ "${CI_BRANCH}" =~ ^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)((a|rc|b)(0|[1-9]\d*))?$ ]]
   then
       echo "YES: this is a release version: \"${CI_BRANCH}\""
       return 0
   else  # otherwise, point to the appropriate docker image tag
       mkdir -p ~/.moldesign
       echo "default_version_tag: ${CI_BRANCH}" >> ~/.moldesign/moldesign.yml
   fi


   if [ "${TESTENV}" == "complete" ]; then
       echo "YES: always run in 'complete' environment"
       return 0
   fi

   case "${CI_BRANCH}" in
       master|deploy|dev)
       echo "YES: always run in branch \"${CI_BRANCH}\""
       return 0
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
    cd ${test_location}


    echo
    echo "Test command running in working dir '$(pwd)':"
    echo "py.test ${PYTESTFLAGS}"
    echo

    py.test ${PYTESTFLAGS} | tee /opt/reports/pytest.${VERSION}.log
    exitstat=${PIPESTATUS[0]}
    statline="$(tail -n1 /opt/reports/pytest.${VERSION}.log)"

    # Make a copy of the coverage report
    mkdir -p /opt/reports/env-coverage/
    cp .coverage /opt/reports/env-coverage/coverage.${VERSION}

    echo 'Test status:'
    echo ${statline}

    send_status_update "${exitstat}" "${statline}"

    exit ${exitstat}
}


check_if_tests_should_run
run_tests
