#!/usr/bin/env bash

# this script expects to run from the moldesign/_tests directory

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
     python ../../deployment/send_test_status.py 0 "Skipped (dev/deploy branches only)"
     exit 0
   fi
}


function run_tests(){
	if [ "${TESTENV}" == "complete" ]; then
	   coverageflags="--cov .. --cov-report=term-missing:/opt/reports/coverage --cov-report="
	fi

	py.test -n 4 --junit-xml=/opt/reports/junit.${TESTENV}.xml $coverageflags | tee /opt/reports/pytest.${TESTENV}.log
	exitstat=${PIPESTATUS[0]}

	statline="$(tail -n1 /opt/reports/pytest.${TESTENV}.log)"

	echo 'Test status:'
	echo ${statline}

	python ../../deployment/send_test_status.py ${exitstat} "${statline}"

	if [ "${TESTENV}" == "complete" ]; then
       if [ ${exitstat} -eq 0 ]; then
          coveralls || echo "Failed to upload code coverage stats"
       else
       	  echo "Skipping coveralls upload: tests not passing."
       	  coverage report
       fi
    fi

	exit ${exitstat}
}


check_if_tests_should_run
run_tests
