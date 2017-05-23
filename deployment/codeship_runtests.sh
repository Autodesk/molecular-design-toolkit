#!/usr/bin/env bash


function shouldruntests(){
   if [ "${TESTENV}" == "complete" ]; then
       echo true
   fi

   case "${CI_BRANCH}" in
       master|deploy|dev)
       echo true
       ;;
   esac

   if [[ "${CI_COMMIT_MESSAGE}" == *" --testall "* ]];  then
       echo true
   fi

   echo false
}

if [ $(shouldruntests) == false ]; then
    echo "SKIPPING minimal environment tests for this commit."
    echo "To run these tests, add \"--testall\" to your commit message"
    echo "(or work in the dev or deploy branches)"
    python ../../deployment/send_test_status.py 0 "Skipped (dev/deploy branches only)"
    exit 0
else
    py.test -n 4 --junit-xml=/opt/reports/junit.${TESTENV}.xml | tee /opt/reports/pytest.${TESTENV}.log
    exitstat=${PIPESTATUS[0]}

    statline="$(tail -n1 /opt/reports/pytest.${TESTENV}.log)"

    echo 'Test status:'
    echo ${statline}

    python ../../deployment/send_test_status.py ${exitstat} "${statline}"

    exit ${exitstat}
fi
