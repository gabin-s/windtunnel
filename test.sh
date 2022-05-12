#!/bin/sh

EXEC_REF_PATH="./wind_seq"
EXEC_PATH="./wind_seq"
TEST_CASES="test_cases"

run_test() {
    echo -n "$test_name: ..."

    start=`date +%s`

    # run the actual test
    res=$(echo $args | xargs "$EXEC_PATH" | grep "Result: ")
    res=${res#Result: }

    end=`date +%s`
    runtime=$((end-start)) # time measurement

    if [ "$res" = "$expected" ]; then
        status="OK"
    else
        status="FAIL"
    fi

    printf "\r$test_name: $status (${runtime}s)\n"
}

echo "--- Running tests --"

while read line || [ -n "$line" ]; do
    
    case $line in
        "A:"*)
            args=${line#A:};;
        "-"*)
            test_name=${line#-};;
        "E:"*) 
            expected=${line#E:}
            run_test $args $expected;;
    esac

done < "$TEST_CASES"
