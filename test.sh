#!/bin/sh


EXEC_REF_PATH="./wind_seq"
EXEC_PATH="./$1"
TEST_CASES="test_cases"

[ "$#" -ne "1" ]      && echo "usage: $0 <binary_path>" && exit 1
[ ! -f "$EXEC_PATH" ] && echo "Binary not found: $EXEC_PATH" && exit 1

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
            run_test;;
    esac

done < "$TEST_CASES"
