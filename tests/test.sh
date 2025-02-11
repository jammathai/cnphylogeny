#!/bin/bash

make
if [ $? -ne 0 ]; then
    exit
fi

for src in *.c; do
    test=$(basename $src .c)

    ./$test
    if [ $? -eq 0 ]; then
        result="pass"
        valgrind -q --error-exitcode=1 --leak-check=full ./$test
        if [ $? -eq 1 ]; then
            result="fail"
        fi
    else
        result="fail"
    fi

    echo $test: $result
done
