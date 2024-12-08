make

for src in *.c; do
    test=$(basename $src .c)

    ./$test
    if [ $? -eq 0 ]; then
        result="pass"
    else
        result="fail"
    fi

    valgrind -q --error-exitcode=1 --leak-check=full ./$test
    if [ $? -eq 1 ]; then
        result="fail"
    fi

    echo $test: $result
done
