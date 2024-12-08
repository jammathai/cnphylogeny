for src in *.c; do
    test=$(basename $src .c)
    ./$test
    if [ $? -eq 0 ]; then
        result="pass"
    else
        result="fail"
    fi
    echo $test: $result
done
