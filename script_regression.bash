#!/bin/bash

rino="/Users/sylvie/RINO/"
ref="/Users/sylvie/RINO_NONREGRESSION/RINO/"
rino_exec=$rino"main"
ref_exec=$ref"main"
rino_output=$rino"output/x"
ref_output=$ref"output/x"

# ODE examples
examples_indexes=(1 2 3 4 5 6 7 18)   # indexes of ODE examples we wish to test for non regression
sysdim=(1 2 4 5 2 4 4 14)


for i in "${!examples_indexes[@]}"
do
echo "******* Running RINO regression on ODE example no ${examples_indexes[i]} **************"
cd $rino
`$rino_exec 0 ${examples_indexes[i]} > /dev/null 2>&1`
cd $ref
`$ref_exec 0 ${examples_indexes[i]} > /dev/null 2>&1`

#for value in {1..2}
for ((value=1;value<=${sysdim[i]};value++));
do
    for suffix in "center.out" "inner.out" "inner_minimal.out" "inner_robust.out" "outer.out" "outer_minimal.out" "outer_robust.out"
    do
        filename=$value$suffix
        file_rino=$rino_output$filename
            file_ref=$ref_output$filename
        diff $file_rino $file_ref > /dev/null 2>&1
        error=$?
        if [ $error -eq 0 ]
        then
        echo "$file_rino and $file_ref are the same file"
        elif [ $error -eq 1 ]
        then
        echo "$file_rino and $file_ref differ"
        else
        echo "There was something wrong with the diff command between $file_rino and $file_ref "
        fi
    done
done

done


# similarly for DDE examples
examples_indexes=(1 2 3 4 5 6 7 8 9 10 11)   # indexes of ODE examples we wish to test for non regression
sysdim=(1 2 7 1 1 2 4 2 1 9 19)


for i in "${!examples_indexes[@]}"
do
echo "******* Running RINO regression on DDE example no ${examples_indexes[i]} **************"
cd $rino
`$rino_exec 1 ${examples_indexes[i]} > /dev/null 2>&1`
cd $ref
`$ref_exec 1 ${examples_indexes[i]} > /dev/null 2>&1`

#for value in {1..2}
for ((value=1;value<=${sysdim[i]};value++));
do
for suffix in "center.out" "inner.out" "inner_minimal.out" "inner_robust.out" "outer.out" "outer_minimal.out" "outer_robust.out"
do
filename=$value$suffix
file_rino=$rino_output$filename
file_ref=$ref_output$filename
diff $file_rino $file_ref > /dev/null 2>&1
error=$?
if [ $error -eq 0 ]
then
echo "$file_rino and $file_ref are the same file"
elif [ $error -eq 1 ]
then
echo "$file_rino and $file_ref differ"
else
echo "There was something wrong with the diff command between $file_rino and $file_ref "
fi
done
done

done

