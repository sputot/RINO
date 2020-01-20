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

compare_to_ref=true  # when true, compare to results stored in output_0_xx ; when false, run and store results of ref_version

for i in "${!examples_indexes[@]}"
do
echo "******* Running RINO regression on ODE example no ${examples_indexes[i]} **************"
cd $rino
start_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
`$rino_exec 0 ${examples_indexes[i]} > /dev/null 2>&1`
end_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
elapsed_ms_rino=$((end_ms - start_ms))
echo "Execution time for RINO: $elapsed_ms_rino milliseconds"

# run RINO_REF only if we do not compare to stored version
if [ "$compare_to_ref" = false ]
then
cd $ref
start_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
`$ref_exec 0 ${examples_indexes[i]} > /dev/null 2>&1`
end_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
elapsed_ms_ref=$((end_ms - start_ms))
echo "Execution time for REF: $elapsed_ms_ref milliseconds"
if [ $(($elapsed_ms_rino)) -gt $((2*$elapsed_ms_ref)) ]
then
    echo "Warning: Execution time getting worse"
elif [ $(($elapsed_ms_ref)) -gt $((2*$elapsed_ms_rino)) ]
then
    echo "Congratulations: execution time improved"
fi
fi

# cp -r $ref"output/" $ref"output_0_"${examples_indexes[i]}"/"

if [ "$compare_to_ref" = false ]
then
    rm -r $ref"output_0_"${examples_indexes[i]}; cp -r $ref"output" $ref"output_0_"${examples_indexes[i]}
else
    rm -r $ref"output"; cp -r $ref"output_0_"${examples_indexes[i]} $ref"output"
fi

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
            : # empty command instead of echo "$file_rino and $file_ref are the same file"
        elif [ $error -eq 1 ]
#        if [ $error -eq 1 ]
        then
        results_changed=true
        tail -1 $file_rino  # print last line of file
        tail -1 $file_ref
        echo "$file_rino and $file_ref differ"
        else
        results_changed=true
        echo "Error: there was something wrong with the diff command between $file_rino and $file_ref "
        fi
    done
done
if [ "$results_changed" = true ]
then
    echo "Warning: results changed when running RINO on 0 ${examples_indexes[i]}"
else
    echo "No change when running RINO on 0 ${examples_indexes[i]}"
fi

done


# similarly for DDE examples
examples_indexes=(1 2 3 4 5 6 7 8 9 10 11)   # indexes of ODE examples we wish to test for non regression
sysdim=(1 2 7 1 1 2 4 2 1 9 19)

for i in "${!examples_indexes[@]}"
do
echo "******* Running RINO regression on DDE example no ${examples_indexes[i]} **************"
cd $rino
start_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
`$rino_exec 1 ${examples_indexes[i]} > /dev/null 2>&1`
end_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
elapsed_ms=$((end_ms - start_ms))
echo "Execution time for RINO: $elapsed_ms milliseconds"

# run RINO_REF only if we do not compare to stored version
if [ "$compare_to_ref" = false ]
then
cd $ref
start_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
`$ref_exec 1 ${examples_indexes[i]} > /dev/null 2>&1`
end_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
elapsed_ms=$((end_ms - start_ms))
echo "Execution time for REF: $elapsed_ms milliseconds"

if [ $(($elapsed_ms_rino)) -gt $((2*$elapsed_ms_ref)) ]
then
echo "Warning: Execution time getting worse"
elif [ $(($elapsed_ms_ref)) -gt $((2*$elapsed_ms_rino)) ]
then
echo "Congratulations: execution time improved"
fi
fi

#cp -r $ref"output/" $ref"output_1_"${examples_indexes[i]}"/"

if [ "$compare_to_ref" = false ]
then
    rm -r $ref"output_1_"${examples_indexes[i]}; cp -r $ref"output" $ref"output_1_"${examples_indexes[i]}
else
    rm -r $ref"output"; cp -r $ref"output_1_"${examples_indexes[i]} $ref"output"
fi

results_changed=false

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
: #echo "$file_rino and $file_ref are the same file"
elif [ $error -eq 1 ]
then
results_changed=true
tail -1 $file_rino  # print last line of file
tail -1 $file_ref
echo "$file_rino and $file_ref differ"
else
results_changed=true
echo "Error: there was something wrong with the diff command between $file_rino and $file_ref "
fi
done
done

if [ "$results_changed" = true ]
then
echo "Warning: results changed when running RINO on 1 ${examples_indexes[i]}"
else
echo "No change when running RINO on 1 ${examples_indexes[i]}"
fi


done

