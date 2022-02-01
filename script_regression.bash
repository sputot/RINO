#!/bin/bash

rino="/Users/sylvie/RINO/"
ref="/Users/sylvie/RINO_NONREGRESSION/RINO/"
rino_exec=$rino"rino"
ref_exec=$ref"rino"
rino_output=$rino"output"
ref_output=$ref"output"

test_ode=true
test_ode_nn=true
test_dde=true
test_discrete=true

if [ $# == 1 ]
then
  if [ "$1" == "all" ]
  then
    test_ode=true
    test_ode_nn=true
    test_dde=true
    test_discrete=true
  elif [ "$1" == "ode" ]
  then
    test_ode=true
    test_ode_nn=false
    test_dde=false
    test_discrete=false
  elif [ "$1" == "ode-nn" ]
  then
    test_ode=true
    test_ode_nn=true
    test_dde=false
    test_discrete=false
  elif [ "$1" == "dde" ]
  then
    test_ode=false
    test_ode_nn=false
    test_dde=true
    test_discrete=false
  elif [ "$1" == "discrete" ]
  then
    test_ode=false
    test_ode_nn=false
    test_dde=false
    test_discrete=true
  fi
fi


yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}


# examples (defined with config files)
examples_ode=(cfg_ode_2.txt cfg_ode_6.txt cfg_ode_18.txt)
examples_ode_nn=(cfg_B1_sigmoid.txt cfg_B1_tanh.txt cfg_B2_sigmoid.txt cfg_B2_tanh.txt cfg_B3_sigmoid.txt cfg_B3_tanh.txt cfg_B4_sigmoid.txt cfg_B4_tanh.txt cfg_B5_sigmoid.txt cfg_B5_tanh.txt cfg_acc_tanh.txt cfg_tora_sigmoid.txt cfg_tora_tanh.txt)
examples_dde=()
examples_discrete=(cfg_discrete_15.txt cfg_discrete_16_1.txt cfg_discrete_16_2.txt cfg_discrete_17.txt)



# when true, compare to results stored in output_0_xx ; when false, run and store results of ref_version (=> set to false when new ref version, true otherwise)
compare_to_ref=false

if [ "$test_discrete" == true ]
then
  echo "******* Running RINO regression on discrete examples **************"
 for i in "${!examples_discrete[@]}"
 do
    cd $rino
    config_file=$rino"Examples/ConfigFiles/"${examples_discrete[i]}
    echo "*** Running: rino -configfile $config_file ***"
    start_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
    `$rino_exec -configfile $config_file > /dev/null 2>&1`
    end_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
    elapsed_ms_rino=$((end_ms - start_ms))
    echo "Execution time for RINO: $elapsed_ms_rino milliseconds"
    
    # run RINO_REF only if we do not compare to stored version
    if [ "$compare_to_ref" = false ]
    then
        cd $ref
        start_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
        `$ref_exec -configfile $config_file > /dev/null 2>&1`
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
    
    if [ "$compare_to_ref" = false ]
    then
        rm -r $ref"output_"${examples_discrete[i]}; cp -r $ref"output" $ref"output_"${examples_discrete[i]}
    else
        rm -r $ref"output"; cp -r $ref"output_"${examples_discrete[i]} $ref"output"
    fi

    results_changed=false
    

    file_rino=${rino_output}/sumup.yaml
    file_ref=${ref_output}/sumup.yaml
    
    zouter_rino=$(yaml $file_rino "['zouter']")
    zouter_ref=$(yaml $file_ref "['zouter']")
    #echo $zouter_rino
    #echo $zouter_ref
    
    diff ${rino_output}/sumup.yaml ${ref_output}/sumup.yaml > /dev/null 2>&1
    error=$?
    if [ $error -eq 0 ]
    then
        : # empty command instead of echo "$file_rino and $file_ref are the same file"
    elif [ $error -eq 1 ]
#       if [ $error -eq 1 ]
    then
        results_changed=true
        echo "$file_rino and $file_ref differ"
        echo "zouter_rino is" $zouter_rino
        echo "zouter_ref is" $zouter_ref
    else
        results_changed=true
        echo "Error: there was something wrong with the diff command between $file_rino and $file_ref "
    fi


    if [ "$results_changed" = true ]
    then
        echo "Warning: results changed when running RINO on ${examples_discrete[i]}"
    else
        echo "No change when running RINO on ${examples_discrete[i]}"
    fi
    
 done
fi







if [ "$test_ode" == true ]
then
for i in "${!examples_indexes[@]}"
do
echo "******* Running RINO regression on ODE example no ${examples_indexes[i]} **************"
cd $rino
start_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
config_file=$rino"examples/config_0_"${examples_indexes[i]}".txt"
if test -f $config_file
then
    `$rino_exec -systype ode -syschoice ${examples_indexes[i]} -configfile $config_file > /dev/null 2>&1`
else
    `$rino_exec -systype ode -syschoice ${examples_indexes[i]} > /dev/null 2>&1`
fi
end_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
elapsed_ms_rino=$((end_ms - start_ms))
echo "Execution time for RINO: $elapsed_ms_rino milliseconds"

# run RINO_REF only if we do not compare to stored version
if [ "$compare_to_ref" = false ]
then
cd $ref
start_ms=$(ruby -e 'puts (Time.now.to_f * 1000).to_i')
if test -f $config_file
then
    `$ref_exec -systype ode -syschoice ${examples_indexes[i]} -configfile $config_file > /dev/null 2>&1`
else
    `$ref_exec -systype ode -syschoice ${examples_indexes[i]} > /dev/null 2>&1`
fi
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

results_changed=false

#for value in {1..2}
for ((value=0;value<${sysdim[i]};value++));
do
    for suffix in "center.out" "inner.out" "inner_minimal.out" "inner_robust.out" "outer.out" "outer_minimal.out" "outer_robust.out"
    do
        filename=$value$suffix
        file_rino=$rino_output$filename
        file_ref=$ref_output$filename
        if test -f $file_rino || test -f $file_ref
        then
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
fi

