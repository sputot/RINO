#!/bin/bash

rino="/Users/sylvie/RINO/"
ref="/Users/sylvie/RINO_NONREGRESSION/RINO/"
rino_exec=$rino"rino"
ref_exec=$ref"rino"
rino_output=$rino"output"
ref_output=$ref"output"


yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}

yaml_bis() {
    python3 -c "import yaml;data=yaml.safe_load(open('$1'));print(data['$2'] if '$2' in data else '')"
}


# examples (defined with config files)
examples_ode=(cfg_ode_2.txt cfg_ode_6.txt cfg_ode_18.txt cfg_ode_50.txt)
examples_ode_nn=(cfg_B1_sigmoid.txt cfg_B1_tanh.txt cfg_B2_sigmoid.txt cfg_B2_tanh.txt cfg_B3_sigmoid.txt cfg_B3_tanh.txt cfg_B4_sigmoid.txt cfg_B4_tanh.txt cfg_B5_sigmoid.txt cfg_B5_tanh.txt cfg_acc_tanh.txt cfg_tora_sigmoid.txt cfg_tora_tanh.txt)
examples_dde=(cfg_dde_1.txt cfg_dde_3.txt cfg_dde_6.txt cfg_dde_8.txt cfg_dde_10.txt cfg_dde_11.txt)
examples_discrete=(cfg_discrete_15.txt cfg_discrete_16_1.txt cfg_discrete_16_2.txt cfg_discrete_17.txt cfg_discrete_18.txt)


#UnixShell=("${Unix[@]}" "${Shell[@]}")


if [ $# == 1 ]
then
  if [ "$1" == "all" ]
  then
    examples=(${examples_ode_nn[@]} ${examples_discrete[@]} ${examples_ode[@]} ${examples_dde[@]})
  elif [ "$1" == "ode" ]
  then
    examples=(${examples_ode[@]})
  elif [ "$1" == "ode-nn" ]
  then
    examples=(${examples_ode_nn[@]})
  elif [ "$1" == "dde" ]
  then
    examples=(${examples_dde[@]})
  elif [ "$1" == "discrete" ]
  then
    examples=(${examples_discrete[@]})
  else
    examples=($1)
  fi
fi

    #echo ${examples[@]}


# when true, compare to results stored in output_0_xx ; when false, run and store results of ref_version (=> set to false when new ref version, true otherwise)
compare_to_ref=true


 echo "******* Running RINO regression on ${examples[@]} **************"
 for i in "${!examples[@]}"
 do
    cd $rino
    config_file=$rino"Examples/ConfigFiles/"${examples[i]}
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
        rm -r $ref"output_"${examples[i]}; cp -r $ref"output" $ref"output_"${examples[i]}
    else
        rm -r $ref"output"; cp -r $ref"output_"${examples[i]} $ref"output"
    fi

    results_changed=false
    

    file_rino=${rino_output}/sumup.yaml
    file_ref=${ref_output}/sumup.yaml
    
    zouter_rino=$(yaml $file_rino "['zouter']")
    zouter_ref=$(yaml $file_ref "['zouter']")
    zinner_rino=$(yaml $file_rino "['zinner']")
    zinner_ref=$(yaml $file_ref "['zinner']")
    elapsed_rino=$(yaml $file_rino "['elapsed-secs']")
    elapsed_ref=$(yaml $file_ref "['elapsed-secs']")

    # same as above but function yaml_bis handles the case when the key is not present
    zouter_rob_rino=$(yaml_bis $file_rino "zouter_rob")
    zouter_rob_ref=$(yaml_bis $file_ref "zouter_rob")
    zinner_rob_rino=$(yaml_bis $file_rino "zinner_rob")
    zinner_rob_ref=$(yaml_bis $file_ref "zinner_rob")

    if [[ "$zouter_rino" != "$zouter_ref" ]] || [[ "$zinner_rino" != "$zinner_ref" ]] || [[ "$zouter_rob_rino" != "$zouter_rob_ref" ]] || [[ "$zinner_rob_rino" != "$zinner_rob_ref" ]]
    then
        results_changed=true
        echo "Reachability results in $file_rino and $file_ref differ"
        if [[ $zouter_rino != $zouter_ref ]]
        then
            echo "zouter_rino is" $zouter_rino
            echo "zouter_ref is" $zouter_ref
        fi
        if [[ "$zinner_rino" != "$zinner_ref" ]]
        then
            echo "zinner_rino is" $zinner_rino
            echo "zinner_ref is" $zinner_ref
        fi
        if [[ $zouter_rob_rino != $zouter_rob_ref ]]
        then
            echo "zouter_rob_rino is" $zouter_rob_rino
            echo "zouter_rob_ref is" $zouter_rob_ref
        fi
        if [[ "$zinner_rob_rino" != "$zinner_rob_ref" ]]
        then
            echo "zinner_rob_rino is" $zinner_rob_rino
            echo "zinner_rob_ref is" $zinner_rob_ref
        fi
    fi
    
    echo "Execution time (sec) for analysis phase: rino=$elapsed_rino ref=$elapsed_ref"
    
    #diff ${rino_output}/sumup.yaml ${ref_output}/sumup.yaml > /dev/null 2>&1
    #error=$?
    #if [ $error -eq 0 ]
    #then
    #    : # empty command instead of echo "$file_rino and $file_ref are the same file"
    #elif [ $error -eq 1 ]
#       if [ $error -eq 1 ]
    #then
    #    results_changed=true
    #    echo "$file_rino and $file_ref differ"
    #    echo "zouter_rino is" $zouter_rino
    #    echo "zouter_ref is" $zouter_ref
    #    echo "zinner_rino is" $zinner_rino
    #    echo "zinner_ref is" $zinner_ref
    #else
    #    results_changed=true
    #    echo "Error: there was something wrong with the diff command between $file_rino and $file_ref "
    #fi

    if [ "$results_changed" = true ]
    then
        echo "Warning: results changed when running RINO on ${examples[i]}"
    else
        echo "No change when running RINO on ${examples[i]}"
    fi
    
 done

