#! /bin/sh
# TORA
./rino -configfile Examples/ConfigFiles/cfg_tora_tanh.txt
mv output output_tora_tanh
./rino -configfile Examples/ConfigFiles/cfg_tora_sigmoid.txt
mv output output_tora_sigmoid

# B1
./rino -configfile Examples/ConfigFiles/cfg_B1_tanh.txt
mv output output_B1_tanh
./rino -configfile Examples/ConfigFiles/cfg_B1_sigmoid.txt
mv output output_B1_sigmoid

# B2
./rino -configfile Examples/ConfigFiles/cfg_B2_sigmoid.txt
mv output output_B2_sigmoid

# B3
./rino -configfile Examples/ConfigFiles/cfg_B3_tanh.txt
mv output output_B3_tanh
./rino -configfile Examples/ConfigFiles/cfg_B3_sigmoid.txt
mv output output_B3_sigmoid

# B4:
./rino -configfile Examples/ConfigFiles/cfg_B4_tanh.txt
mv output output_B4_tanh
./rino -configfile Examples/ConfigFiles/cfg_B4_sigmoid.txt
mv output output_B4_sigmoid

# B5:
./rino -configfile Examples/ConfigFiles/cfg_B5_tanh.txt
mv output output_B5_tanh
./rino -configfile Examples/ConfigFiles/cfg_B5_sigmoid.txt
mv output output_B5_sigmoid

# ACC:
./rino -configfile Examples/ConfigFiles/cfg_acc_tanh.txt
mv output output_ACC_tanh

# Continuous-time Mountain Car
./rino -configfile Examples/ConfigFiles/cfg_MC_sigmoid.txt
mv output output_MC_sigmoid

# Discrete-time Mountain Car
./rino -configfile Examples/ConfigFiles/cfg_discrete_mc.txt
mv output output_discrete_mc

