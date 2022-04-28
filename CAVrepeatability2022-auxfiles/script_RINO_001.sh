# Now with time step 0.01s

# TORA
./rino -configfile Examples/ConfigFiles/cfg_tora_tanh_001.txt
mv output output_tora_tanh_001
./rino -configfile Examples/ConfigFiles/cfg_tora_sigmoid_001.txt
mv output output_tora_sigmoid_001

# B1
./rino -configfile Examples/ConfigFiles/cfg_B1_tanh_001.txt
mv output output_B1_tanh_001
./rino -configfile Examples/ConfigFiles/cfg_B1_sigmoid_001.txt
mv output output_B1_sigmoid_001

# B2
./rino -configfile Examples/ConfigFiles/cfg_B2_sigmoid_001.txt
mv output output_B2_sigmoid_001

# B3
./rino -configfile Examples/ConfigFiles/cfg_B3_tanh_001.txt
mv output output_B3_tanh_001
./rino -configfile Examples/ConfigFiles/cfg_B3_sigmoid_001.txt
mv output output_B3_sigmoid_001

# B4:
./rino -configfile Examples/ConfigFiles/cfg_B4_tanh_001.txt
mv output output_B4_tanh_001
./rino -configfile Examples/ConfigFiles/cfg_B4_sigmoid_001.txt
mv output output_B4_sigmoid_001

# B5:
./rino -configfile Examples/ConfigFiles/cfg_B5_tanh_001.txt
mv output output_B5_tanh_001
./rino -configfile Examples/ConfigFiles/cfg_B5_sigmoid_001.txt
mv output output_B5_sigmoid_001

# ACC:
./rino -configfile Examples/ConfigFiles/cfg_acc_tanh_001.txt
mv output output_acc_tanh_001

