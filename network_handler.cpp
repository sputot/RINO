#include <iostream>
#include <vector>

#include "network_handler.h"

network_handler NH;
vector<Layer> L(max_nb_layers);

network_handler :: network_handler(const char* name)
{
    build_from_file(name);
}




void network_handler :: build_from_file(const char* name){
    ifstream file;
    file.open(name);
    
    // The expcted input format goes something like this,
    // no_of_inputs
    // no_of_outputs
    // no_of_layers
    // Activation function input layer
    // configuration no_1 :  no_of_neurons in layer 1
    // Activation function layer 1
    // configuration no_2 :  no_of_neurons in layer 2
    // configuration no_3 :  no_of_neurons in layer 3
    // configuration no_4 :  no_of_neurons in layer 4
    // ....
    // weight into the first neuron from the 1st input
    // weight into the first neuron from the 2nd input
    // ....
    // bias for the first neuron
    // weight into the 2nd neuron from the 1st input
    // ....
    // ....
    
    
    unsigned int i, j , k, size_1, size_2, buffer_integer ;
    file >> n_inputs;
    file >> n_outputs;
    file >> n_hidden_layers;
    double buff;
    string act;
    network_configuration.clear();
    i = 0;
    while(i < n_hidden_layers + 1)
    {
        file >> act;
        if(act.compare("Tanh") == 0){
            activations.push_back(ACT_TANH);
        }
        else if(act.compare("Sigmoid") == 0){
            activations.push_back(ACT_SIGMOID);
        }
        else if(act.compare("Linear") == 0){
            activations.push_back(ACT_LINEAR);
        }
        else{
            std::cout << "Failed to load network. Unknown activation : " << act << std::endl;
            exit(1);
        }
        if(i < n_hidden_layers){
            file >> buff;
            buffer_integer = (int)buff;
            network_configuration.push_back(buffer_integer);
        }
        i++;
    }
    
    double data;
    // reserving space for the number of hidden layers + 1
    
    vector<vector< double > > input_buffer_mat(network_configuration[0],vector< double > (n_inputs,0));
    vector< double > input_buffer_vec(network_configuration[0]);
    
    actual_weights.reserve(n_hidden_layers + 1);
    actual_biases.reserve(n_hidden_layers + 1);
    
    // reading the input matrix
    int count = 0;
    for(i = 0 ; i < network_configuration[0] ; i++)
    {
        for( j = 0; j < n_inputs; j++)
        {
            file >> data;
            input_buffer_mat[i][j] = data;
        }
        
        file >> data;
        input_buffer_vec[i] = data;
    }
    
    actual_weights.push_back(input_buffer_mat);
    actual_biases.push_back(input_buffer_vec);
    
    vector< double > buffer_vec;
    vector<vector< double > > buffer_weight_mat;
    vector< double > buffer_bias_vec;
    
    // Reading the inner matrices
    if(n_hidden_layers > 1)
    {
        for(i = 1; i < n_hidden_layers; i++)
        {
            buffer_weight_mat.clear();
            buffer_bias_vec.clear();
            
            for(j = 0 ; j < network_configuration[i]; j++)
            {
                buffer_vec.clear();
                for(k = 0; k < network_configuration[i-1] ; k++)
                {
                    file >> data;
                    buffer_vec.push_back(data);
                }
                file >> data;
                buffer_bias_vec.push_back(data);
                buffer_weight_mat.push_back(buffer_vec);
            }
            
            actual_weights.push_back(buffer_weight_mat);
            actual_biases.push_back(buffer_bias_vec);
        }
        
    }
    
    
    vector<vector< double > > outer_buffer_mat(n_outputs,
                                               vector< double > (network_configuration[n_hidden_layers-1],0));
    
    vector< double > outer_buffer_vec(n_outputs);
    
    
    // Reading the  output matrix
    for(i  = 0; i < n_outputs ; i++)
    {
        for(j = 0; j < network_configuration[n_hidden_layers - 1]; j++)
        {
            file >> data;
            outer_buffer_mat[i][j] = data;
        }
        file >> data;
        outer_buffer_vec[i] = data;
    }
    
    actual_weights.push_back(outer_buffer_mat);
    actual_biases.push_back(outer_buffer_vec);
    
    
    data = -100;
    file >> data;
    if(data != (-100))
    {
        cout << "Network input file probably has an error ! " << endl;
    }
    
    file.close();
}
