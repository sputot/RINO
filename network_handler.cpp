#include <iostream>
#include <vector>
#include <mutex>
#include <thread>

#include "network_handler.h"

std::mutex mtx2;
network_handler NH;
//vector<Layer> L(max_nb_layers);
computation_graph CG;

network_handler :: network_handler(const char* name)
{
    build_from_file(name);
}



// parser for the older sherlock format
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
        if(act.compare("Tanh") == 0 || act.compare("tanh") == 0){
            activations.push_back(ACT_TANH);
        }
        else if(act.compare("Sigmoid") == 0 || act.compare("sigmoid") == 0){
            activations.push_back(ACT_SIGMOID);
        }
        else if(act.compare("Linear") == 0 || act.compare("linear") == 0){
            activations.push_back(ACT_LINEAR);
        }
        else if(act.compare("Relu") == 0 || act.compare("relu") == 0){
            activations.push_back(ACT_RELU);
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
    
    
    // building the layers
    for (int i=0 ; i<n_hidden_layers+1 ; i++ )
        L.push_back(build_layer(i));
        //        L[i] = Layer(NH,i);
    
}

Layer network_handler::build_layer(int no_layer)
{
  //  cout << actual_weights[no_layer][0].size() << " " << actual_biases[no_layer].size() << endl;
    Layer res = Layer(actual_weights[no_layer][0].size(),actual_biases[no_layer].size());
    res.activation = activations[no_layer];
    res.nb_inputs = actual_weights[no_layer][0].size();
    res.nb_outputs = actual_biases[no_layer].size();
    cout << "nb_inputs layer "<< no_layer << " = " << res.nb_inputs << " nb_outputs = " << res.nb_outputs << " activation function = " << res.activation <<  endl;
    
    for (int i=0 ; i<res.nb_outputs; i++) {
        for (int j=0 ; j<res.nb_inputs; j++) {
            res.weights[i][j] = actual_weights[no_layer][i][j];
        }
    }
    for (int i=0 ; i<res.nb_outputs; i++)
        res.biases [i] = actual_biases[no_layer][i];

    return res;
}


// parser for ONNX : copied from evaluate_graph in computation_graph.cpp  (but external to the class => just redundant so useful only as a basis for the abstract version which is in the header file as it operates on templates)
void computation_graph_evaluate_graph( computation_graph CG, map < uint32_t, double > input_node_and_value,
                                      map < uint32_t, double > & output_node_and_value )
{
    // Making sure you received some non empty input map
    assert(!(input_node_and_value.empty()));
    // Making sure that the counts are consistent
    assert(CG.no_of_input_nodes == input_node_and_value.size());
    
    map< uint32_t, double > memoized_table;
    
    // Set the value of the input nodes and make sure those are constant nodes
    
    for(auto input_node : input_node_and_value)
    {
        assert( CG.all_nodes[input_node.first].return_node_type() == const_string );
        CG.all_nodes[input_node.first].set_node_val(input_node.second);
        memoized_table[input_node.first] = input_node.second;
        
        cout << "Setting node val of node id " << CG.all_nodes[input_node.first].get_node_number() << " as " << CG.all_nodes[input_node.first].return_current_output() << endl;
        
        }
        
        // Creating the space for the outputs of the network
        output_node_and_value.clear();
        for(auto output_node : CG.output_nodes)
        {
            output_node_and_value.insert(make_pair(output_node, 0.0));
        }
        
        
        for(auto & output_values : output_node_and_value)
        {
            computation_graph_evaluate_node_seq(CG, output_values.first, memoized_table, output_values.second);
        }
        
        
}






// copied from evaluate_node from computation_graph.cc as a starting point to implement the abstract version (not udeful in itself)
void computation_graph_evaluate_node(computation_graph & c_graph, uint32_t node_id ,
                   map< uint32_t , double > & table,
                   int & available_threads, double & ret_val , int thread_id)
{
    
    
    while(!mtx2.try_lock());
    auto & current_node = c_graph.all_nodes[node_id];
    map< uint32_t , pair< node * , datatype > > backward_connections_, backward_connections;
    current_node.get_backward_connections(backward_connections_);
    backward_connections = backward_connections_;
    mtx2.unlock();
    
    
    vector< thread > vector_of_threads_created;
    
    
    bool direction;
    if(thread_id % 2)
    {
        direction = true;
    }
    else
    {
        direction = false;
    }
    
    // for(auto some_connection : backward_connections)
    if(direction == true)
    {
        
        for(auto it = backward_connections.begin(); it != backward_connections.end(); it ++)
        {
            // auto input_node_ptr = some_connection.second.first;
            auto input_node_ptr = it->second.first;
            
            if(  input_node_ptr->return_node_type() == const_string) // If a constant node then just get the value and store in the table
            {
                // pair< uint32_t , double > node_and_value = make_pair(some_connection.first, input_node_ptr->return_current_output()) ;
                pair< uint32_t , double > node_and_value = make_pair(it->first, input_node_ptr->return_current_output()) ;
                while(!mtx2.try_lock());
                table.insert(node_and_value);
                mtx2.unlock();
                continue;
            }
            
            // check if the value is already in the table
            while(!mtx2.try_lock());
            bool value_in_the_table = ( ( table.find(input_node_ptr->get_node_number())  == table.end() ) ? (false) : (true) ) ;
            mtx2.unlock();
            if( value_in_the_table )
            {
                while(!mtx2.try_lock());
                // pair< uint32_t , double > node_and_value = make_pair(some_connection.first, table[input_node_ptr->get_node_number()] ) ;
                pair< uint32_t , double > node_and_value = make_pair(it->first, table[input_node_ptr->get_node_number()] ) ;
                mtx2.unlock();
            }
            else // make  a recursive call to the inputs, get the value and compute
            {
                if(available_threads == 1)
                {
                    double buffer;
                    // evaluate_node(c_graph, some_connection.first, table, available_threads, buffer, thread_id);
                    computation_graph_evaluate_node(c_graph, it->first, table, available_threads, buffer, thread_id);
                    // pair< uint32_t, double > node_and_value = make_pair(some_connection.first, buffer);
                    pair< uint32_t, double > node_and_value = make_pair(it->first, buffer);
                    while(!mtx2.try_lock());
                    table.insert(node_and_value);
                    mtx2.unlock();
                }
                else
                {
                    
                    while(!mtx2.try_lock());
                    available_threads--;
                    mtx2.unlock();
                    
                    // ci-dessous etait dans un debug eval, enlever ?
                        while(!mtx2.try_lock());
                        // cout << "Starting a thread from node number = " << some_connection.first << endl;
                        cout << "Starting a thread with node number = " << it->first << " from thread " << thread_id << endl;
                        cout << "Available threads = " << available_threads << endl;
                        mtx2.unlock();
                    
                    
                    double buffer;
                    thread current_thread(computation_graph_evaluate_node,
                                          ref(c_graph),
                                          // some_connection.first,
                                          it->first,
                                          ref(table),
                                          ref(available_threads),
                                          ref(buffer),
                                          // some_connection.first) ;
                                          it->first) ;
                    
                    vector_of_threads_created.push_back(move(current_thread));
                }
                
            }
            
        }
        
        
    }
    else
    {
        for(auto it = backward_connections.rbegin(); it != backward_connections.rend(); it ++)
        {
            // auto input_node_ptr = some_connection.second.first;
            auto input_node_ptr = it->second.first;
            
            if(  input_node_ptr->return_node_type() == const_string) // If a constant node then just get the value and store in the table
            {
                // pair< uint32_t , double > node_and_value = make_pair(some_connection.first, input_node_ptr->return_current_output()) ;
                pair< uint32_t , double > node_and_value = make_pair(it->first, input_node_ptr->return_current_output()) ;
                while(!mtx2.try_lock());
                table.insert(node_and_value);
                mtx2.unlock();
                continue;
            }
            
            // check if the value is already in the table
            while(!mtx2.try_lock());
            bool value_in_the_table = ( ( table.find(input_node_ptr->get_node_number())  == table.end() ) ? (false) : (true) ) ;
            mtx2.unlock();
            if( value_in_the_table )
            {
                while(!mtx2.try_lock());
                // pair< uint32_t , double > node_and_value = make_pair(some_connection.first, table[input_node_ptr->get_node_number()] ) ;
                pair< uint32_t , double > node_and_value = make_pair(it->first, table[input_node_ptr->get_node_number()] ) ;
                mtx2.unlock();
            }
            else // make  a recursive call to the inputs, get the value and compute
            {
                if(available_threads == 1)
                {
                    double buffer;
                    // evaluate_node(c_graph, some_connection.first, table, available_threads, buffer, thread_id);
                    computation_graph_evaluate_node(c_graph, it->first, table, available_threads, buffer, thread_id);
                    // pair< uint32_t, double > node_and_value = make_pair(some_connection.first, buffer);
                    pair< uint32_t, double > node_and_value = make_pair(it->first, buffer);
                    while(!mtx2.try_lock());
                    table.insert(node_and_value);
                    mtx2.unlock();
                }
                else
                {
                    while(!mtx2.try_lock());
                    available_threads--;
                    mtx2.unlock();
                    
                    // ci-dessous etait dasn un if debug-eval (enlever?)
                        while(!mtx2.try_lock());
                        // cout << "Starting a thread from node number = " << some_connection.first << endl;
                        cout << "Starting a thread with node number = " << it->first << " from thread " << thread_id << endl;
                        cout << "Available threads = " << available_threads << endl;
                        mtx2.unlock();
                    
                    
                    double buffer;
                    thread current_thread(computation_graph_evaluate_node,
                                          ref(c_graph),
                                          // some_connection.first,
                                          it->first,
                                          ref(table),
                                          ref(available_threads),
                                          ref(buffer),
                                          // some_connection.first) ;
                                          it->first) ;
                    
                    vector_of_threads_created.push_back(move(current_thread));
                }
                
            }
            
        }
        
    }
    
    
    
    int threads_counter = 0;
    
    for (thread & some_thread : vector_of_threads_created)
    {
        if (some_thread.joinable())
        {
            some_thread.join();
            while(!mtx2.try_lock());
            available_threads++;
            mtx2.unlock();
            
            // ci-dessous etait dasn un if debug-eval (enlever ou laisser, cf les autres blocs?)
            /*
            if(debug_eval)
            {
                while(!mtx2.try_lock());
                cout << "Some thread ended" << endl;
                cout << "Available threads = " << available_threads << endl;
                mtx2.unlock();
            }
             */
        }
        
    }
    
    // Since the input to all the nodes is ready, compute the output now
    
    map< uint32_t, double > inputs_to_the_node;
    for(auto some_connection : backward_connections)
    {
        auto input_node_ptr = some_connection.second.first;
        
        while(!mtx2.try_lock());
        bool value_in_the_table = ( ( table.find(input_node_ptr->get_node_number())  == table.end() ) ? (false) : (true) ) ;
        assert( value_in_the_table );
        pair< uint32_t , double > node_and_value = make_pair(some_connection.first, table[input_node_ptr->get_node_number()] ) ;
        mtx2.unlock();
        
        inputs_to_the_node.insert(node_and_value);
    }
    
    
    
    while(!mtx2.try_lock());
    
    current_node.set_inputs(inputs_to_the_node);
    double result = current_node.return_current_output();
    table.insert( make_pair ( node_id , result ) );
   // if(debug_eval)
   // {
        cout << "Computed value of node_id : " << node_id << " as " << result << " in thread id = " << thread_id <<  endl;
        cout << "Current table : " << " [ " ;
        for(auto each_entry : table)
        {
            cout << each_entry.first << " --- " << each_entry.second << " , ";
        }
        cout << " ] " << endl;
        
   // }
    mtx2.unlock();
    
    ret_val = result;
    
    return;
    
}

// just a simpler version of computation_graph_evaluate_node without using threads (not used in itslef, just as a test version to derive the abstract version)
void computation_graph_evaluate_node_seq(computation_graph & c_graph, uint32_t node_id ,
                                     map< uint32_t , double > & table, double & ret_val)
{
    
    auto & current_node = c_graph.all_nodes[node_id];
    map< uint32_t , pair< node * , datatype > > backward_connections_, backward_connections;
    current_node.get_backward_connections(backward_connections_);
    backward_connections = backward_connections_;
    
    
        for(auto it = backward_connections.rbegin(); it != backward_connections.rend(); it ++)
        {
            // auto input_node_ptr = some_connection.second.first;
            auto input_node_ptr = it->second.first;
            
            if(  input_node_ptr->return_node_type() == const_string) // If a constant node then just get the value and store in the table
            {
                // pair< uint32_t , double > node_and_value = make_pair(some_connection.first, input_node_ptr->return_current_output()) ;
                pair< uint32_t , double > node_and_value = make_pair(it->first, input_node_ptr->return_current_output()) ;
                table.insert(node_and_value);
                continue;
            }
            
            // check if the value is already in the table
            bool value_in_the_table = ( ( table.find(input_node_ptr->get_node_number())  == table.end() ) ? (false) : (true) ) ;
            if( value_in_the_table )
                pair< uint32_t , double > node_and_value = make_pair(it->first, table[input_node_ptr->get_node_number()] ) ;
            else // make  a recursive call to the inputs, get the value and compute
            {
                    double buffer;
                    // evaluate_node(c_graph, some_connection.first, table, available_threads, buffer, thread_id);
                    computation_graph_evaluate_node_seq(c_graph, it->first, table, buffer);
                    // pair< uint32_t, double > node_and_value = make_pair(some_connection.first, buffer);
                    pair< uint32_t, double > node_and_value = make_pair(it->first, buffer);
                    table.insert(node_and_value);
            }
        }
    
    
    // Since the input to all the nodes is ready, compute the output now
    
    map< uint32_t, double > inputs_to_the_node;
    for(auto some_connection : backward_connections)
    {
        auto input_node_ptr = some_connection.second.first;
        
        bool value_in_the_table = ( ( table.find(input_node_ptr->get_node_number())  == table.end() ) ? (false) : (true) ) ;
        assert( value_in_the_table );
        pair< uint32_t , double > node_and_value = make_pair(some_connection.first, table[input_node_ptr->get_node_number()] ) ;
        
        inputs_to_the_node.insert(node_and_value);
    }
    
    
    
    
    current_node.set_inputs(inputs_to_the_node);
    double result = current_node.return_current_output();
    table.insert( make_pair ( node_id , result ) );
    // if(debug_eval)
    // {
    cout << "Computed value of node_id : " << node_id << " as " << result <<  endl;
    cout << "Current table : " << " [ " ;
    for(auto each_entry : table)
    {
        cout << each_entry.first << " --- " << each_entry.second << " , ";
    }
    cout << " ] " << endl;
    
    // }
    
    ret_val = result;
    
    return;
    
}



void test_network_sigmoid_cav(computation_graph & CG)
{
    CG.clear();
    // The two input nodes to the graph declared as constants
    node node_1(1, "constant");
    CG.add_new_node(1, node_1);
    node node_2(2, "constant");
    CG.add_new_node(2, node_2);
    
    // The internal nodes
    node node_3(3, "sigmoid");
    CG.add_new_node(3, node_3);
    node node_4(4, "sigmoid");
    CG.add_new_node(4, node_4);
    
    // The output node
    node node_5(5, "sigmoid");
    CG.add_new_node(5, node_5);
    
    
    // First let's mark some of the nodes as inputs and outputs
    CG.mark_node_as_input(1);
    CG.mark_node_as_input(2);
    CG.mark_node_as_output(5);
    
    // Now let's create the connections:
    
    // first layer connections and bias
    CG.connect_node1_to_node2_with_weight(1,3,2.0);
    CG.connect_node1_to_node2_with_weight(1,4,-2.0);
    CG.connect_node1_to_node2_with_weight(2,3,-1.0);
    CG.connect_node1_to_node2_with_weight(2,4,3.);
    CG.set_bias_of_node(3, -0.4);
    CG.set_bias_of_node(4, -0.5);
    
    CG.connect_node1_to_node2_with_weight(3,5,2.0);
    CG.connect_node1_to_node2_with_weight(4,5,1.0);
    CG.set_bias_of_node(5,-1.0);
    
}



/*
abstract_node::abstract_node()
{
   
    neuron_bounds.clear();
    node_bounds.clear();
}
 */
