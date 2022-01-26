#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <map>
#include <queue>

#include "aa_aaf.h"

#define ONNX_active 0

#if ONNX_active
    #include "sherlock.h"
#endif



using namespace std;

#ifndef NETWORK_HANDLER_H
#define NETWORK_HANDLER_H


enum Activation {ACT_RELU, ACT_SIGMOID, ACT_TANH, ACT_LINEAR};

#define max_nb_layers 10

#define testmode false


class Layer
{
public:
    Activation activation;
    int nb_inputs;
    int nb_outputs;
    
    vector<double> biases; // dimension (nb_of_outputs);
    vector<vector<double>> weights; // dimension (nb_of_outputs,vector<double>(nb_of_inputs));
    
    Layer() {}
    
    Layer(int _nb_inputs, int _nb_outputs) : biases(_nb_outputs), weights(_nb_outputs,vector<double>(_nb_inputs))
    {
        nb_inputs = _nb_inputs;
        nb_outputs = _nb_outputs;
    }
    
    template <class C> vector<C> eval_layer(vector<C> x) {
        vector<C> z(nb_outputs);
        assert(nb_inputs == x.size());
        
        for (int i=0 ; i<nb_outputs; i++) {
            z[i] = 0;
            for (int j=0 ; j<nb_inputs; j++)
                z[i] += weights[i][j]*x[j];
            z[i] = eval_activation(activation,z[i]+biases[i]);
        }
        return z;
    }
    
    template <class C> vector<C> eval_linear_layer(vector<C> x) {
        vector<C> z(nb_outputs);
        assert(nb_inputs == x.size());
        
        for (int i=0 ; i<nb_outputs; i++) {
            z[i] = 0;
            for (int j=0 ; j<nb_inputs; j++)
                z[i] += weights[i][j]*x[j];
            z[i] = z[i]+biases[i]; // eval_activation(activation,z[i]+biases[i]);
        }
        return z;
    }
    
};





class network_handler
{
    public :
    
    char * name_of_file;
    // data structures for the neural net information
    
    vector < vector < vector < double > > > actual_weights;
    vector < vector < double > > actual_biases;
    // contains all the biases, the size of biases is one more than
    // the number of sets of weights
    
    // unsigned int no_of_inputs,no_of_outputs,no_of_hidden_layers;
    
    vector< unsigned int > network_configuration;
    std::vector< Activation > activations;
    
    vector<Layer> L;
    
    
    int n_inputs, n_outputs, n_hidden_layers;
    
    network_handler() {};
    
    network_handler(const char* name);
    // the constructor which takes in the information file
    void build_from_file(const char*);
    
    Layer build_layer(int no_layer);
    
    
    template <class C> vector<C> eval_network(vector<C> x) {
        
        if (testmode)
        {
            vector<C> res(x.size());
            for (int i=0 ; i<x.size(); i++)
                res[i] = x[i];
         //   cout << "res=" << res[0] << endl;
            return res;
        }
        
        vector<vector<C>> net_outputs(n_hidden_layers+2);
        net_outputs[0] = x;
        for (int i=0 ; i<n_hidden_layers+1 ; i++ ) {
            net_outputs[i+1] = L[i].eval_layer(net_outputs[i]);
        }
        return net_outputs[n_hidden_layers+1];
        
       
    }
    
        
    
};


extern network_handler NH;


template <class C> C act_sigmoid(const C &x) { return 1./(1.+exp(-x));}  // f'(x) = f(x) (1 - f(x))

template <class C> C act_tanh(C x) { return 2.0/(1.+exp(-2.0*x)) - 1.;}  // f'(x) = 1 - f(x)^2  -- a coder en dur si necessaire oour ameliorer la differentiation ?

#define beta_swish 2.

template <class C> C act_relu(C x) {
// approximation by Swish function = xσ(βx), 10 >= beta >= 1
    // a modifier pour rendre correct
    return x*act_sigmoid(beta_swish*x);
}

template <class C> C eval_activation(Activation a, C x)
{
    switch(a){
            case ACT_SIGMOID:
                return act_sigmoid(x);
                break;
            case ACT_TANH:
                return act_tanh(x);
                break;
            case ACT_LINEAR:
                return x;
                break;
            case ACT_RELU:
                return act_relu(x);
                break;
            //    std::cout << "Error in \"compute activation\" : Relu not implemented" << std::endl;
            //    exit(1);
    }
}



// Below for ONNX FILES - relying on Sherlock code

#if ONNX_active
extern computation_graph CG;


// copied from evaluate_graph in computation_graph.cpp  (but external to the class => just redundant so useful only as a basis for the abstract version)
void computation_graph_evaluate_graph( computation_graph CG, map < uint32_t, double > input_node_and_value,
                                      map < uint32_t, double > & output_node_and_value );

// same but for abstract values
// template function has to be in header file...
template <class C> void computation_graph_evaluate_graph_abstract( computation_graph CG, map < uint32_t, C > input_node_and_value,
                                                                           map < uint32_t, C > & output_node_and_value )
{
    // Making sure you received some non empty input map
    assert(!(input_node_and_value.empty()));
    // Making sure that the counts are consistent
    assert(CG.no_of_input_nodes == input_node_and_value.size());
    
    map< uint32_t, C > memoized_table;
    
    // Set the value of the input nodes and make sure those are constant nodes
    
    for(auto input_node : input_node_and_value)
    {
        assert( CG.all_nodes[input_node.first].return_node_type() == const_string );
        //    CG.all_nodes[input_node.first].set_node_val(input_node.second);
        memoized_table[input_node.first] = input_node.second;
        
      //  cout << "Setting node val of node id " << CG.all_nodes[input_node.first].get_node_number() << " as " << memoized_table[input_node.first] << endl;
        
    }
    
    // Creating the space for the outputs of the network
    output_node_and_value.clear();
    for(auto output_node : CG.output_nodes)
    {
        output_node_and_value.insert(make_pair(output_node, C(0.0)));
    }
    
    
    for(auto & output_values : output_node_and_value)
    {
        computation_graph_evaluate_node_abstract(CG, output_values.first, memoized_table, output_values.second);
    }
    
    
}


// copied from evaluate_node from computation_graph.cc as a starting point to implement the abstract version (not udeful in itself)
void computation_graph_evaluate_node(computation_graph & c_graph, uint32_t node_id ,
                                     map< uint32_t , double > & table,
                                     int & available_threads, double & ret_val , int thread_id);

// just a simpler version of computation_graph_evaluate_node without using threads (not used in itslef, just as a test version to derive the abstract version)
void computation_graph_evaluate_node_seq(computation_graph & c_graph, uint32_t node_id ,
                                         map< uint32_t , double > & table, double & ret_val);

// the abstract version
template <class C> void computation_graph_evaluate_node_abstract(computation_graph & c_graph, uint32_t node_id ,
                                                                          map< uint32_t , C > & table, C & ret_val)
{
    
    auto & current_node = c_graph.all_nodes[node_id];
    map< uint32_t , pair< node * , datatype > > backward_connections_, backward_connections;
    current_node.get_backward_connections(backward_connections_);
    backward_connections = backward_connections_;
    
    
    for(auto it = backward_connections.rbegin(); it != backward_connections.rend(); it ++)
    {
        
        auto input_node_ptr = it->second.first;
        
        if(  input_node_ptr->return_node_type() == const_string) // If a constant node then just get the value and store in the table
        {
            
            pair< uint32_t , C > node_and_value = make_pair(it->first, input_node_ptr->return_current_output()) ;
            table.insert(node_and_value);
            continue;
        }
        
        // check if the value is already in the table
        bool value_in_the_table = ( ( table.find(input_node_ptr->get_node_number())  == table.end() ) ? (false) : (true) ) ;
        if( value_in_the_table )
            pair< uint32_t , C > node_and_value = make_pair(it->first, table[input_node_ptr->get_node_number()] ) ;
        else // make  a recursive call to the inputs, get the value and compute
        {
            C buffer;
            computation_graph_evaluate_node_abstract(c_graph, it->first, table, buffer);
            pair< uint32_t, C > node_and_value = make_pair(it->first, buffer);
            table.insert(node_and_value);
        }
    }
    
    
    // Since the input to all the nodes is ready, compute the output now
    
    map< uint32_t, C > inputs_to_the_node;
    for(auto some_connection : backward_connections)
    {
        auto input_node_ptr = some_connection.second.first;
        
        bool value_in_the_table = ( ( table.find(input_node_ptr->get_node_number())  == table.end() ) ? (false) : (true) ) ;
        assert( value_in_the_table );
        pair< uint32_t , C > node_and_value = make_pair(some_connection.first, table[input_node_ptr->get_node_number()] ) ;
        
        inputs_to_the_node.insert(node_and_value);
    }
    
    C result = return_abstract_output(current_node,inputs_to_the_node);
    
    table.insert( make_pair ( node_id , result ) );
   
    ret_val = result;
    
    return;
    
}



template <class C> C return_abstract_output(node current_node, map< uint32_t, C > &inputs_to_the_node)
{
    //  interval result;
    
    type node_type =  current_node.get_node_type();
    
    if (node_type == constant)
    {
        cout << "constant_val=" << current_node.return_current_output() << endl;
        return current_node.return_current_output();
    }
    
    datatype bias;
    current_node.get_bias(bias);
    C argument = bias;
    
    map < uint32_t, pair< node * , datatype > > backward_nodes;
    current_node.get_backward_connections(backward_nodes);
    
    // input_node.second = weight, backward_nodes[input_node.first].second = input
    for(auto input_node : inputs_to_the_node)
    {
        argument += (backward_nodes[input_node.first].second * inputs_to_the_node[input_node.first] ) ;
    }
    
    if(node_type == _tanh_ )
    {
        return act_tanh(argument);
    }
    else if(node_type == _sigmoid_ )
    {
        
        return act_sigmoid(argument);
    }
    else if (node_type == _relu_)
    {
        // last_signature[node_id] = ((argument > 0) ? true : false) ;
        //return (argument > 0) ? argument : 0 ;
        cout << "ReLU not handled ! Exiting....  " << endl;
        exit(0);
    }
    else if(node_type == _none_)
    {
        return argument;
    }
    else
    {
        cout << "Node type not included in the list of evaluation functions ! Exiting....  " << endl;
        exit(0);
    }
    
}


void test_network_sigmoid_cav(computation_graph & CG);
#endif  // endif ONNX_active

#endif
