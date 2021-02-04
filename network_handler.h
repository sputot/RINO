#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <map>
#include <queue>

#include "aa_aaf.h"

using namespace std;

#ifndef NETWORK_HANDLER_H
#define NETWORK_HANDLER_H


enum Activation {ACT_RELU, ACT_SIGMOID, ACT_TANH, ACT_LINEAR};



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
    
    
    int n_inputs, n_outputs, n_hidden_layers;
    
    network_handler() {};
    
    network_handler(const char* name);
    // the constructor which takes in the information file
    void build_from_file(const char*);
};


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
    
    Layer(network_handler NH, int no_layer): biases(NH.actual_biases[no_layer].size()), weights(NH.actual_biases[no_layer].size(),vector<double>(NH.actual_weights[no_layer][0].size()))
    {
        activation = NH.activations[no_layer];
        nb_inputs = NH.actual_weights[no_layer][0].size();
        nb_outputs = NH.actual_biases[no_layer].size();
        cout << "nb_inputs = " << nb_inputs << " nb_outputs = " << nb_outputs << endl;
        
        for (int i=0 ; i<nb_outputs; i++) {
            for (int j=0 ; j<nb_inputs; j++) {
                weights[i][j] = NH.actual_weights[no_layer][i][j];
            }
        }
        for (int i=0 ; i<nb_outputs; i++)
            biases [i] = NH.actual_biases[no_layer][i];
        cout << "nb_inputs = " << nb_inputs << " nb_outputs = " << nb_outputs << endl;
    }
    
};

#define max_nb_layers 10

extern network_handler NH;
extern vector<Layer> L;


template <class C> C sigmoid(C x) { return 1./(1.+exp(-x));}  // f(x) = f(x) (1 - f(x))

template <class C> C act_tanh(C x) { return 2.0/(1.+exp(-2.0*x)) - 1.;}  // f'(x) = 1 - f(x)^2  -- a coder en dur si necessaire oour ameliorer la differentiation ?


template <class C> C eval_activation(Activation a, C x)
{
    switch(a){
            case ACT_SIGMOID:
                return sigmoid(x);
                break;
            case ACT_TANH:
                return act_tanh(x);
                break;
            case ACT_LINEAR:
                break;
            case ACT_RELU:
                std::cout << "Error in \"compute activation\" : Relu not implemented" << std::endl;
                exit(1);
    }
}

// class EvalLayer {
// public:
    template <class C> vector<C> eval_layer(Layer L, vector<C> x) {
        vector<C> z(L.nb_outputs);
        assert(L.nb_inputs == x.size());
        
        for (int i=0 ; i<L.nb_outputs; i++) {
            z[i] = 0;
            for (int j=0 ; j<L.nb_inputs; j++)
                z[i] += L.weights[i][j]*x[j];
            z[i] = eval_activation(L.activation,z[i]+L.biases[i]);
        }
        return z;
    }








#endif
