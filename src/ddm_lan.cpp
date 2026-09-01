
#include "ddm_lan.h" 

using namespace Rcpp;
using namespace arma;

arma::vec ddm4_lanll_dnn_rcpp( const arma::vec& rt, const arma::ivec& x, const double& a, 
     const double& t0, const double& z, const double& v, const int& K, SEXP dnn ) {
     Rcpp::Environment myEnv = Rcpp::Environment::global_env();
     Rcpp::Function myFun = myEnv["lan_loglik"];
     Rcpp::NumericVector result = myFun( 
         Rcpp::Named("dnn", dnn),
         Rcpp::Named("a",   a),
         Rcpp::Named("v",   v),
         Rcpp::Named("t0",  t0),
         Rcpp::Named("z",   z),
         Rcpp::Named("rt",  rt),
         Rcpp::Named("xs",  x),
         Rcpp::Named("K",  K)
     );
     return Rcpp::as<arma::vec>( result );
}

void lan_load_weights_rcpp( const std::string& path, const std::string& ddm ) {
    
    // load und reshape:
    W1.load(path + "W1.bin", arma::raw_binary); 
    if ( ddm == "seven" ) {
        W1.reshape(9, 100);
    } else {  // "four"
        W1.reshape(6, 100);
    }
    W2.load(path + "W2.bin", arma::raw_binary); W2.reshape(100, 100);
    W3.load(path + "W3.bin", arma::raw_binary); W3.reshape(100, 120);
    W4.load(path + "W4.bin", arma::raw_binary); W4.reshape(120,   1);
    b1.load(path + "b1.bin", arma::raw_binary);
    b2.load(path + "b2.bin", arma::raw_binary);
    b3.load(path + "b3.bin", arma::raw_binary);
    b4.load(path + "b4.bin", arma::raw_binary);
    Rcpp::Rcout << "Weights loaded from " << path << std::endl;

}

double lan_forward_rcpp( const arma::vec& input ) {
    
    arma::vec x = input;
    x = arma::tanh( W1.t() * x + b1 );
    x = arma::tanh( W2.t() * x + b2 );
    x = arma::tanh( W3.t() * x + b3 );
    x = W4.t() * x + b4;  // linear
    return x(0);
}

arma::vec lan_forward_backward_rcpp( const arma::vec& input, const int& idx ) {
    
    //- list to save activations:
    std::vector<arma::vec> activations(5);
    
    //- make forward pass:
    activations[0] = input;
    arma::vec x = input; 
    x = arma::tanh( W1.t() * x + b1 );  
    activations[1] = x;
    x = arma::tanh( W2.t() * x + b2 );  
    activations[2] = x;
    x = arma::tanh( W3.t() * x + b3 );
    activations[3] = x;
    x = W4.t() * x + b4;
    activations[4] = x;
    double output = x(0);
    
    //- make backward Pass
    arma::mat d = arma::ones(1, 1);  // = matrix(1.0) in R
    d = (W4 * d) % (1.0 - arma::square( activations[3] ));
    d = (W3 * d) % (1.0 - arma::square( activations[2] ));
    d = (W2 * d) % (1.0 - arma::square( activations[1] ));
    
    //- gradient for the input:
    arma::vec grad_input = W1 * d;
    
    //- all for the output:
    arma::vec result(1 + idx);
    result(0) = output;
    result.subvec(1, idx) = grad_input.head(idx);
    return result;
}

arma::vec ddm4_lanll_weights_rcpp( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v ) {
    
    int n = rt.size();
    arma::vec ll(n);
    arma::vec input(6);
    
    for ( int i = 0; i < n; i++ ) {
        input(0) = a;
        input(1) = v;
        input(2) = t0;
        input(3) = z;
        input(4) = x(i);
        input(5) = rt(i);
        ll(i) = lan_forward_rcpp( input );
    }
    
    return ll;
}

arma::vec ddm4_lanll_weights_batch_rcpp( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v ) {
    
    int n = rt.size();

    // Rcpp::Rcout << "Hier bin ich:" << n << std::endl;
    
    // Input-Matrix bauen: (n x 6)
    arma::mat X(n, 6);
    X.col(0).fill(a);
    X.col(1).fill(v);
    X.col(2).fill(t0);
    X.col(3).fill(z);
    X.col(4) = arma::conv_to<arma::vec>::from(x);
    X.col(5) = rt;
    
    // Forward-Pass als Batch (eine Matrixmultiplikation pro Layer statt n einzelne):
    arma::mat H1 = X * W1;
    H1.each_row() += b1.t();
    // H1 = arma::tanh(H1);
    
    arma::mat H2 = H1 * W2;
    H2.each_row() += b2.t();
    H2 = arma::tanh(H2);
    
    arma::mat H3 = H2 * W3;
    H3.each_row() += b3.t();
    H3 = arma::tanh(H3);
    
    arma::mat out = H3 * W4;
    out.each_row() += b4.t();   // linear, keine Aktivierung
    
    return out.col(0);   // (n x 1) -> Vektor der Länge n
}

arma::vec ddm7_lanll_weights_rcpp( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v,
    const double& sv, const double& sz, const double& st0 ) {
    
    int n = rt.size();
    arma::vec ll(n);
    arma::vec input(9);  
    
    for ( int i = 0; i < n; i++ ) {
        input(0) = a;
        input(1) = v;
        input(2) = t0;
        input(3) = z;
        input(4) = sv;   
        input(5) = sz;   
        input(6) = st0;  
        input(7) = x(i);
        input(8) = rt(i);
        ll(i) = lan_forward_rcpp( input );
    }
    
    return ll;
}