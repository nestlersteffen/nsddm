
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

    // Copy to Eigen format for fast batched forward pass.
    // arma::mat is column-major (same as Eigen::MatrixXd default), so we can
    // use Eigen::Map to wrap the existing memory and then assign (= copy).
    // W1_e = Eigen::Map<Eigen::MatrixXd>(W1.memptr(), W1.n_rows, W1.n_cols);
    // W2_e = Eigen::Map<Eigen::MatrixXd>(W2.memptr(), W2.n_rows, W2.n_cols);
    // W3_e = Eigen::Map<Eigen::MatrixXd>(W3.memptr(), W3.n_rows, W3.n_cols);
    // W4_e = Eigen::Map<Eigen::MatrixXd>(W4.memptr(), W4.n_rows, W4.n_cols);
    
    // b1_e = Eigen::Map<Eigen::VectorXd>(b1.memptr(), b1.n_elem).transpose();
    // b2_e = Eigen::Map<Eigen::VectorXd>(b2.memptr(), b2.n_elem).transpose();
    // b3_e = Eigen::Map<Eigen::VectorXd>(b3.memptr(), b3.n_elem).transpose();
    // b4_e = Eigen::Map<Eigen::VectorXd>(b4.memptr(), b4.n_elem).transpose();
}

double lan_forward_rcpp( const arma::vec& input ) {
    
    arma::vec x = input;
    x = arma::tanh( W1.t() * x + b1 );
    x = arma::tanh( W2.t() * x + b2 );
    x = arma::tanh( W3.t() * x + b3 );
    x = W4.t() * x + b4;  // linear
    return x(0);
}

// arma::vec lan_forward_batch_rcpp( const arma::mat& X ) {
    
//     const int n = X.n_rows;
//     const int in_dim = X.n_cols;
    
//     // Zero-copy view of the armadillo input as an Eigen matrix
//     Eigen::Map<const Eigen::MatrixXd> X_e( X.memptr(), n, in_dim );
    
//     // Forward pass:
//     // Layer 1: (n x in_dim) * (in_dim x 100) -> (n x 100), add bias, tanh
//     Eigen::MatrixXd A = X_e * W1_e;
//     A.rowwise() += b1_e;
//     A = A.array().tanh().matrix();
    
//     // Layer 2: (n x 100) * (100 x 100) -> (n x 100)
//     Eigen::MatrixXd B = A * W2_e;
//     B.rowwise() += b2_e;
//     B = B.array().tanh().matrix();
    
//     // Layer 3: (n x 100) * (100 x 120) -> (n x 120)
//     A = B * W3_e;
//     A.rowwise() += b3_e;
//     A = A.array().tanh().matrix();
    
//     // Layer 4: (n x 120) * (120 x 1) -> (n x 1), linear (no tanh)
//     B = A * W4_e;
//     B.rowwise() += b4_e;
    
//     // Result: copy column 0 into arma::vec
//     arma::vec result(n);
//     Eigen::Map<Eigen::VectorXd>( result.memptr(), n ) = B.col(0);
//     return result;
// }

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

// arma::vec ddm4_lanll_weights_rcpp_new( const arma::vec& rt, const arma::ivec& x,
//     const double& a, const double& t0, const double& z, const double& v ) {
    
//     // size:
//     int n = rt.size();
//     // build input matrix: 
//     arma::mat X(n, 6);
//     X.col(0).fill(a);
//     X.col(1).fill(v);
//     X.col(2).fill(t0);
//     X.col(3).fill(z);
//     X.col(4) = arma::conv_to<arma::vec>::from(x);  // ivec -> vec
//     X.col(5) = rt;
//     return lan_forward_batch_rcpp(X);
// }

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