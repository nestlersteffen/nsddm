
namespace ddmTMB{

	const int MAX = INT_MAX;
	constexpr double PI = acos(-1.0);

	template<class Type>
    Type sumexplog_tmb(vector<Type> x) {
        Type max_x = max(x);
        Type sum_exp = 0;
        for(int i = 0; i < x.size(); i++) {
            sum_exp += exp(x(i) - max_x);
        }
        return max_x + log(sum_exp);
    }

    template<class Type>
    matrix<Type> outer_product( const vector<Type> &x ) {
      int n = x.size();
      matrix<Type> result(n, n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          result(i, j) = x(i) * x(j);
        }
      }
      return result;
    }

    // Rekursives asDouble für beliebig tiefe AD-Typen:
	double asDouble_deep(double x) {
	    return x;
	}

	template<class Type>
	double asDouble_deep(Type x) {
	    return asDouble_deep(CppAD::Value(x));
	}

	// -------------------------
	//       Tuerlinckx 
	// -------------------------

	// Denisty ( CppAD implementation with fixed kmax )

	template<class Type>
	Type ddm4_dstand_tmb( const Type& t, const Type& a, const Type& w, const Type& v, 
		const double& s2, const int& kmax ){  

		// compute multiplicative term:
		Type log_mult = log( s2*PI ) - 2*log(a) + ( -w*a*v - 0.5*v*v*t )/s2;
		
		// we approximate the infinite sum:
		Type tmp_sum = Type(0);

		for ( int k = 1; k < kmax; k++ ){

			//* compute elements of the sum:
			Type pt_1 = sin( PI*k*w );
			Type log_pt_2 = (-0.5*PI*PI*k*k*s2*t)/(a*a);

			//* compute the sum and save:
			Type log_hist = log( Type( k ) ) + log_pt_2;
			tmp_sum +=  pt_1*exp( log_hist );
			
		}
		
		return CppAD::CondExpGt( tmp_sum, Type(0), tmp_sum*exp( log_mult ), Type(1e-29));
		
	}

	// -------------------------
	//      Navarro & Fuss
	// -------------------------

	// Limits

	int get_ks( const double& t, const double& err )
	{
		if (err * 2 * sqrt(2 * PI * t) < 1) {
			double ks  = 2 + sqrt(-2 * t * log(2 * err * sqrt(2 * PI * t)));
			double bc  = sqrt(t) + 1;
			if (ks > MAX || bc > MAX) return MAX;
			return (int) ceil( std::max( ks, bc ) );
		}
		return 2;
	}

	int get_kl( const double& t, const double& err )
	{
		double bc  = 1.0 / ( PI * sqrt(t));
		if (bc > MAX) return MAX;
		if (err * PI * t < 1) {
			double kl = sqrt(-2 * log( PI * t * err) / (PI*PI*t));
			if (kl > MAX) return MAX;
			return (int) ceil( std::max( kl, bc ) );
		}
		return (int) ceil( bc );
	}

	template<class Type>
	Type ddm4_dnavfuss_tmb( const Type& t, const Type& a, const Type& w,
	    const Type& v, const double& err, const int& k_terms = 20 )
	{
	    Type tt   = t / (a * a);
	    Type mult = exp( -w*a*v - Type(0.5)*v*v*t ) / (a*a);

	    // Auswahl auf double-Ebene, außerhalb des AD-Graphen:
	    double tt_d = asDouble_deep(t) / (asDouble_deep(a) * asDouble_deep(a));
	    int ks_nav   = get_ks(tt_d, err);
	    int kl_nav   = get_kl(tt_d, err);
	    bool use_small = (ks_nav < kl_nav);

	    Type dens = Type(0);

	    if( use_small ) {
	        int k_lower = -(k_terms - 1) / 2;
	        int k_upper =   k_terms      / 2;
	        // log-sum-exp Trick wie bisher:
	        Type c = Type(-1e29);
	        for (int k = k_lower; k <= k_upper; k++) {
	            Type wk = w + Type(2*k);
	            Type r  = -wk*wk / (Type(2)*tt);
	            c = CppAD::CondExpGt(r, c, r, c);
	        }
	        Type tmp_sum = Type(0);
	        for (int k = k_lower; k <= k_upper; k++) {
	            Type wk = w + Type(2*k);
	            Type r  = -wk*wk / (Type(2)*tt);
	            tmp_sum += wk * exp(r - c);
	        }
	        dens = exp(c) * tmp_sum / sqrt(Type(2)*PI*tt*tt*tt);
	    } else {
	        Type gamma = Type(-0.5) * PI*PI * tt;
	        for (int k = 1; k <= k_terms; k++) {
	            dens += sin(PI * Type(k) * w) * exp(gamma * Type(k*k)) * Type(k);
	        }
	        dens *= PI;
	    }

	    return CppAD::CondExpGt(dens, Type(0), dens * mult, Type(1e-29));
	}	

	// ------------------------------
	//    7-Parameter model
	// ------------------------------

	template<class Type>
	Type ddm7_integrand_tmb( const Type& tx, const Type& a, const Type& z, const Type& v, const Type& t0, const Type& sv, const int& kmax ) 
	{  

	    // compute true reaction time:
	    Type td = tx - t0;
	    Type td_safe = CppAD::CondExpGt( td, Type(1e-10), td, Type(1e-10) );

	    // compute multiplicative term:
	    Type tdsv2 = 1 + td_safe*sv*sv;
	    Type log_mult1 = log( PI ) - 2*log(a) -0.5*log( tdsv2 );
	    Type log_mult2 = -0.5*( v*v*td_safe + 2*v*z - z*z*sv*sv )/tdsv2;

	    // we approximate the infinite sum:
	    Type tmp_sum = Type(0);

	    for ( int k = 1; k < kmax; k++ ){

	        //* compute elements of the sum:
	        Type pt_1 = sin( (PI*k*z)/a );
	        Type log_pt_2 = (-0.5*PI*PI*k*k*td_safe)/(a*a);

	        //* compute the sum and save:
	        Type log_hist = log( Type( k ) ) + log_pt_2;
	        tmp_sum +=  pt_1*exp( log_hist );
	            
	    }
	    
	    Type out = tmp_sum * exp( log_mult1 ) * exp( log_mult2 );
	    return CppAD::CondExpGt( out, Type(1e-29), out, Type(1e-29) );
	        
	}

	template<class Type>
	Type ddm7_density_dao1_tmb( const Type& tx, const Type& a, const Type& z, const Type& v, const Type& t0, const Type& sv,
		const Type& sz, const Type& st0, const vector<Type> &pts, const vector<Type> &wgh, const int& kmax ) 
	{  

	    // compute true reaction time:
	    Type npts = pts.size();

	    // compute multiplicative term:
	    Type total = Type(0);

	    for ( int i = 0; i < npts; i++ ){
	    	
	        for ( int j = 0; j < npts; j++ ) {

	            //* get current points
	            Type z_ij = z + 0.5*sz*pts(i);
	            Type t0_ij = t0 + 0.5*st0*pts(j);

	            //* check conditions:
	            if (z_ij <= 0 || z_ij >= a) continue;
	            if (t0_ij <= 0) continue;
	            if (tx - t0_ij <= 0) continue;

	            //* compute the integrand:
	            Type f_ij = ddm7_integrand_tmb( tx, a, z_ij, v, t0_ij, sv, kmax );
	            
	            Type valid = CppAD::CondExpGt(z_ij, Type(0),
		            CppAD::CondExpLt(z_ij, a,
		            CppAD::CondExpGt(t0_ij, Type(0),
		            CppAD::CondExpGt(tx - t0_ij, Type(0),
		            f_ij, Type(1e-29)),
		            Type(1e-29)),
		            Type(1e-29)),
		            Type(1e-29));

	            total +=  valid*wgh(i)*wgh(j);

	        }
	            
	    }
	        
	    return 0.25*total;
	        
	}

	template<class Type>
	Type ddm7_density_dao2_tmb( const Type& tx, const Type& a, const Type& z, const Type& v, const Type& t0, const Type& sv,
		const Type& sz, const Type& st0, const vector<Type> &pts, const vector<Type> &wgh, const int& kmax ) 
	{  

	    // compute true reaction time:
	    Type npts = pts.size();

	    //- tau-Grenzen: L_t0 = t0 - st0/2, U_t0 = min(tx, t0 + st0/2)
    	Type L_t0 = t0 - 0.5 * st0;
    	Type U_t0_full = t0 + 0.5 * st0;
    	Type U_t0 = CppAD::CondExpLt( tx, U_t0_full, tx, U_t0_full );

    	//- ...
    	Type width_t = U_t0 - L_t0;
    	Type width_t_safe = CppAD::CondExpGt( width_t, Type(1e-10), width_t, Type(1e-10) );

    	//- helpers:
    	Type half_range_t = 0.5 * width_t_safe;
    	Type mid_t        = 0.5 * (U_t0 + L_t0);

    	// compute multiplicative term:
	    Type total = Type(0);

	    for ( int i = 0; i < npts; i++ ){
	    	
	    	//- current z:
	    	Type z_ij = z + 0.5*sz*pts(i);
	    	Type valid_z = CppAD::CondExpGt( z_ij, Type(0),
                           CppAD::CondExpLt( z_ij, a, Type(1), Type(0) ),
                           Type(0) );

	        for ( int j = 0; j < npts; j++ ) {

	            //* current t0
	            Type t0_ij = half_range_t * pts(j) + mid_t;
				Type valid_t = CppAD::CondExpGt( t0_ij, Type(0), Type(1), Type(0) );
	            
	            //* another check:
	            Type valid_w = CppAD::CondExpGt( width_t, Type(1e-10), Type(1), Type(0) );

	            //* compute the integrand:
	            Type f_ij = ddm7_integrand_tmb( tx, a, z_ij, v, t0_ij, sv, kmax );
	            total += valid_z * valid_t * valid_w * f_ij * wgh(i) * wgh(j);

	        }
	            
	    }
	        
	    return 0.25 * width_t_safe / st0 * total;
	        
	}

	// -----------------------------------
	//    Generic Laplace Approximation
	// -----------------------------------

	template<class Type, class Functor>
	struct laplace_t {
	    Functor f;       
	    vector<Type>& u; 
	    int niter;       
	    Type tolerance;  
	    laplace_t(Functor f_, vector<Type> &u_, int niter_, Type tolerance_ = 1e-6 ) :
	        f(f_), u(u_), niter(niter_), tolerance(tolerance_) {}
	    Type operator()(){
	        vector<Type> u_old = u;
	        int iter = 0;
	        Type max_change = tolerance + 1;
	        while ( iter < niter && max_change > tolerance ) {
	            u_old = u;
	            vector<Type> g = autodiff::gradient(f, u_old);
	            matrix<Type> H = autodiff::hessian(f, u_old);
	            u = u_old - atomic::matinv(H) * g;
	            vector<Type> diff = (u - u_old)*(u - u_old);
	            max_change = diff.sum();
	            iter += 1;
	        }
	        matrix<Type> H = autodiff::hessian(f, u);
	        Type ans = .5 * atomic::logdet(H) + f(u);
	        ans -= .5 * Type(u.size()) * log(2.0 * M_PI);
	        return ans;
	    }
	};

	template<class Type, class Functor>
	Type laplace(Functor f, vector<Type> &u, int niter){
	    laplace_t<Type, Functor> L(f, u, niter);
	    return L();
	}

	// ------------------------------------------------
	//    Generic extended variational approximation
	// ------------------------------------------------

    template<class Type, class Functor>
    struct eva_t {
        Functor f;       // User's implementation of the function
        vector<Type>& u; // Variational mean vector
        matrix<Type>& invSIGMA; // a matrix
        eva_t(Functor f_, vector<Type> &u_, matrix<Type> &invSIGMA_ ) :
            f(f_), u(u_), invSIGMA(invSIGMA_) {}
        Type operator()(){
            // EVA approximation
            matrix<Type> H = autodiff::hessian(f, u);
            matrix<Type> invUistar = (invSIGMA - H);
            // for(int k = 0; k < 4; k++) {
    		// 	invUistar(k,k) += Type(1e-6);
			// }
            Type logdet_val = atomic::logdet(invUistar);
			//if(!std::isfinite(asDouble(logdet_val))) {
    		// 	return Type(-1e6);  // großer Strafwert
			//}
            //Type ans = -.5 * atomic::logdet( invUistar ) + f(u);
            Type ans = logdet_val + f(u);
            return ans;
        }
    };

    template<class Type, class Functor>
    Type eva(Functor f, vector<Type> &u, matrix<Type> &invSIGMA ){
        eva_t<Type, Functor> L(f, u, invSIGMA);
        return L();
    }

}