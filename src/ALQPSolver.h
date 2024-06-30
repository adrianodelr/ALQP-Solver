#ifndef ALQPS_H
#define ALQPS_H

#include <Arduino.h>
#include <BasicLinearAlgebra.h>

using namespace BLA;

namespace ALQPS {

template<int nx, int m, int p>

class QP{
    public:
        QP(int n_x, int n_eq, int n_ineq); // Constructor for empty qp

        // updates the qp
        void update(Matrix<nx,nx> Q, Matrix<nx,1> q, Matrix<m,nx> A, Matrix<m,1> b, Matrix<p,nx> G, Matrix<p,1> h); 
        
        // solving the qp
        Matrix<nx,1> solve();

    private: 
        int _n_x;                   // dimension of the search space 
        int _n_eq;                  // number of equalities 
        int _n_ineq;                // number of inequalities 
        Matrix<nx,nx> _Q;           // quadratic coefficient matrix  
        Matrix<nx,1>  _q;           // linear coefficient vector 
        Matrix<m,nx>  _A;           // equality constraint matrix 
        Matrix<m,1>   _b;           // equality constraint vector 
        Matrix<p,nx>  _G;           // inequality constraint matrix 
        Matrix<p,1>   _h,           // inequality constraint vector 

        // returns left hand side of the equality constraints i.e. A*x-b  
        Matrix<m,1> c_eq(Matrix<nx,1> x);
        // returns left hand side of the inequality constraints  i.e. G*x-h 
        Matrix<p,1> c_in(Matrix<nx,1> x);
        Matrix<nx, 1> primal_residual(Matrix<nx> x, Matrix<m> lambda, Matrix<p> mu); 
        Matrix<m+p,1> dual_residual(Matrix<nx> x, Matrix<m> lambda, Matrix<p> mu); 
        Matrix<nx, 1> newton_solve(Matrix<nx> x, Matrix<m> lambda, Matrix<p> mu, int rho);
        Matrix<nx, 1> ALgradient(Matrix<nx> x, Matrix<m> lambda, Matrix<p> mu, int rho)
        Matrix<nx,nx> ALhessian(Matrix<nx> x, Matrix<m> lambda, Matrix<p> mu, int rho);
        void dual_update(Matrix<nx> x, Matrix<m> &lambda, Matrix<p> &mu, int rho);
        Matrix<p,p> coinactfilt(Matrix<nx> x, Matrix<p> mu, int rho);
};

} // namespace ALQPS

#endif // ALQPS_H