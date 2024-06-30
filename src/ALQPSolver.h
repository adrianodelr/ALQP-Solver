#ifndef ALQPS_H
#define ALQPS_H

#include <Arduino.h>
#include <BasicLinearAlgebra.h>

using namespace BLA;

namespace ALQPS {

template<int nx, int m, int p>
class QP{
    public:
        // updates the qp
        void update(Matrix<nx,nx> Q, Matrix<nx,1> q, Matrix<m,nx> A, Matrix<m,1> b, Matrix<p,nx> G, Matrix<p,1> h){
            _Q=Q;
            _q=q;
            _A=A;
            _b=b;                        
            _G=G;
            _h=h;            
        }; 
        // Default constructor for empty qp
        QP(){
            Matrix<nx,nx> Q;
            Matrix<nx,1>  q;
            Matrix<m ,nx> A;
            Matrix<m ,1>  b;                        
            Matrix<p ,nx> G;
            Matrix<p ,1>  h;   
            Q.Fill(0.);
            q.Fill(0.);
            A.Fill(0.);
            b.Fill(0.);            
            G.Fill(0.);           
            h.Fill(0.);
            update(Q, q, A, b, G, h);

            // fixed solver settings 
            max_iter_newton = 20;
            max_iter_outer  = 50;

            precision_newton = 1e-5;
            penalty_initial  = 10.0;
            penalty_scaling  = 10.0;
            precision_primal = 1e-6;
        }; 
        // solving the qp
        Matrix<nx,1> solve(){
            Matrix<nx> x;
            Matrix<m> lambda;
            Matrix<p> mu;     
            x.Fill(0);
            lambda.Fill(0);
            mu.Fill(0);

            double rho = penalty_initial;
            double Phi = penalty_scaling;

            for (int i = 0; i < max_iter_outer; i++){
                x = newton_solve(x, lambda, mu, rho);
                // update parameters lambda,mu
                dual_update(x, lambda, mu, rho); 
                rho = Phi*rho;
                Matrix<m+p> rd = dual_residual(x, lambda, mu);
                Matrix<1,1> innerprod = ~rd * rd;
                double normdr = sqrt(innerprod(0));
                if (normdr < precision_primal){
                    return x;
                }
            }
            return x;
        };

    private: 
        // QP matrices 
        Matrix<nx,nx> _Q;           // quadratic coefficient matrix  
        Matrix<nx,1>  _q;           // linear coefficient vector 
        Matrix<m ,nx> _A;           // equality constraint matrix 
        Matrix<m ,1>  _b;           // equality constraint vector 
        Matrix<p ,nx> _G;           // inequality constraint matrix 
        Matrix<p ,1>  _h;           // inequality constraint vector 
        
        // Solver settings
        size_t max_iter_newton;
        size_t max_iter_outer;

        double precision_newton;
        double penalty_initial;
        double penalty_scaling;
        double precision_primal;

        // returns left hand side of the equality constraints   
        Matrix<m,1> c_eq(Matrix<nx,1> x){
            return _A*x - _b;
        };
        // returns left hand side of the inequality constraints   
        Matrix<p,1> c_in(Matrix<nx,1> x){
            return _G*x - _h;
        };

        Matrix<nx,1> primal_residual(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu){
            if (m != 0){
                return x + _q + ~_A*lambda +  ~_G*mu;
            }
            else{
                return _Q*x + _q + ~_G*mu;
            }            
        }; 

        Matrix<m+p,1>  dual_residual(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu){
            Matrix<m,1> c = c_eq(x);
            Matrix<p,1> h = c_in(x);
            for (int i = 0; i < p; i++){
                h(i) = max(h(i), 0);
            }
            if (m != 0){
                return c && h;
            }
            else{
                return h;
            }            
        };

        void dual_update(Matrix<nx,1> x, Matrix<m,1> &lambda, Matrix<p> &mu, int rho){
            Matrix<p,1> c = c_in(x);
            Matrix<m,1> h = c_eq(x);
            for (int i = 0; i < p; i++){
                mu(i) = max(0, mu(i)+rho*c(i));
            }
            if (m != 0){
                for (int i = 0; i < p; i++){
                    lambda(i) = lambda(i)+rho*h(i);
                }
            }            
        };
        
        Matrix<p,p> coinactfilt(Matrix<nx,1> x, Matrix<p,1> mu, int rho){
            Matrix<p,p> Ip;
            Ip.Fill(0);
            Matrix<p,1> h = c_in(x);
            for (int i = 0; i < p; i++){
                if (h(i) < 0 && mu(i) == 0){
                    Ip(i,i) = 0;
                }
                else {
                    Ip(i,i) = rho;
                }
            }
            return Ip;            
        };

        Matrix<nx, 1> ALgradient(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p> mu, int rho){
            Matrix<nx> Nabla_x_L = primal_residual(x, lambda, mu);
            Matrix<p,p> Ip = coinactfilt(x, mu, rho);
            Matrix<nx> Nabla_x_g; 
            if (m != 0){
                Nabla_x_g = ~_A || ~_G*Ip * dual_residual(x, lambda, mu);
            }
            else {
                Nabla_x_g = ~_G*Ip * dual_residual(x, lambda, mu);
            }    
            return Nabla_x_L+Nabla_x_g;
        };

        Matrix<nx,nx> ALhessian(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu, int rho){
            Matrix<nx,nx> Nabla_xx_L = _Q;
            Matrix<p,p> Ip = coinactfilt(x, mu, rho);
            Matrix<nx,nx> Nabla_xx_g; 
            if (m != 0){
                Nabla_xx_g = (~_A || ~_G*Ip) * (_A && _G);
            }
            else {
                Nabla_xx_g = (~_G*Ip) * (_G);
            }    
            return Nabla_xx_L+Nabla_xx_g;
        };

        Matrix<nx, 1> newton_solve(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu, int rho){
            Matrix<nx,1> x_sol = x;
            for (int i = 0; i < max_iter_newton; i++){

                Matrix<nx,1> g = ALgradient(x_sol, lambda, mu, rho);
                Matrix<1,1> innerprod = ~g * g;
                double normg = sqrt(innerprod(0));

                if (normg < precision_newton){
                    return x_sol;
                }

                Matrix<nx,nx> H = ALhessian(x_sol, lambda, mu, rho);
                Matrix<nx,nx> Hinv;
                Hinv=Invert(H);
                Matrix<nx> Deltax = -Hinv*g;
                x_sol = x_sol+Deltax;    
            }
            return x_sol;
        };
    };
}; // namespace ALQPS

#endif // ALQPS_H