#ifndef ALQPS_H
#define ALQPS_H

#include <Arduino.h>
#include <BasicLinearAlgebra.h>

using namespace BLA;

namespace ALQPS {

template<int nx, int m, int p>
class QPsol{
    public:
        // updates the qp solution
        void update(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu, float objective, bool solved, bool pinf, bool dinf){
            _x=x;
            _lambda=lambda;                        
            _mu=mu;
            _obj_val=objective;
            _solved=solved;
            _pinf=pinf;
            _dinf=dinf;            
        }; 

        QPsol(){
            // initialize primal, dual variables as zero and objective as NaN 
            Matrix<nx,1> x;
            Matrix<m ,1> lambda;                        
            Matrix<p ,1> mu;   
            x.Fill(0.);
            lambda.Fill(0.);            
            mu.Fill(0.);

            update(x,lambda,mu,0.0,false,false,false);    
        }; 
        // primal and dual variables
        Matrix<nx,1>  _x;
        Matrix<m ,1>  _lambda;                        
        Matrix<p ,1>  _mu;
        // objective value and status variables   
        float _obj_val;         
        bool _solved;
        bool _pinf;
        bool _dinf;            
};

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
            max_iter_newton = 10;
            max_iter_outer  = 10;

            precision_newton = 1e-8;
            penalty_initial  = 1.00;
            penalty_scaling  = 10.0;
            precision_primal = 1e-8;
        }; 
        // solving the qp
        QPsol<nx,m,p> solve(){
            Matrix<nx,1> x;
            Matrix<m ,1> lambda;
            Matrix<p ,1> mu;     
            x.Fill(0);
            lambda.Fill(0);
            mu.Fill(0);

            double rho = penalty_initial;
            double Phi = penalty_scaling;

            Matrix<m+p> rp;
            Matrix<1,1> innerprod;
            double normrp;

            for (int i = 0; i < max_iter_outer; i++){
                x = newton_solve(x, lambda, mu, rho);
                // update dual variables lambda, mu
                dual_update(x, lambda, mu, rho); 
                rho = Phi*rho;
                rp = primal_residual(x, lambda, mu);
                innerprod = ~rp * rp;
                normrp = sqrt(innerprod(0));
                // verify conditions if subset of KKT conditions are satisfied 
                if (normrp < precision_primal && dual_feasibility(mu)){
                
                    sol.update(x,lambda,mu,objective(x),true,false,false);                    
                    // Serial.println("Found optimal solution");
                    return sol;
                }
            }
            if (normrp > precision_primal){
                sol.update(x,lambda,mu,objective(x),false,true,false);                    
                Serial.println("Problem is primal infeasible");
            }
            if (!dual_feasibility(mu)){
                sol.update(x,lambda,mu,objective(x),false,false,true);                    
                Serial.println("Problem is dual infeasible");
            }
            return sol;
        };
        QPsol<nx,m,p> sol;
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
        // computes the objective value 
        float objective(Matrix<nx,1> x){
            return 0.5*(~x*_Q*x)(0) + (~_q*x)(0);             
        };


        Matrix<nx,1> stationarity(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu){
            // Since basic linear algebra seems not to perform addition well with empty matrices use below approach 
            // return _Q*x + _q + ~_A*lambda +  ~_G*mu; 
            Matrix<nx,1> pr = _Q*x + _q; 
            if (m != 0) {
                pr += ~_A*lambda;
            }
            if (p != 0){
                pr += ~_G*mu;
            }
            return pr;
        }; 
        
        Matrix<m+p,1>  primal_residual(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu){
            Matrix<m,1> c = c_eq(x);
            Matrix<p,1> h = c_in(x);
            for (int i = 0; i < p; i++){
                h(i) = max(h(i), 0);
            }
            return c && h;
        };

        bool dual_feasibility(Matrix<p,1> mu){
            for (int i = 0; i < p; i++){
                if(mu(i) < 0.0) return false;
            }
            return true;
        }

        void dual_update(Matrix<nx,1> x, Matrix<m,1> &lambda, Matrix<p> &mu, float rho){
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
        
        Matrix<p,p> active_set(Matrix<nx,1> x, Matrix<p,1> mu, float rho){
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

        Matrix<nx, 1> algradient(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p> mu, float rho){
            Matrix<nx> Nabla_x_L = stationarity(x, lambda, mu);
            Matrix<p,p> Ip = active_set(x, mu, rho);
            Matrix<nx> Nabla_x_g; 
            Nabla_x_g = (~_A*rho || ~_G*Ip) * primal_residual(x, lambda, mu);
            return Nabla_x_L+Nabla_x_g;
        };

        Matrix<nx,nx> alhessian(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu, float rho){
            Matrix<nx,nx> Nabla_xx_L = _Q;
            Matrix<p,p> Ip = active_set(x, mu, rho);
            Matrix<nx,nx> Nabla_xx_g; 
            if (m != 0){
                Nabla_xx_g = (~_A*rho || ~_G*Ip) * (_A && _G);
            }
            else {
                Nabla_xx_g = (~_G*Ip) * (_G);
            }    
            return Nabla_xx_L+Nabla_xx_g;
        };

        Matrix<nx, 1> newton_solve(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu, float rho){
            Matrix<nx,1> x_sol = x;
            for (int i = 0; i < max_iter_newton; i++){

                Matrix<nx,1> g = algradient(x_sol, lambda, mu, rho);

                Matrix<1,1> innerprod = ~g * g;
                double normg = sqrt(innerprod(0));

                if (normg < precision_newton){
                    return x_sol;
                }

                Matrix<nx,nx> H = alhessian(x_sol, lambda, mu, rho);
                // build in something for feasibility checking 
                Matrix<nx> Deltax = -Inverse(H)*g;
                x_sol = x_sol+Deltax;    
            }
            return x_sol;
        };
    };
}; // namespace ALQPS

#endif // ALQPS_H