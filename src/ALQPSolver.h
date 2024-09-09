#pragma once

#include <Arduino.h>
#include <BasicLinearAlgebra.h>

using namespace BLA;

namespace ALQPS {

template<int nx, int m, int p>
class QPsol{
    public:
        QPsol(Matrix<nx,1> x_, Matrix<m,1> lambda_, Matrix<p,1> mu_, float obj_val_, bool solved, bool pinf, bool dinf, bool degh, bool verbose) : x(x_), lambda(lambda_), mu(mu_), obj_val(obj_val_), _solved(solved), _pinf(pinf), _dinf(dinf), _degh(degh) {
            if (verbose) print_report();
        }; 

        // primal and dual variables 
        Matrix<nx,1> x;
        Matrix<m ,1> lambda;                        
        Matrix<p ,1> mu;
        // objective value and status variables   
        const float obj_val;         
        const bool _solved;
        const bool _pinf;
        const bool _dinf;
        const bool _degh;   

        void print_report(){
            
            Serial.print("status: ");
            if (_solved) Serial.println("solved");
            else if (_pinf) Serial.println("primal infeasible");
            else if (_dinf) Serial.println("dual infeasible");  
            else if (_degh) Serial.println("singular Hessian");  

            Serial.print("optimal objective: ");
            Serial.println(obj_val);

            Serial.print("primal value (solution): ");
            Serial.println(x);

            Serial.print("dual value (equalities): ");
            Serial.println(lambda);

            Serial.print("dual value (inequalities): ");
            Serial.println(mu);
        }          
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
            _max_iter_newton = 25;
            _max_iter_outer  = 25;

            _precision_newton = 1e-6;
            _penalty_initial  = 1.00;
            _penalty_scaling  = 10.0;
            _precision_primal = 1e-6;
        }; 

        // solving the qp
        QPsol<nx,m,p> solve(String mode = ""){
            
            bool verb = false;
            if (mode=="verbose") verb=true;

            Matrix<nx,1> x;
            Matrix<m ,1> lambda;
            Matrix<p ,1> mu;     
            x.Fill(0);
            lambda.Fill(0);
            mu.Fill(0);

            double rho = _penalty_initial;
            double Phi = _penalty_scaling;

            Matrix<1,1> innerprod;
            Matrix<m+p> rp;
            double normrp;

            for (int i = 0; i < _max_iter_outer; i++){
                x = newton_solve(x, lambda, mu, rho, verb);
                if (isnan(x(0))){
                    // condition for rank deficient AL Hessian                     
                    QPsol<nx,m,p> sol(x,lambda,mu,0.0/0.0,false,false,false,true,verb);
                    return sol;
                }                
                if (m+p == 0){
                    // reconsider if stationarity condition is satisfied to determine if x is truly the unconstrained optimum                                            
                    // otherwise re enter newton solver in next iteration 
                    Matrix<nx,1> Nabla_x_L = stationarity(x, lambda, mu); 
                    innerprod = ~Nabla_x_L * Nabla_x_L;
                    double normg = sqrt(innerprod(0));                    
                    if (normg < _precision_newton){
                        QPsol<nx,m,p> sol(x,lambda,mu,objective(x),true,false,false,false,verb);
                        return sol;
                    }
                } 
                                
                // update dual variables lambda, mu
                dual_update(x, lambda, mu, rho); 
                rho = Phi*rho;
                
                rp = primal_residual(x, lambda, mu);
                innerprod = ~rp * rp;
                normrp = sqrt(innerprod(0));
                
                // verify if subset of KKT conditions (primal+dual feasibility) are satisfied 
                if (normrp < _precision_primal && dual_feasibility(mu)){
                    QPsol<nx,m,p> sol(x,lambda,mu,objective(x),true,false,false,false,verb);                    
                    return sol;
                }
            }
            x.Fill(0.0/0.0);
            // primal infeasible
            if (normrp > _precision_primal){
                QPsol<nx,m,p> sol(x,lambda,mu,0.0/0.0,false,true,false,false,verb);
                return sol;                     
            }
            // dual infeasible 
            else if (!dual_feasibility(mu)){
                QPsol<nx,m,p> sol(x,lambda,mu,0.0/0.0,false,false,true,false,verb);
                return sol;
            }
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
        size_t _max_iter_newton;      
        size_t _max_iter_outer;

        double _precision_newton;
        double _penalty_initial;
        double _penalty_scaling;
        double _precision_primal;

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
        // gradient of Lagrangian 
        Matrix<nx,1> stationarity(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu){
            Matrix<nx,1> Nabla_x_L = _Q*x + _q; 
            if (m != 0) {
                Nabla_x_L += ~_A*lambda;
            }
            if (p != 0){
                Nabla_x_L += ~_G*mu;
            }
            return Nabla_x_L;
        }; 
        // measures how much constraint violation
        Matrix<m+p,1>  primal_residual(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu){
            Matrix<m,1> c = c_eq(x);
            Matrix<p,1> h = c_in(x);
            for (int i = 0; i < p; i++){
                h(i) = max(h(i), 0);
            }
            return c && h;
        };
        // dual problem requires non negative duals 
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
            for (int i = 0; i < m; i++){
                lambda(i) = lambda(i)+rho*h(i);
            }  
        };
        // matrix indicating active inequalites 
        Matrix<p,p> active_ineq(Matrix<nx,1> x, Matrix<p,1> mu, float rho){
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
        // gradient of augmented Lagrangian 
        Matrix<nx, 1> algradient(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p> mu, float rho){
            Matrix<nx> Nabla_x_L = stationarity(x, lambda, mu);
            if (m+p == 0){
                return Nabla_x_L;
            }
            else {
                Matrix<nx> Nabla_x_g; 
                Nabla_x_g.Fill(0);
                if(m != 0){
                    Nabla_x_g += (~_A*rho) * c_eq(x);    
                }
                if(p != 0){
                    Matrix<p,p> Ip = active_ineq(x, mu, rho);
                    Nabla_x_g += (~_G*Ip) * c_in(x);
                }
                return Nabla_x_L+Nabla_x_g;
            }
            // Nabla_x_g = (~_A*rho || ~_G*Ip) * primal_residual(x, lambda, mu); 
        };
        // Hessian of augmented Lagrangian 
        Matrix<nx,nx> alhessian(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu, float rho){
            Matrix<nx,nx> Nabla_xx_L = _Q;
            if (m+p == 0){
                return Nabla_xx_L;
            }
            else {
                Matrix<nx,nx> Nabla_xx_g; 
                Nabla_xx_g.Fill(0);
                if(m != 0){
                    Nabla_xx_g += (~_A*rho) * (_A);    
                }
                if(p != 0){
                    Matrix<p,p> Ip = active_ineq(x, mu, rho);
                    Nabla_xx_g += (~_G*Ip) * (_G);
                }
                return Nabla_xx_L+Nabla_xx_g;
            }
        };
        // 'Inner' newton solver 
        Matrix<nx, 1> newton_solve(Matrix<nx,1> x, Matrix<m,1> lambda, Matrix<p,1> mu, float rho, bool verb){
            Matrix<nx,1> x_sol = x;
            for (int i = 0; i < _max_iter_newton; i++){
                
                Matrix<nx,1> g = algradient(x_sol, lambda, mu, rho);
                Matrix<1,1> innerprod = ~g * g;
                double normg = sqrt(innerprod(0));
                
                if (normg < _precision_newton){
                    return x_sol;
                }

                Matrix<nx,nx> H = alhessian(x_sol, lambda, mu, rho);

                // If Hessian of augmented Lagrangian is rank defficient abort procedure 
                if (Determinant(H)==0.0){
                    x_sol.Fill(0.0/0.0);
                    return x_sol;                    
                }
                Matrix<nx> Deltax = -Inverse(H)*g;
                x_sol = x_sol+Deltax;    
            }
            return x_sol;
        };
    };
}; // namespace ALQPS
