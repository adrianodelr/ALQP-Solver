#pragma once

#include <Arduino.h>
#include <BasicLinearAlgebra.h>

using namespace BLA;

namespace ALQPS {

// parameter setting for solver 
class QPparams {
public:
    QPparams()
        : max_iter_newton(10),
          max_iter_outer(25),
          max_iter_backtrack(10),          
          precision_newton(1e-5),
          precision_primal(1e-4),
          penalty_initial(1.0),
          penalty_scaling(10.0), 
          backtrack_beta(0.8)           
          {}  // Default values

    size_t max_iter_newton;
    size_t max_iter_outer;
    size_t max_iter_backtrack;
    double precision_newton;
    double precision_primal;
    double penalty_initial;
    double penalty_scaling;
    double backtrack_beta; 
};

// status indicating solver success or infeasibility
struct SolverStatus{
    bool solved = false; 
    bool pinf = false; 
    bool dinf = false;
    bool degh = false;
};

// solution class 
template<int nx, int m, int p, typename DType = float>
class QPsol{
    public:
        QPsol(const Matrix<nx,1,DType>& x_, const Matrix<m,1,DType>& lambda_, const Matrix<p,1,DType>& mu_, 
              float obj_val_, const SolverStatus& status, bool verbose) 
            : x(x_), 
              lambda(lambda_), 
              mu(mu_), 
              obj_val(obj_val_), 
              _solved(status.solved), 
              _pinf(status.pinf), 
              _dinf(status.dinf), 
              _degh(status.degh) {
            if (verbose) print_report();
        }; 

        // primal and dual variables 
        Matrix<nx,1,DType> x;
        Matrix<m ,1,DType> lambda;                        
        Matrix<p ,1,DType> mu;

        // objective value and status variables   
        const DType obj_val;         
        const bool _solved;
        const bool _pinf;
        const bool _dinf;
        const bool _degh;   

        void print_report(){
            
            Serial.print(F("status: "));
            if (_solved) Serial.println(F("solved"));
            else if (_pinf) Serial.println(F("primal infeasible"));
            else if (_dinf) Serial.println(F("dual infeasible"));  
            else if (_degh) Serial.println(F("singular Hessian"));  

            Serial.print(F("optimal objective: "));
            Serial.println(obj_val);
            Serial.print(F("primal value (solution): "));
            Serial.println(x);
            Serial.print(F("dual value (equalities): "));
            Serial.println(lambda);
            Serial.print(F("dual value (inequalities): "));
            Serial.println(mu);
        }          
};

template<int nx, int m, int p, typename DType=float>
class QP{
    public:
        // updates the qp
        void update(const Matrix<nx,nx, DType>& Q, const Matrix<nx,1,DType>& q, 
                    const Matrix<m, nx, DType>& A, const Matrix<m, 1,DType>& b, 
                    const Matrix<p, nx, DType>& G, const Matrix<p, 1,DType>& h){
            _Q = Q;
            _q = q;
            _A = A;
            _b = b;                        
            _G = G;
            _h = h;            
        }; 
        // Default constructor for empty qp
        QP() 
            : _params(new QPparams()), _alloc(true) {            
            init_QPmat();
        }; 

        QP(const QPparams& parameters) 
            : _params(&parameters), _alloc(false){
            init_QPmat();
        }

        // destructor in case of taking default parameters 
        ~QP(){
            if(_alloc) delete _params;
        }

        // solving the qp
        QPsol<nx,m,p,DType> solve(const String& mode = ""){
            
            bool verb = false;
            if (mode=="verbose") verb=true;

            Matrix<nx,1,DType> x;
            Matrix<m ,1,DType> lambda;
            Matrix<p ,1,DType> mu;     

            x.Fill(static_cast<DType>(0.0));
            lambda.Fill(static_cast<DType>(0.0));
            mu.Fill(static_cast<DType>(0.0));

            DType rho = _params -> penalty_initial;
            DType Phi = _params -> penalty_scaling;

            Matrix<m+p,1,DType> pr;
            Matrix<nx,1,DType> Nabla_x_L;
            double prnorm;

            DType obj_val; 
            SolverStatus status;

            for (int i = 0; i < _params -> max_iter_outer; i++){
                x = newton_solve(x, lambda, mu, rho, verb);
                if (isnan(x(0))){
                    // condition for rank deficient AL Hessian
                    status.degh = true;
                    obj_val = 0.0/0.0;
                    break;
                }                
                if (m+p == 0){
                    // reconsider if stationarity condition is satisfied to determine if x is truly the unconstrained optimum                                            
                    // otherwise re enter newton solver in next iteration
                    Nabla_x_L = stationarity(x, lambda, mu); 
                    double gnorm = residual_norm(Nabla_x_L);                    

                    if (gnorm < _params -> precision_newton){
                        obj_val = objective(x);
                        status.solved = true; 
                        break;
                    }
                } 
                else {
                    // update dual variables lambda, mu
                    dual_update(x, lambda, mu, rho); 
                    rho = Phi*rho;

                    pr = primal_residual(x, lambda, mu);
                    prnorm = sqrt((~pr * pr)(0));
                    
                    // verify if subset of KKT conditions (primal+dual feasibility) are satisfied 
                    if (prnorm <  _params -> precision_primal && dual_feasibility(mu)){
                        obj_val = objective(x);
                        status.solved = true;
                        break; 
                    }
                }
            }

            // check for constraint violation
            if(m+p != 0){
                // primal infeasible
                if (prnorm > _params -> precision_primal){
                    status.pinf = true;
                    obj_val = 0.0/0.0; 
                    x.Fill(0.0/0.0);
                }

                // dual infeasible 
                if (!dual_feasibility(mu)){
                    status.dinf = true;
                    obj_val = 0.0/0.0;
                    x.Fill(0.0/0.0);
                }
            }

            QPsol<nx,m,p,DType> sol(x,lambda,mu,obj_val,status,verb);
            return sol;            
        };

    private: 
        // QP matrices 
        Matrix<nx,nx, DType> _Q;           // quadratic coefficient matrix  
        Matrix<nx, 1, DType> _q;           // linear coefficient vector 
        Matrix<m, nx, DType> _A;           // equality constraint matrix 
        Matrix<m,  1, DType> _b;           // equality constraint vector 
        Matrix<p, nx, DType> _G;           // inequality constraint matrix 
        Matrix<p,  1, DType> _h;           // inequality constraint vector 
        
        // Solver settings
        const QPparams* _params;
        bool _alloc; 
        
        // initializes all matrices of given dimensions with zeros 
        void init_QPmat(){
            Matrix<nx,nx, DType> Q;
            Matrix<nx, 1, DType> q;
            Matrix<m, nx, DType> A;
            Matrix<m,  1, DType> b;                        
            Matrix<p, nx, DType> G;
            Matrix<p,  1, DType> h; 

            Q.Fill(static_cast<DType>(0.0));
            q.Fill(static_cast<DType>(0.0));
            A.Fill(static_cast<DType>(0.0));
            b.Fill(static_cast<DType>(0.0));            
            G.Fill(static_cast<DType>(0.0));           
            h.Fill(static_cast<DType>(0.0));
            update(Q, q, A, b, G, h);
        }

        // returns 'left hand side' of the equality constraints   
        inline Matrix<m,1, DType> c_eq(const Matrix<nx,1,DType>& x){
            return _A*x - _b; 
        };
        // returns 'left hand side' of the inequality constraints   
        inline Matrix<p,1, DType> c_in(const Matrix<nx,1,DType>& x){
            return _G*x - _h; 
        };
        // returns norm of a residual vector   
        inline DType residual_norm(const Matrix<nx,1,DType>& g){
            return sqrt((~g * g)(0));
        }

        // computes the objective value 
        DType objective(const Matrix<nx,1,DType>& x){
            return 0.5*(~x*_Q*x)(0) + (~_q*x)(0);             
        };
        // gradient of Lagrangian 
        Matrix<nx,1,DType> stationarity(const Matrix<nx,1,DType>& x, const Matrix<m,1,DType>& lambda, const Matrix<p,1,DType>& mu){
            Matrix<nx,1,DType> Nabla_x_L = _Q*x + _q; 
            if (m != 0) {
                Nabla_x_L += ~_A*lambda;
            }
            if (p != 0){
                Nabla_x_L += ~_G*mu;
            }
            return Nabla_x_L;
        }; 

        // measures how much constraint violation
        Matrix<m+p,1,DType>  primal_residual(const Matrix<nx,1,DType>& x, const Matrix<m,1,DType>& lambda, const Matrix<p,1,DType>& mu){
            Matrix<m,1,DType> c = c_eq(x);
            Matrix<p,1,DType> h = c_in(x);
            for (int i = 0; i < p; i++){
                h(i) = max(h(i,0), static_cast<DType>(0));
            }
            return c && h;
        };

        // dual problem requires non negative duals 
        bool dual_feasibility(const Matrix<p,1,DType>& mu){
            for (int i = 0; i < p; i++){
                if(mu(i,0) < 0.0) return false;
            }
            return true;  
        }

        void dual_update(const Matrix<nx,1,DType>& x, Matrix<m,1,DType> &lambda, Matrix<p,1,DType> &mu, const DType& rho){
            Matrix<p,1,DType> c = c_in(x);
            Matrix<m,1,DType> h = c_eq(x);
            for (int i = 0; i < p; i++){
                mu(i,0) = max(static_cast<DType>(0.0), mu(i)+rho*c(i));
            }
            for (int i = 0; i < m; i++){
                lambda(i,0) = lambda(i,0)+rho*h(i,0);
            }  
        };

        // matrix indicating active inequalites 
        Matrix<p,p,DType> active_ineq(const Matrix<nx,1,DType>& x, const Matrix<p,1,DType>& mu, const float& rho){
            Matrix<p,p,DType> Ip;
            Ip.Fill(static_cast<DType>(0));
            Matrix<p,1,DType> h = c_in(x);
            for (int i = 0; i < p; i++){
                if (h(i,0) < 0 && mu(i,0) == 0){
                    Ip(i,i) = static_cast<DType>(0);
                }
                else {
                    Ip(i,i) = rho;
                }
            }
            return Ip;           
        };

        // gradient of augmented Lagrangian 
        void algradient(Matrix<nx,1,DType>& g, const Matrix<nx,1,DType>& x, const Matrix<m,1,DType>& lambda, const Matrix<p,1,DType>& mu, const DType& rho){
            g = stationarity(x, lambda, mu);
            if (m+p != 0){
                if(m != 0){
                    g += (~_A*rho) * c_eq(x);    
                }
                if(p != 0){
                    Matrix<p,p,DType> Ip = active_ineq(x, mu, rho);
                    g += (~_G*Ip) * c_in(x);
                }
            }
            // Nabla_x_g = (~_A*rho || ~_G*Ip) * primal_residual(x, lambda, mu); 
        };

        // Hessian of augmented Lagrangian 
        void alhessian(Matrix<nx,nx,DType>& H, const Matrix<nx,1,DType>& x, const Matrix<m,1,DType>& lambda, const Matrix<p,1,DType>& mu, const DType& rho ){
            H = _Q;
            if(m != 0){
                H += (~_A*rho) * (_A);    
            }
            if(p != 0){
                Matrix<p,p,DType> Ip = active_ineq(x, mu, rho);
                H += (~_G*Ip) * (_G);
            }
        };        

        // 'Inner' newton solver 
        Matrix<nx, 1, DType> newton_solve(const Matrix<nx,1,DType>& x, const Matrix<m,1,DType>& lambda, const Matrix<p,1,DType>& mu, const DType& rho, const bool& verb){
            
            Matrix<nx,1,DType> x_sol = x;            
            Matrix<nx,1,DType> Deltax; 

            // gradient and hessian of augmented Lagrangian 
            Matrix<nx,1,DType> g; 
            Matrix<nx,nx,DType> H; 

            for (int i = 0; i < _params -> max_iter_newton; i++){
                
                algradient(g, x_sol, lambda, mu, rho);
                DType gnorm = residual_norm(g);

                if (gnorm < _params -> precision_newton){
                    return x_sol;
                }

                alhessian(H, x_sol, lambda, mu, rho);

                // If Hessian of augmented Lagrangian is rank defficient abort procedure 
                if (Determinant(H)==0.0){
                    x_sol.Fill(0.0/0.0);
                    return x_sol;                    
                }
                Deltax = -Inverse(H)*g;
                
                // simple back tracking line search on the AL gradient residual 
                DType alpha = 1.0;                  // scaling parameter 
                Matrix<nx,1,DType> x_backtrack;     // solution after taking reduced newton step

                for (int i = 0; i < _params -> max_iter_backtrack; i++){
                    // solution with reduced stepsize
                    x_backtrack = x_sol+alpha*Deltax;
                    
                    // compute residual with reduced stepsize 
                    algradient(g, x_backtrack, lambda, mu, rho);
                    DType gnorm_red = residual_norm(g);
                        
                    // if residual before newton step is higher than after, reduce the scaling parameter alpha   
                    if (gnorm_red > gnorm)
                        alpha *= _params -> backtrack_beta; 
                    else 
                        break;
                }
                x_sol = x_backtrack;
            }
            return x_sol;
        };
    };
}; // namespace ALQPS
