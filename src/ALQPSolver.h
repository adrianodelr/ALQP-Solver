#pragma once

#include <Arduino.h>
#include <BasicLinearAlgebra.h>

using namespace BLA;

namespace ALQPS {


template<int n, int m, typename DType>
void printMatrix(const Matrix<n, m, DType>& mat) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            Serial.print(mat(i, j));
            Serial.print(" ");
        }
        Serial.println();
    }
}

// parameter setting for solver 
class QPparams {
public:
    QPparams()
        : max_iter_newton(10),
          max_iter_outer(25),
          precision_newton(1e-6),
          precision_primal(1e-6),
          penalty_initial(1.0),
          penalty_scaling(10.0) {}  // Default values

    size_t max_iter_newton;
    size_t max_iter_outer;
    double precision_newton;
    double precision_primal;
    double penalty_initial;
    double penalty_scaling;
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
        QPsol(Matrix<nx,1,DType> x_, Matrix<m,1,DType> lambda_, Matrix<p,1,DType> mu_, 
              float obj_val_, SolverStatus status, bool verbose) 
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
        QPsol<nx,m,p,DType> solve(String mode = ""){
            
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

            Matrix<1,1,DType> innerprod;
            Matrix<m+p,1,DType> rp;
            double normrp;

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
                    Matrix<nx,1,DType> Nabla_x_L = stationarity(x, lambda, mu); 
                    innerprod = ~Nabla_x_L * Nabla_x_L;
                    double normg = sqrt(innerprod(0));

                    if (normg < _params -> precision_newton){
                        obj_val = objective(x);
                        status.solved = true; 
                        break;
                    }
                } 
                // update dual variables lambda, mu
                dual_update(x, lambda, mu, rho); 
                rho = Phi*rho;

                rp = primal_residual(x, lambda, mu);
                innerprod = ~rp * rp;
                normrp = sqrt(innerprod(0));
                
                // verify if subset of KKT conditions (primal+dual feasibility) are satisfied 
                if (normrp <  _params -> precision_primal && dual_feasibility(mu)){
                    obj_val = objective(x);
                    status.solved = true;
                    break; 
                }
            }

            // check for constraint violation
            if(m+p != 0){
                // primal infeasible
                if (normrp > _params -> precision_primal){
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

        // returns left hand side of the equality constraints   
        Matrix<m,1, DType> c_eq(const Matrix<nx,1,DType>& x){
            return _A*x - _b; 
        };
        // returns left hand side of the inequality constraints   
        Matrix<p,1, DType> c_in(const Matrix<nx,1,DType>& x){
            return _G*x - _h; 
        };
        // computes the objective value 
        float objective(const Matrix<nx,1,DType>& x){
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
        Matrix<nx, 1, DType> algradient(const Matrix<nx,1,DType>& x, const Matrix<m,1,DType>& lambda, const Matrix<p,1,DType>& mu, const DType& rho){
            Matrix<nx,1,DType> Nabla_x_L = stationarity(x, lambda, mu);
            if (m+p == 0){
                return Nabla_x_L;
            }
            else {
                Matrix<nx,1,DType> Nabla_x_g; 
                Nabla_x_g.Fill(static_cast<DType>(0));
                if(m != 0){
                    Nabla_x_g += (~_A*rho) * c_eq(x);    
                }
                if(p != 0){
                    Matrix<p,p,DType> Ip = active_ineq(x, mu, rho);
                    Nabla_x_g += (~_G*Ip) * c_in(x);
                }
                return Nabla_x_L+Nabla_x_g;
            }
            // Nabla_x_g = (~_A*rho || ~_G*Ip) * primal_residual(x, lambda, mu); 
        };
        // Hessian of augmented Lagrangian 
        Matrix<nx,nx,DType> alhessian(const Matrix<nx,1,DType>& x, const Matrix<m,1,DType>& lambda, const Matrix<p,1,DType>& mu, const DType& rho){
            Matrix<nx,nx,DType> Nabla_xx_L = _Q;
            if (m+p == 0){
                return Nabla_xx_L;
            }
            else {
                Matrix<nx,nx,DType> Nabla_xx_g; 
                Nabla_xx_g.Fill(static_cast<DType>(0));
                if(m != 0){
                    Nabla_xx_g += (~_A*rho) * (_A);    
                }
                if(p != 0){
                    Matrix<p,p,DType> Ip = active_ineq(x, mu, rho);
                    Nabla_xx_g += (~_G*Ip) * (_G);
                }
                return Nabla_xx_L+Nabla_xx_g;
            }
        };
        // 'Inner' newton solver 
        Matrix<nx, 1, DType> newton_solve(const Matrix<nx,1,DType>& x, const Matrix<m,1,DType>& lambda, const Matrix<p,1,DType>& mu, const DType& rho, const bool& verb){
            Matrix<nx,1,DType> x_sol = x;
            for (int i = 0; i < _params -> max_iter_newton; i++){
                
                Matrix<nx,1,DType> g = algradient(x_sol, lambda, mu, rho);
                Matrix<1,1,DType> innerprod = ~g * g;
                double normg = sqrt(innerprod(0));
                
                if (normg < _params -> precision_newton){
                    return x_sol;
                }

                Matrix<nx,nx,DType> H = alhessian(x_sol, lambda, mu, rho);

                // If Hessian of augmented Lagrangian is rank defficient abort procedure 
                
                if (Determinant(H)==0.0){
                    x_sol.Fill(0.0/0.0);
                    return x_sol;                    
                }
                Matrix<nx,1,DType> Deltax = -Inverse(H)*g;
                x_sol = x_sol+Deltax;
            }
            return x_sol;
        };
    };
}; // namespace ALQPS
