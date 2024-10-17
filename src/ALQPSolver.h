#pragma once

#include <Arduino.h>
#include <BasicLinearAlgebra.h>

extern unsigned int __heap_start;
extern void *__brkval;

using namespace BLA;

namespace ALQPS {

int freeMemory() {
  int free_memory;
  if ((int)__brkval == 0) {
    free_memory = ((int)&free_memory) - ((int)&__heap_start);
  } else {
    free_memory = ((int)&free_memory) - ((int)__brkval);
  }
  return free_memory;
}

// parameter setting for solver 
class QPparams {
public:
    QPparams()
        : max_iter_outer(5),
          max_iter_newton(10),
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

// Set all matrix elements to zero 
template<int Rows, int Cols, typename DType>
void AllZero(Matrix<Rows, Cols, DType>& mat) {
    mat.Fill(static_cast<DType>(0)); 
}

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
        QP(const String& mode = "") 
            : _params(new QPparams()), _alloc(true) {
            QPZero();
            if (mode=="D") _debugging=true;
            else _debugging=false;       
        };

        // Destructor
        ~QP(){
            if(_alloc) delete _params;
        }

        // Constructor accepting a QPparams object
        QP(const QPparams& parameters, const String& mode = "") 
            : _params(&parameters), _alloc(false){
            QPZero();
            if (mode=="D") _debugging=true;
            else _debugging=false;                                             
        }

        // Move constructor 
        QP(QP&& other) noexcept
            : _params(other._params),
              _alloc(other._alloc),
              _debugging(other._debugging),
              _Q(other._Q),
              _q(other._q),
              _A(other._A),
              _b(other._b),
              _G(other._G),
              _h(other._h){
            
            other._params = nullptr;
            other._alloc = false;
            // reset other matrices to zero
            QPZero(other);
        };

        // Move assignment operator
        QP& operator=(QP&& other) noexcept {
            if (this != &other) {

                if (_alloc && _params != nullptr) {
                    delete _params; 
                }

                _params = other._params;
                _alloc = other._alloc;
                _debugging = other._debugging;

                _Q = other._Q;
                _q = other._q;
                _A = other._A;
                _b = other._b;
                _G = other._G;
                _h = other._h;

                other._params = nullptr;
                other._alloc = false;
                // reset other matrices to zero
                QPZero(other);
            }
            return *this; 
        }

        // Copy constructor
        QP(const QP& other)
            : _alloc(other._alloc),  
            _debugging(other._debugging),  
            _Q(other._Q),  
            _q(other._q),
            _A(other._A),
            _b(other._b),
            _G(other._G),
            _h(other._h){
            _params = new QPparams(*other._params);  
        }

        // Copy assigment operator
        QP& operator=(const QP& other) {
            if (this != &other) { 

                if (_params != nullptr) {
                    delete _params; 
                }

                _alloc = other._alloc;
                _debugging = other._debugging;
                _Q = other._Q;  
                _q = other._q;
                _A = other._A;
                _b = other._b;
                _G = other._G;
                _h = other._h;

                _params = new QPparams(*other._params); 
            }
            return *this;  
        }

        // solving the qp
        QPsol<nx,m,p,DType> solve(const String& mode = ""){

            if (_debugging){
                Serial.print("Free RAM when solving the QP: ");
                Serial.print(freeMemory());            
                Serial.println(" kB");
            }

            bool verb = false;
            if (mode=="verbose") verb=true;

            Matrix<nx,1,DType> x;
            Matrix<m ,1,DType> lambda;
            Matrix<p ,1,DType> mu;     

            AllZero(x);
            AllZero(lambda);
            AllZero(mu);

            DType rho = _params -> penalty_initial;
            DType Phi = _params -> penalty_scaling;

            Matrix<m+p,1,DType> pr;
            Matrix<nx,1,DType> Nabla_x_L;
            double prnorm;

            DType obj_val; 
            SolverStatus status;

            for (int i = 0; i < _params -> max_iter_outer; i++){
                x = newton_solve(x, lambda, mu, rho);
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
        bool _debugging; 
        
        // initializes all matrices of 'this' with zeros 
        void QPZero(){
            AllZero(_Q);
            AllZero(_q);
            AllZero(_A);
            AllZero(_b);
            AllZero(_G);
            AllZero(_h);
        }

        // initializes all matrices of qp with zeros 
        void QPZero(QP& qp){
            AllZero(qp._Q);
            AllZero(qp._q);
            AllZero(qp._A);
            AllZero(qp._b);
            AllZero(qp._G);
            AllZero(qp._h);
        }

        // returns residual vector of the equality constraints   
        inline Matrix<m,1, DType> c_eq(const Matrix<nx,1,DType>& x){
            return _A*x - _b; 
        };
        
        // returns residual vector of the inequality constraints   
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
            if(m != 0){
                g += (~_A*rho) * c_eq(x);    
            }
            if(p != 0){
                Matrix<p,p,DType> Ip = active_ineq(x, mu, rho);
                g += (~_G*Ip) * c_in(x);
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
        Matrix<nx, 1, DType> newton_solve(const Matrix<nx,1,DType>& x, const Matrix<m,1,DType>& lambda, const Matrix<p,1,DType>& mu, const DType& rho){
            
            Matrix<nx,1,DType> x_sol = x;            
            Matrix<nx,1,DType> Deltax; 

            // preallocate gradient and Hessian of augmented Lagrangian 
            Matrix<nx,1,DType> g; 
            Matrix<nx,nx,DType> H; 

            for (int i = 0; i < _params -> max_iter_newton; i++){

                algradient(g, x_sol, lambda, mu, rho);
                DType gnorm = residual_norm(g);

                if (gnorm < _params -> precision_newton){
                    return x_sol;
                }

                alhessian(H, x_sol, lambda, mu, rho);
                
                auto Hdecomp = LUDecompose(H);
                Deltax = LUSolve(Hdecomp, -g);

                // If Hessian of augmented Lagrangian is rank defficient abort procedure 
                if (Hdecomp.singular){
                    x_sol.Fill(0.0/0.0);
                    return x_sol;                    
                }

                // simple back tracking line search on the AL gradient residual 
                DType alpha = 1.0;                  // scaling parameter 
                Matrix<nx,1,DType> x_backtrack;     // solution after taking reduced newton step
                DType gnorm_red;                    // norm of AL gradient residual after taking reduced newton step 
                for (int i = 0; i < _params -> max_iter_backtrack; i++){
                    // solution with reduced stepsize
                    x_backtrack = x_sol+alpha*Deltax;
                    
                    // compute residual with reduced stepsize 
                    algradient(g, x_backtrack, lambda, mu, rho);
                    gnorm_red = residual_norm(g);
                        
                    // if residual before newton step is higher than after, reduce the scaling parameter alpha   
                    if (gnorm_red > gnorm)
                        alpha *= _params -> backtrack_beta; 
                    else 
                        break;
                }
                if (_debugging){
                    Serial.print(F("Newton iteration "));
                    Serial.print(i);
                    Serial.print(F(": norm of gradient residual: ")); 
                    Serial.print(gnorm_red, 10); 
                    Serial.print(F(" vs ")); 
                    Serial.print(_params -> precision_newton, 10); 
                    Serial.println(F(" (threshold) ")); 
                }        
                x_sol = x_backtrack;
            }
            return x_sol;
        };
    };
}; // namespace ALQPS
