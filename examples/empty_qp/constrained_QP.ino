#include "ALQPSolver.h"
using namespace ALQPS;

void setup() {

    Serial.begin(9600); // Starting Serial Terminal
    
    // QP dimensions
    const int n = 2;
    const int m = 2;
    const int p = 1;

    // quadratic objective
    Matrix<n,n> Q = {1, 0, 0, 1};
    Matrix<n,1> q = {1, -1};
    // equality constraints 
    Matrix<m,n> A = {0,1,1,0};
    Matrix<m,1> b = {0,1};
    // inequality constraints  
    Matrix<p ,n> G = {0.5,0};
    Matrix<p ,1> h = {0.5};


    // build and update the quadratic program 
    QP<n,m,p> qp; 
    qp.update(Q,q,A,b,G,h); 

    // solve QP 
    auto sol = qp.solve("verbose"); 
    // access solution and value of objective at optimum 
    Matrix<n,1> x_opt = sol.x;
    float obj_val = sol.obj_val;
}
    

void loop(){
    
}