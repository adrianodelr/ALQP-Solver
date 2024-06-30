#include "ALQPSolver.h"
using namespace ALQPS;

void setup() {

   Serial.begin(9600); // Starting Serial Terminal

    const int nx = 5;
    const int m = 5;
    const int p = 5;

    Matrix<nx,nx> Q;
    Matrix<nx,1>  q;
    Matrix<p ,nx> G;
    Matrix<p ,1>  h;
    
    Q.Fill(0.);
    q.Fill(0.);
    G.Fill(0.);
    h.Fill(0.);

    QP<nx,m,p> qp; 
    qp.update(Q,q,{},{},G,h); 

    auto xsol = qp.solve();   
    Serial.print("Solution: ");
    Serial.println(xsol);    
}

void loop(){
    
}