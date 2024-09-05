using OSQP, SparseArrays, LinearAlgebra 

# Cost function 
Q = [1 0; 0 1] 
# Q = Q'Q     # ensure matrix is SPD 
q = [1;-1]  
# # equality constraints 
# A = [0 1; 1 0]
# b = [0;1] 

# # inequality constraints 
# G = [0.5 0]
# # h = [1;0]
# h = [0.5]
 

G = []
h = []

A = []
b = [] 

 
lu = [b; ] 
uu = [b; h] 
    
solver = OSQP.Model() 
OSQP.setup!(solver, P=sparse(Q), q=Float64.(q), A=SparseMatrixCSC([A;G]), l=Float64.(lu), u=Float64.(uu), verbose=true,warm_start=false)
u = OSQP.solve!(solver).x     
     
sol_ = OSQP.solve!(solver)

solver = OSQP.Model() 
OSQP.setup!(solver, P=sparse(Q), q=Float64.(q), A=SparseMatrixCSC(G), l=Float64.([-Inf; -Inf]), u=Float64.(h), verbose=true,warm_start=false)
u = OSQP.solve!(solver).x  
  

solver = OSQP.Model() 
OSQP.setup!(solver, P=sparse(Q), q=Float64.(q),verbose=true,warm_start=false)
u = OSQP.solve!(solver).x     
