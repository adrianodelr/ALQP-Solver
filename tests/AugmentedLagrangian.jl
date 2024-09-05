using Random, LinearAlgebra, Plots

mutable struct QPsol
    x::Vector
    λ::Vector
    μ::Vector
    solved::Bool
    pinf::Bool
    dinf::Bool
    objective::Float64
    function QPsol(x::Vector,λ::Vector,μ::Vector)
        return new(x,λ,μ,false,false,false,0.0)
    end
end 

struct QPdata
    Q::Matrix
    q::Vector
    A::Matrix
    b::Vector
    G::Matrix
    h::Vector
    sol::QPsol
    function QPdata(Q::Matrix,q::Vector,A::Matrix,b::Vector,G::Matrix,h::Vector)
        return new(Q,q,A,b,G,h,QPsol(zero(q),zero(b),zero(h)));      
    end
end; 


function QPdata(n::Int,m::Int,p::Int)
    """
    n=state dimension
    m=number of equality constraints
    p=number of inequality constraints
    """
    QPdata(zeros(n,n),zeros(n),zeros(m,n),zeros(m),zeros(p,n),zeros(p));      
end

num_eq(qp::QPdata) = length(qp.b);
num_ineq(qp::QPdata) = length(qp.h);
Base.size(qp::QPdata) = (length(qp.q), num_eq(qp), num_ineq(qp)); 
objective(qp::QPdata, x) = 0.5 * x'qp.Q*x + qp.q'x;
ceq(qp::QPdata, x) = qp.A * x - qp.b;
cin(qp::QPdata, x) = qp.G * x - qp.h;

function update_solution!(qp::QPdata,x,λ,μ; pinf=false,dinf=false,solved=true)
    qp.sol.solved = solved
    qp.sol.λ = λ
    qp.sol.μ = μ
    if solved
        qp.sol.objective = objective(qp, x)
        qp.sol.x = x
    else
        qp.sol.objective = NaN
        fill!(qp.sol.x,NaN)
    end 
    qp.sol.pinf = pinf
    qp.sol.dinf = dinf
end 


function stationarity(qp::QPdata, x, λ, μ)
    # println("qp.Q*x: $(qp.Q*x), qp.q: $(qp.q); qp.A'λ: $(qp.A'λ), qp.G'μ: $(qp.G'μ)")
    p = qp.Q*x + qp.q + qp.A'λ + qp.G'μ;
    return p
end 

function primal_residual(qp::QPdata, x, λ, μ)
    c = ceq(qp, x)
    h = cin(qp, x)
    return [c; max.(0, h)]
end 

function dual_feasibility(qp::QPdata, μ)
    return all(μ .>= 0)      
end  

function complimentarity(qp::QPdata, x, λ, μ)
    return [min.(0, μ); μ .* cin(qp, x)]
end 
function active_set(qp::QPdata, x, ρ, μ)
    p = size(qp.h,1) 
    Ip = zeros(p,p)
    h = cin(qp, x)
    for (hi,μi,i) in zip(h,μ,eachindex(h))
        if hi < 0 && μi == 0
            Ip[i,i]=0
        else 
            Ip[i,i]=ρ 
        end      
    end 
    return Ip
end 
 
"""
algrad(qp, x, λ, μ, ρ)
Compute the gradient of the augmented Lagrangian, provided the QP data `qp`, penalty parameter `ρ`,
primal variables `x`, equality Lagrange multipliers `λ` and inequality Lagrange multipliers `μ`
"""
function algrad(qp, x, λ, μ, ρ)
    ∇xL = stationarity(qp, x, λ, μ)
    # println("stationarity: $(∇xL)")
    Ip = active_set(qp, x, ρ, μ)
    # println("primal residual: $(primal_residual(qp, x, λ, μ))")
    ∇xg = [ρ*qp.A' qp.G'Ip]*primal_residual(qp, x, λ, μ)    
    return ∇xL+∇xg
end  
"""
alhess(qp, x, λ, μ, ρ)
Compute the Hessian of the augmented Lagrangian, provided the QP data `qp`, penalty parameter `ρ`,
primal variables `x`, equality Lagrange multipliers `λ` and inequality Lagrange multipliers `μ`
"""
function alhess(qp, x, λ, μ, ρ)
    ∇xxL = qp.Q
    Ip = active_set(qp, x, ρ, μ)
    ∇xxg = [ρ*qp.A' qp.G'Ip]*[qp.A;qp.G]
    return ∇xxL+∇xxg
end
"""
newton_solve(qp, x, λ, μ, ρ; kwargs...)
Minimize the augmented Lagranginan given the current values of the dual variables `λ` and `μ` and the penalty parameter `ρ`.
"""
function newton_solve(qp::QPdata, x, λ, μ, ρ; eps_inner=1e-6)
    for i = 1:10
    # Compute the gradient and the Hessian of the augmented Lagrangian
        # println(" Iteration: $i")
        g = algrad(qp, x, λ, μ, ρ)
        # println("g: $(g)")
        if norm(g) < eps_inner
            return x
        end
        H = alhess(qp, x, λ, μ, ρ)
        # println("H: $(H)")
        # println("inner loop: ", i)
        # println("Hinv: $(inv(H))")
        Δx = -H\g
        x += Δx
        # println("Newton ",i,":", "normg: ", norm(g), " precision_newton: ", eps_inner , " Deltax ", Δx, " x: ", x)
        # println("x_sol $x")
    end
    @warn "Inner solve max iterations"
    return x
end
"""
dual_update(qp, x, λ, μ, ρ)
Update the dual variables `λ` and `μ` give the primal variables `x`, QP data `qp` and penalty parameter `ρ`.
"""
function dual_update(qp, x, λ, μ, ρ)
    μnext = max.(0, μ + ρ*cin(qp, x))
    λnext = λ + ρ.*ceq(qp, x)
    return λnext, μnext
end
"""
    solve_qp(qp::QPData, x0, [λ0, μ0]; kwargs...)

Solve the quadratic program (QP) specified by `qp::QPData`, given initial guess `x` for the primal variables, 
and optionally the Lagrange multipliers for the equality `λ` and inequality `μ` constraints.

Returns the optimized solution of primal and dual variables, `xstar,λstar,μstar`.

# Optional Keyword Arguments
* `penalty_initial` initial value of the penalty parameter
* `penalty_scaling` geometric scaling factor for the penalty updates
* `eps_primal` tolerance for primal feasiblity (constraint violation)
* `eps_inner` tolerance for inner Newton solve
* `max_iters` maximum number of outer loop iterations
"""
function solve_qp(qp::QPdata; 
        penalty_initial=1.0, 
        penalty_scaling=10.0, 
        eps_primal=1e-8,
        eps_inner=1e-8,
        max_iters=10
    )
    n,m,p = size(qp)
    x=zeros(n) 
    λ=zeros(m) 
    μ=zeros(p)

    ρ = penalty_initial
    ϕ = penalty_scaling
    
    # Start outer loop
    for i = 1:max_iters
                
        x  = newton_solve(qp, x, λ, μ, ρ, eps_inner=eps_inner)

        λ, μ = dual_update(qp, x, λ, μ, ρ)
        ρ = ϕ*ρ
        
        if norm(primal_residual(qp, x, λ, μ)) < eps_primal && dual_feasibility(qp,μ)
            if all(x->isapprox(x,0,atol=1e-6),stationarity(qp, x, λ, μ))  
                println("found globally optimal solution")
            end
            update_solution!(qp,x,λ,μ)
            return qp.sol
        end        
    end
    if norm(primal_residual(qp, x, λ, μ)) > eps_primal 
        println("Problem is primal infeasible!")
        update_solution!(qp,x,λ,μ,pinf=true,solved=false)
    end 
    if !dual_feasibility(qp,μ)
        println("Problem is dual infeasible!")
        update_solution!(qp,x,λ,μ,dinf=true,solved=false)
    end
    return qp.sol 
end  
  
 
n,m,p=2,2,2;  
qp = QPdata(n,m,p);    
qp.Q.=Q; 
qp.q.=q;
qp.A.=A;
qp.b.=b;
qp.G.=G;
qp.h.=h;  
sol = solve_qp(qp, max_iters=10)  
sol.x 
sol.λ  
sol.μ 

# eps_primal=1e-6 
# eps_inner=1e-4
# max_iters=50

# penalty_initial=1.0 
# penalty_scaling=10.0 
 
# ρ = penalty_initial
# ϕ = penalty_scaling 



# n,m,p = size(qp) 
# x=zeros(n) 
# λ=zeros(m) 
# μ=zeros(p)
 
# g = algrad(qp, x, λ, μ, ρ)   


# if norm(g) < eps_inner
#     λ, μ = dual_update(qp, x, λ, μ, ρ)
#     ρ = ϕ*ρ
# end  

# H = alhess(qp, x, λ, μ, ρ) 
 
# Δx = -H\g   
# x += Δx   

n,m,p=2,0,0;  
qp = QPdata(n,m,p);    
qp.Q.=Q; 
qp.q.=q;
qp.A.=zeros(m,n);  
qp.b.=zeros(m,1); 
qp.G.=zeros(p,n); 
qp.h.=zeros(p,1);    

sol = solve_qp(qp, max_iters=25)               
x = sol.x  
 
# x, λ, μ, qp = let  
#     n,m,p=2,2,2;  
#     qp = QPdata(n,m,p);  
#     qp.Q.=randn(n,n);
#     qp.Q.=qp.Q'qp.Q;  # make sure is symmetric positive definite
#     isposdef(qp.Q) 
#     qp.q.=randn(n,1);
#     qp.A.=randn(m,n);
#     qp.b.=randn(m,1);
#     qp.G.=randn(p,n);
#     qp.h.=randn(p,1);
#     # example solving qp 
#     x, λ, μ = solve_qp(qp, max_iters=25)
#     x, λ, μ, qp
# end     
# n,m,p=2,0,2;  
# qp = QPdata(n,m,p);  
# qp.Q.=randn(n,n);
# qp.Q.=qp.Q'qp.Q;  # make sure is symmetric positive definite
# isposdef(qp.Q) 
# qp.q.=randn(n,1);
# # # qp.A.=randn(m,n);
# # # qp.A.=[];
# # qp.b.=randn(m,1);
# qp.G.=randn(p,n);
# qp.h.=randn(p,1);
# # example solving qp 
# x, λ, μ = solve_qp(qp)



# let 
#     Nsamp = 45
#     Xsamp = LinRange(-100,100,Nsamp)
#     Ysamp = LinRange(-100,100,Nsamp)
#     Zsamp = zeros(Nsamp,Nsamp)

#     for (ix,x) in enumerate(Xsamp)
#         for (iy,y) = enumerate(Ysamp)
#             Zsamp[ix,iy] = objective(qp, [x; y])
#             #f([x; y])
#         end
#     end
#     #contour(Xsamp,Ysamp,Zsamp)
#     h = heatmap(Xsamp,Ysamp,Zsamp, alpha = 0.8)

#     # xc = LinRange(-10,10,Nsamp)
#     # plot!(xc,xc.+2, linewidth=3, alpha = 1, colour=:green, label = "cin")
#     # plot!(xc,-xc.-1, linewidth=3, alpha = 1, colour=:blue, label = "ceq")
# end 
