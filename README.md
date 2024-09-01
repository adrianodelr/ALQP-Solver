# Augmented Lagrangian Quadratic Program Solver
This library provides an Augmented Lagrangian solver for Arduino projects. It lets you solve quadratic programs (QPs) with optional linear equality and inequality constraints, and can even warn you about infeasible problems. The library is super lightweight since all matrix operations are handled by the [BasicLinerAlgebra library](https://github.com/tomstewart89/BasicLinearAlgebra), so it should work on any Arduino board. 

The Augmented Lagrangian formulation is based on one of the homeworks I did when following online [the CMU optimal control and reinforcement learning course](https://optimalcontrol.ri.cmu.edu/) by Prof. Zachary Manchester et al., which I highly recommend to anyone - it's an awesome course! But if you're more interested in jumping straight into solving constrained quadratic programs, I've laid out the basics of the problem formulation below, along with a guide on how to use the library. 

For more help, you can also check out [this tutorial given by Kevin Tracy](https://www.youtube.com/watch?v=0x0JD5uO_ZQ) as supplementary material to the homework. 

## Constrained QP formulation
$$\begin{align}
\min_{\mathbf{x}} \quad & \frac{1}{2}\mathbf{x}^T\mathbf{Q}\mathbf{x} + \mathbf{q}^T\mathbf{x} \\ 
\mbox{s.t.}\quad &  \mathbf{A}\mathbf{x} -\mathbf{b} = \mathbf{0} \\ 
&  \mathbf{G}\mathbf{x} - \mathbf{h} \leq \mathbf{0} 
\end{align}$$

Here $\mathbf{Q} \in \mathbb{R}^{n \times n}$ and $\mathbf{q}\in \mathbb{R}^{n}$ are quadratic and linear coefficient matrix, respectively, which determine our quadratic objective. Equality constraints are represented by matrix $\mathbf{A}\in \mathbb{R}^{m \times n}$ and vector $\mathbf{b}\in \mathbb{R}^m$, inequality constraints by matrix $\mathbf{G} \in \mathbb{R}^{p \times n}$ and vector $\mathbf{h}\in \mathbb{R}^p$. (This notation also is used throughout the code, and is conform with the notation used in the tutorial linked above.)

## Usage
Defining a QP is fairly simple. Matrices can be built according to the documentation of the Basic Linear Algebra library linked above:
```cpp
// QP dimensions
const int n = 2;
const int m = 2;
const int p = 2;

// quadratic objective
Matrix<n,n> Q = {2, 0, 2, 0};
Matrix<n,1> q = {-4, -6};
// equality constraints 
Matrix<m ,n> A = {1, 1, 0.3, -1};
Matrix<m ,1> b = {3, 1};
// inequality constraints  
Matrix<p ,n> G = {-1, 1, 0, -1};
Matrix<p ,1> h = {1, 1};
```
The QP is built as follows and can updated sequentially at runtime, as e.g. required by a model predictive controller:  

```cpp
// build and update the quadratic program 
QP<nx,m,p> qp; 
qp.update(Q,q,A,b,G,h); 
```
Solving the QP will return a solution object which, alongside the solution, contains values of the according dual variables and information about success of the solver. 
```cpp
// solve QP 
auto sol = qp.solve(); 
// access solution and value of objective at optimum 
Matrix<n,1> x_opt = sol.x;
float obj_val = sol.obj_val;
// 
```
It is also possible to print a solver status report by setting a  default argument to verbose":        
```cpp
// solve QP 
auto sol = qp.solve("verbose"); 
```
This will result in following output printed to the console 
```bash
status: solved
optimal objective: -2.62
primal value (solution): [[3.08],[-0.08]]
dual value (equalities): [[-1.69],[-1.54]]
dual value (inequalities): [[0.00],[0.00]]
```

