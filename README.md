# Augmented Lagrangian Quadratic Program Solver
This library provides an Augmented Lagrangian solver for Arduino projects. It lets you solve quadratic programs (QPs) with optional linear equality and inequality constraints, and can even warn you about infeasible problems. The library is super lightweight since all matrix operations are handled by the [BasicLinerAlgebra library](https://github.com/tomstewart89/BasicLinearAlgebra), so it should work on any Arduino board. 

The Augmented Lagrangian formulation is very much based on one of the homeworks I did when following online [the CMU optimal control and reinforcement learning course](https://optimalcontrol.ri.cmu.edu/) by Prof. Zachary Manchester et al., which I highly recommend to anyone - it's an awesome course! But if you're more interested in jumping straight into solving constrained quadratic programs, I've laid out the basics of the problem formulation below, along with a guide on how to use the library. 

For more help, you can also check out [this tutorial given by Kevin Tracy](https://www.youtube.com/watch?v=0x0JD5uO_ZQ) as supplementary material to the homework. 

## Constrained QP formulation
$$\begin{align}
\min_x \quad & \frac{1}{2}x^TQx + q^Tx \\ 
\mbox{s.t.}\quad &  Ax -b = 0 \\ 
&  Gx - h \leq 0 
\end{align}$$



