# Augmented Lagrangian Quadratic Program Solver
This library provides an Augmented Lagrangian solver for Arduino projects. It allows to solve 
quadratic programs with linear equality and inequality constraints (both optional) and optionally can warn the user about infeasible problems. It is extremely lightweight as all linear algebra operations are done only with help of the [BasicLinerAlgebra library](https://github.com/tomstewart89/BasicLinearAlgebra), it therefore should be able run on any Arduino board. The Augmented Lagrangian formulation is very much based on one of the homeworks I did when following online [the CMU optimal control and reinforcement learning course](https://optimalcontrol.ri.cmu.edu/) by Prof. Zachary Manchester et al.. I strongly recommend the course to anyone, it is really awesome! However, for those who like to jump straightforward into solving constraint quadratic programs, below the basics of the problem formulation are given along with a guide on how to use the library. A further very useful video is [this tutorial given by Kevin Tracy](https://www.youtube.com/watch?v=0x0JD5uO_ZQ) as supplementary material to the homework. 






