# Exercise 1
Consider the initial value problem

- *y*'(*t*) = 1 - *y*<sup>2</sup>, *t* ∈ (0, 20),  
- *y*(0) = 0.

The exact solution is *y*(*t*) = tanh(*t*) = (e<sup>2*t*</sup> - 1)/(e<sup>2*t*</sup> + 1).

The solution can be approximated on a set of *N*+1 equispaced nodes {*t*<sub>0</sub>, *t*<sub>1</sub>, **...** ,*t*<sub>*N*</sub>}, where *t*<sub>0</sub> = 0 and *t*<sub>*N*</sub> = 20, *h* = 20/*N* is the step size, and *y*(*t*<sub>*n*</sub>) is approximated by *u*<sub>*n*</sub>.

**(a)** First, we solve the problem using the forward Euler method, for *N* ∈ {19, 21, 40, 80, 160}, and then plotting the 5 approximations along with the exact solution for comparison. We'll also store the maximum error for each *N* to compare later.

**(b)** Doing the same thing as in **(a)**, but using Heun's method.

**(c)** We'll now plot the errors of each method against *h* to compare.

**(d)** Some comments on the stability of the methods used.

# Exercise 2
For this exercise, we'll be implementing our own forward Euler, backward Euler and Crank-Nicolson methods.

Consider the initial boundary value problem:

Given *T* > 0, find *u* : [0, 1] x [0, *T*] → ℝ such that

- ∂*u*/∂*t* - ∂<sup>2</sup>*u*/∂*t*<sup>2</sup> = 1, in [0, 1] x [0, *T*],  
- *u*(0, *t*) = *u*(1, *t*) = 0 for all *t* ∈ [0, *T*],  
- *u*(*x*, 0) = sin(π*x*) + *x*(1 - *x*)/2.

The exact solution is *u*(*x*) = *e*<sup>-π<sup>2</sup>*t*</sup>sin(π*x*) + *x*(1 - *x*)/2.

Similar to exercise 1, we choose a set of *N*+1 equispaced nodes {*x*<sub>0</sub>, *x*<sub>1</sub>, **...** ,*x*<sub>*N*</sub>}, where *x*<sub>0</sub> = 0, *x*<sub>*N*</sub> = 1, and *h* = 1/*N*, and approximate *u*(*x*<sub>*n*</sub>, *t*) as *u*<sub>*n*</sub>(*t*).

We obtain the system of ordinary differential equations

*d***u**/*dt* = **f**(*t*, **u**(*t*)) = -*A***u** + **b**, *t* ∈ [0, *T*], where **u** = (*u*<sub>1</sub>, **...** , *u*<sub>*N*-1</sub>)<sup>*T*</sup>.

But this time, we choose another set of *M*+1 equispaced nodes {*t*<sub>0</sub>, *t*<sub>1</sub>, **...** ,*t*<sub>*M*</sub>}, where *t*<sub>0</sub> = 0, *t*<sub>*M*</sub> = *T*, and *k* = *T*/*M*. We now approximate *u*<sub>*n*</sub>(*t*<sub>*m*</sub>) as *u*<sub>*n*,*m*</sub>.

**(a)** We'll start off by solving the problem using each of the three methods, with *N* = 40, *k* = *h*, and *T* = 0.2, and then plotting four graphs. The first three will contain each approximation, the fourth will be the exact solution on its own.

**(b)** Now, using *N* ∈ {10, 20, 40, 80, 160}, *k* = *h* for each *N*, and *T* = 0.2, we'll plot the root-mean-squared error of each method on three separate graphs, and examine the convergence of the backward Euler and Crank-Nicolson methods.

**(c)** Looking further into the error of the two methods, we'll plot the root-mean-squared error for when *h* is small compared to *k*, and then when *k* is small compared to *h*.

**(d)** Returning to the forward Euler method, we'll examine the spectral radius of *A* for each *N* ∈ {10, 20, 40, 80, 160, 320}.

# Concepts
Here are some brief explanations of the concepts used in these exercises.

Suppose we're faced with the problem

- *y*'(*t*) = *f*(*t*, *y*(*t*)), *t* > 0,
- *y*(0) = *y*<sub>0</sub>

### Forward Euler method
This method uses the following iteration

- *u*<sub>0</sub> = *y*<sub>0</sub>
- *u*<sub>*n*+1</sub> = *u*<sub>*n*</sub> + *h*<sub>*n*</sub>*f*(*t*<sub>*n*</sub>, *u*<sub>*n*</sub>)

Note: *h*<sub>*n*</sub> = *t*<sub>*n*+1</sub> - *t*<sub>*n*</sub>, but in these exercises we use a constant value for *h* so the labelling is omitted.

### Backward Euler method
Similar to the forward Euler method but with a small difference

- *u*<sub>0</sub> = *y*<sub>0</sub>
- *u*<sub>*n*+1</sub> = *u*<sub>*n*</sub> + *h*<sub>*n*</sub>*f*(*t*<sub>*n*+1</sub>, *u*<sub>*n*+1</sub>)

### Crank-Nicolson method
The Crank-Nicolson method uses a mix of the two Euler methods

- *u*<sub>0</sub> = *y*<sub>0</sub>
- *u*<sub>*n*+1</sub> = *u*<sub>*n*</sub> + *h*<sub>*n*</sub>(*f*(*t*<sub>*n*</sub>, *u*<sub>*n*</sub>)+*f*(*t*<sub>*n*+1</sub>, *u*<sub>*n*+1</sub>))/2

### Heun's method
A variation on the Crank-Nicolson method

- *u*<sub>0</sub> = *y*<sub>0</sub>
- *u*<sub>*n*+1</sub> = *u*<sub>*n*</sub> + *h*<sub>*n*</sub>(*f*(*t*<sub>*n*</sub>, *u*<sub>*n*</sub>)+*f*(*t*<sub>*n*+1</sub>, *h*<sub>*n*</sub>*f*(*t*<sub>*n*+1</sub>, *u*<sub>*n*+1</sub>)))/2

Notice that this simply replaces the *u*<sub>*n*+1</sub> term from the Crank-Nicolson method with the definition in the forward Euler method.
