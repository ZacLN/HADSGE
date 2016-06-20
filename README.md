# HADSGE.jl

A package for solving Heterogeneous Agent Dynamic Stochastic General Equilibrium Models.

## Algorithm
This packages uses a coupled iteration scheme to solve nonlinear economic models that allows for a dependancy between agent's policy functions (i.e. households subject to idiosycratic shocks) and aggregate state variables' laws of motion. The inner loop/iteration solves the model equations using a backwards time iteration scheme that takes advantage of Julias function ASTs to provide exact and efficient calculation of the model Jacobian matrix. The outer loop uses Kolmogorov forward equation to find exact representations of the aggregate distribution of variables and thus aggregate laws of motion. This removes the need for the time consuming simulation required by most solution methods for economic models with heterogeneous agents.

## Model
Model specification requires 4 components.
1. A set of equations (laws of motion and equilibrium conditions).
2. The parameters used in the model.
3. The state space variables, both endogenous and exogenous/stochastic.
4. The definition/bounds of policy and other variables used in the model.

### Equations
n equations F() written such that F = 0. Variables may have time indices [-1],[0] or [+1] which indicates the period in which the variables value is first known. Future variables (e.g. x[+1]) must be enclosed in an `Expect()` function to indicate that expectations of that expression are to be calculated. Variables indexed by `[-1]` are endogenous state variables. The equations may contain the following mathematical functions:
* +
* -
* *
* /
* ^
* exp
* log
* max
* min

### Parameters
A simple list of parameters used in the model. Unicode characters may be used, (e.g. `α = 0.33`).

### State Variables
These can be one of three types: endogenous, autoregressive or finite state Markov.
* **Endogenous:** `s=(lb,ub,n)` defines an endogenous state variable with lower/upper bounds given by lb/ub and the density of grid nodes in this dimension controlled by `n`. `n` indicates the level of the sparse grid on which the model is solved in this dimension.
* **Autoregressive:** `e=(mu,rho,sigma,n)` defines a process such that `e[t+1] = (1-rho)*mu + rho*e[t] + sigma*epsilon` with grid density again specified by `n`. `epsilon` is a white noise process.
* **Markov:** `e = ([lb,ub],T,n)` defines a finite state Markov process with equally spaced states between `lb` and `ub`. Transition matrix `T` must be of size `[m,m]` where `m` is given as:
  * `m = 1     if n = 0`
  * `m = 3     if n = 1`
  * `m = 2^n+1 if n > 1`

Variables are defined as aggregate by placing a colon before the equals symbol:

`x:=(mu,rho,sigma,n)`
### Other variables

* **Policy:** `p=(lb,ub,init)` defines a policy variable bounded by `[lb,ub]` with the initial policy function defined by `init`. This expression can only contain state variables.
* **Dependant:** `x = expr` defines a placeholder variable that will be substituted into the model equations. For example `lambda = c^-2.5` leads to all instances of `lambda` being replaced in the set of equations. All instances of `lambda[+1]` are replaced by `c[+1]^-2.5`. Dependant variables must be defined *in order* so that one may define `lambda2 = lambda*2` after but not before the previous declaration.
* **Exogenous:** `x = 0.4` initialises the variable `x` with the same value over the entire grid. Following the declaration of model `M` this variable can be altered using `M[:x] = V` where V is either a number or vector of size `length(M)` and the model resolved with this new value.
* **Aggregate:** `X=∫(x,v)` defines an aggregate variable `X` whose value is initialised as number `v`. `x` defines the expression or variable that is integrated to give the aggregate variable. Following the solution of the model `M`, `X` can be updated by calling `updateA(M)`.


## Solve

`solve(M,maxiter=1000,ϕ=0.8;crit = 1e-6)`

Solves the equation of the model until the maximum equation errors are smaller than `crit` or the number of iterations exceeds the limit. `ϕ` is a dampening parameter that can be set between `0` and `1` where the lower value results in no changes to policy functions.

`updateA(M)`

Updates aggregate variables of agents over the state space. First it constructs a transition matrix over a tensor grid of the state grid.  The transition matrix is then used to try to find a stable distribution over the tensor grid.

## Model Inspection
Model variables can be accessed using the name and time index,`M[:x,ts]`. The distribution of agents over the state space is located in `M.distribution.d`








## Example
```
M=Model(:[
        1-R*β*Expect(λ[+1])/λ
        Uh-λ*η*W
],:[b       = (-2,10.,8)
    η       = (1,0.9,0.1,1)
    W       = (1,0.9,0.01,1)
],:[b       = (-2,10.,b*0.95)
    h       = (0,1,0.7)
    c       = W*h*η+R*b[-1]-b
    λ       = c^-σc
    Uh      = ϕh*(1-h)^-σh
    B       = ∫(b,0.0)
    H       = ∫(h*η,0.3)
],:[β       = 0.98
    σc      = 2.5
    ϕh      = 2.0
    σh      = 2.0
    R       = 1.0166])
solve(M,disp=100)
updateA(M)
scatter(M[:b,-1],M[:c,0])
plot(M[:b].x,sum(M.distribution.d,2))
```
