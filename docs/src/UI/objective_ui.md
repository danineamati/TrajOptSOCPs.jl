# User Interface to Present the Objective Function

```@meta
CurrentModule = TrajOptSOCPs
```

### Linear Quadratic Regulator
To avoid generating trajectories too different from the initial trajectory,
we can use a Linear Quadratic Regulator as an objective. This has the form:

```math
\begin{aligned}
& \frac{1}{2} (x_N - x_{N, ref})^\top Q (x_N - x_{N, ref}) + \\
& \frac{1}{2} \Sigma_{n = 1}^{N-1} (x_n - x_{n, ref})^\top Q (x_n - x_{n, ref})
+ (u_n - u_{n, ref})^\top R (u_n - u_{n, ref}))  \\
\end{aligned}
```

Currently the module only supports Q and R matrices that are identical at
every timestep, but this can be changed in the future. However, the framework
does allow a QP that only uses the Q matrix on the last timestep and ignores
it elsewhere.

In the example below, we build a simple diagonal Q and R matrix.
```@example
using LinearAlgebra #hide
const rocketStart = [2.0; 20.0; 0.0; -15.0] # hide
lqrQMat = 0.0001 * Diagonal(I, size(rocketStart, 1))
lqrRMat = 0.0025 * Diagonal(I, Int64(size(rocketStart, 1) / 2))
nothing #hide
```

We reference the trajectory to the initial trajectory.
```@example
using LinearAlgebra, TrajOptSOCPs #hide
const rocketStart = [2.0; 20.0; 0.0; -15.0] # hide
const rocketEnd = [0.0; 0.0; 0.0; 0.0] # hide
uHover = [0.0, -10.0] # hide

lqrQMat = 0.0001 * Diagonal(I, size(rocketStart, 1)) # hide
lqrRMat = 0.0025 * Diagonal(I, Int64(size(rocketStart, 1) / 2)) # hide
const NSteps = 60 # hide

initTraj = initializeTraj(rocketStart, rocketEnd, uHover, uHover, NSteps)

costFun = makeLQR_TrajReferenced(lqrQMat, lqrRMat, NSteps, initTraj)
```

For more information, see
```@docs
LQR_QP_Referenced
```
