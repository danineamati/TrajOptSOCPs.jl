# User Interface to Solve the Optimization Problem

```@meta
CurrentModule = TrajOptSOCPs
```

The solver requires a set of parameters. Thereafter, we can just run the solver
using the problem defined by the augmented lagrangian.

### Solver Parameters
The solver has many knobs. The recommended settings are below:
```@example
using TrajOptSOCPs; #hide
currSolveParams = solverParams(0.1, 0.5,
                                8, 10,
                                10^-4,
                                10, 10^6,
                                2.5, 2, 0.2, 0.2, 0.4)
TrajOptSOCPs.solParamPrint(currSolveParams)
```

For more information, see
```@docs
solverParams
solParamPrint
```

### Solve it!
At last, we can solve the problem!

```@example
# Initialize the primal-dual vector
# initTraj is from the rocket guess trajectory
# lambdaInit is from the dynamics constraints
initTrajPD = [initTraj; lambdaInit]

trajLambdaSolved, resArr = ALPrimalNewtonMain(initTrajPD, alRocket, currSolveParams)
```

For more information, see
```@docs
ALPrimalNewtonMain
```
