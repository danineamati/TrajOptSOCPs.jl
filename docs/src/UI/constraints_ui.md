# User Interface to Establish the Constraints

The constraints are at the heart of this project. So far, the package includes
three different types of constraints:
1. Dynamic Constraints (Affine Equality Constraints)
2. Ground Constraints (Affine Inequality Constraints)
3. Max Thrust Constraints (Second Order Cone Constraints)

We will explore each of these in this order. All of the constraints are fed
into a `constraintManager`.

### Dynamic Constraints
In continuous time, we have a differential equation of the form:

```math
\frac{dx}{dt} = \tilde Ax + \tilde Bu + \tilde C
```

We solve this via Matrix exponentiation which yields a discretized equation
of the form:

```math
x_{k+1} = A x_k + B u_k + C
```

Where ``\tilde A \neq A``, ``\tilde B \neq B``, and ``\tilde C \neq C``. To
formulate this as an Affine Equality constraint, we want to write this in the
form ``Ax = b``. Thus we rearrange this as

```math
A x_k + B u_k - I x_{k+1} = -C
```

The left hand side defines the matrix ``A_{DYN}`` and the right hand side
defines ``b_{DYN}``. Generally, we also want to include constraints on the
initial and final states specifically. This is all handled together in the
below call:

```@example
# Create the Dynamics Constraints
const ADyn, BDyn = rocketDynamicsFull(rocket, rocketStart, rocketEnd, NSteps)
```

Then the dynamics constraint and associated dual variable are created:
```@example
dynConstraint = AL_AffineEquality(ADyn, BDyn)
lambdaInit = -1 * ones(size(BDyn))
```

Note that the two are separated just in case the you (1) don't care about the
initial or final points or (2) want the ``A_{DYN}`` matrix and ``b_{DYN}``
vector after the fact.

### Ground Constraints
Currently the ground constraints only handle flat ground level with the origin.
This is generally acceptable for the rocket soft landing problem since we can
choose the coordinate system and rockets generally land on flat areas. However,
it is totally possible to include more complicated terrain in future updates.

Since this is a simple constraint, `makeGroundConstraint` handles it in one go.
We only need to provide the dimensionality and normal vector to the ground. As
before, we initialize the associated dual variables.

```@example
# Create the Ground Constraints
groundConstraint = makeGroundConstraint(NSteps, size(grav, 1), size(grav, 1))
groundLambda = zeros(size(groundConstraint.b, 1))
```

### Max Thrust Constraints
This is the crux of the novel addition enabled by this package. Most of the
work is handled behind the scenes. The user only needs to provide the maximum
thrust (as a positive `Float64`). As before, initialize the associated dual
variables.

```@example
thrustMax = 20.0
maxThrustConstraint = makeMaxThrustConstraint(NSteps, size(grav, 1), thrustMax)
maxThrustLambda = zeros(size(maxThrustConstraint.indicatorList))
```

### Stitch it all together!
Now we stitch it all together. The dynamics are put as a special case, the rest
of the constraints are added as an array:

```@example
cMRocket = constraintManager_Dynamics([groundConstraint, maxThrustConstraint],
                                      [groundLambda, maxThrustLambda],
                                      dynConstraint, lambdaInit)
```

For more information see:
```@docs
constraintManager_Dynamics
```
