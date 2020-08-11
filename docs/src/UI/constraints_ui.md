# User Interface to Establish the Constraints

The constraints are at the heart of this project. So far, the package includes
three different types of constraints:
1. Dynamic Constraints (Affine Equality Constraints)
2. Ground Constraints (Affine Inequality Constraints)
3. Max Thrust Constraints (Second Order Cone Constraints)

We will explore each of these in this order.

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

The left hand side defines the matrix ``
