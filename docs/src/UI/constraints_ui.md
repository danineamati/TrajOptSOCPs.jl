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
