# User Interface to Synthesize an Augmented Lagrangian

The last step is to form the Augmented Lagrangian. We generate the struct in
two steps (assuming that the objective and constraints have been instantiated).

```@example
# Equiped with the constraint term and the objective term, I now build the
# Augmented Lagrangian
penaltyStart = 1.0
alRocket = augLag(costFun, cMRocket, penaltyStart)
```

For more information, see
```@docs
augLag
```
