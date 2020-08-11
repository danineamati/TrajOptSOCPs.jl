# User Interface to Make a Model Rocket

```@meta
CurrentModule = TrajOptSOCPs
```

### Let's make a rocket!
For the sake of clarity, assume all units are in base SI units. The main simple struct is the `rocket_simple` struct:

```@docs
rocket_simple
```

Below is an example of how you make the most basic rocket. To obtain the best
results, we set the mass to 1 (hence all relevant units would be in "per
rocket mass" units). The gravity vector should match the problem dimensionality.
In the example below, the gravity vector is two dimensional and pointed in the
negative "y" direction.

```@example
mass = 1
isp = 1
grav = [0; -9.81]
deltaTime = 1
rocket = rocket_simple(mass, isp, grav, deltaTime)
```


### Guess a trajectory!

Now that we have a rocket, we specify the starting point (where the rocket is
right now) and the end point (where you want the rocket to land). Note that the
state `x` denotes the combined position and velocity vector: `[s; v]`.

In the example below, we start the rocket at `(2.0, 20.0)` and expect it to land
at `(0.0, 0.0)`. Moreover, we start the rocket in a freefall with velocity of
`(0.0, -15.0)` and expect it to softly land with zero velocity `(0.0, 0.0)`.
```@example
const rocketStart = [2.0; 20.0; 0.0; -15.0]
const rocketEnd = [0.0; 0.0; 0.0; 0.0]
nothing #hide
```

We can also initialize a base thrust, which is most simply the hover thrust.
```@example
mass = 1; grav = 1; #hide
uHover = mass * grav
nothing #hide
```

Lastly, we specify the number of steps in the simulation (combined with the `deltaTime` parameter above, this implies an a priori flight time).
```@example
# Number of time steps to discretize the trajectory
const NSteps = 60
nothing #hide
```

With the information above, we now have the information to build a guess
trajectory:
```@example
initTraj = initializeTraj(rocketStart, rocketEnd, uHover, uHover, NSteps)
```

For more details, see below:
```@docs
initializeTraj
```
