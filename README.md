# Non-Linear Pendulum
Simulate a non-linear pendulum in time domain without using the small-angle approximation.

Let:
- `theta` be the angular position
- `omega` be the angular velocity
- `alpha` be the angular acceleration

An immediate discretization method would look as follows:
- `theta(t + dt) = omega(t) * dt`
- `omaga(t + dt) = alpha(t) * dt`

However, I was initially unable to calculate omega, and hence I found an alternative solution.

## Method 1 - Taylor expansion
The discretized equation for the `theta` can be derived as follows:
- Compute two Taylor expansions for `theta(t +/- dt)` truncated at the 2nd order wherein
- I add the two expansions and solve for `theta(t + dt)`

```
theta(t + dt) - theta(t) =~ +omega(t) * dt + 1/2 alpha(t) * dt^2  +
theta(t - dt) - theta(t) =~ -omega(t) * dt + 1/2 alpha(t) * dt^2  =
--------------------------------------------------------------------------------
theta(t + dt) + theta(t - dt) =~ alpha(t) * dt^2 =>
theta(t + dt) =~ alpha(t) * dt^2 - theta(t - dt)
```

> NOTE 1: The result does not involve `omega`

> NOTE 2: Expanding to the 3rd order would produce the same result -> this approximation is
          very good!

I don't trust this approach because I made it up. However, I verified that:
- The angular frequency for small angles
- The energy is conserved when no perturbation is applied

### Example 1 - impulse response - no perturbations
This demo shows the validity of this method when no perturbation is applied and the energy is conserved.
```
python3 non-linear-pendulum.py
```
With the default values, the output is
![./doc/angle_no-damping.png](doc/angle_no-damping.png)


### Example 2 - damping force + input force 
In this demo, external forces are applied and the result looks ok.
```
python3 non-linear-pendulum-with-perturbations.py
```

![./doc/input_damped-plus-sinusoidal-input.png](doc/input_damped-plus-sinusoidal-input.png)
![./doc/angle_damped-plus-sinusoidal-input.png](doc/angle_damped-plus-sinusoidal-input.png)
![./doc/state_damped-plus-sinusoidal-input.png](doc/state_damped-plus-sinusoidal-input.png)

> WARNING! Don't stare at the pendulum for too long! :)


## Method 2 - Calculate `omega`
I eventually figured out that `omega` can be immediately derived from the conservation of energy.
```
m * v^2
------- - m * g * l * (1 - cos(theta)) = E
   2
```

Problems:
- Only the absolute value of `omega` can be calculated from the above equation. Some logic to detect a direction change is needed.
- The total energy deviation does not improve compared to the previous method.

For now, the velocity sign is inverted when the kinetic energy is almost 0 and is less than the kinetic energy at the previous step.

For testing, run
```
python3 non-linear-pendulum-explicit-omega.py
```

# Resources:
- [Equations](https://en.wikipedia.org/wiki/Pendulum_(mechanics))
