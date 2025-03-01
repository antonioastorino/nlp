'''
Attempt to implement a non-linear pendulum with perturbations using:

    alpha = - g / l * sin (theta) + (input_force + damping_force) / m / l

This scripts generates and plots the behavior of `theta` with respect to time and a rudimental
animation of the pendulum.

The initial condition can be assigned by setting `forced_torque_vec` as a vector of torque applied
during the first N time steps

'''

import time
import numpy as np
import matplotlib.pyplot as plt

g = 9.81  # Gravitational acceleration [m / s^2]
dt = 0.001  # Time step [s]

######################################## Parameters ###############################################
animation = True
DURATION = 10  # Simulation duration [s]
l = 0.5  # Pendulum length [m]
m = 0.5  # [kg]
damping_k = 0.8  # Damping coefficient [kg / s^2]
theta_d = np.pi / 2  # Desired angle [rad]
MAX_TORQUE = 50  # Max controlling torque [Nm]
K_p = 100
K_i = 45
K_d = 15 
disturbance = [0 for _ in range(1000)] + [-100 for _ in range(10)]
#################################### End of parameters ############################################

NUM_OF_SAMPLES = round(DURATION / dt)


def myPlotter(xVals, yVals, xLabel, yLabel, title):
    plt.figure(figsize=(8, 6))
    plt.plot(xVals, yVals, color='blue', linewidth=4)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.grid(visible=True)


def theta_t_plus_dt(dt, theta_t, theta_t_minus_dt, forced_torque=0, damping_k=0):
    alpha_f = damping_k * (theta_t_minus_dt - theta_t) / dt / m
    alpha_in = forced_torque / l / l / m
    alpha_g = -g / l * np.sin(theta_t)
    alpha = alpha_g + alpha_in + alpha_f
    return alpha * (dt * dt) - theta_t_minus_dt + 2 * theta_t


t = np.linspace(0, NUM_OF_SAMPLES * dt, NUM_OF_SAMPLES + 1)
theta_curr = 0
theta_prev = 0
theta_next = 0
theta_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
omega_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
torque_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
theta_err = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
disturbance_vec = [0 for _ in range(0, NUM_OF_SAMPLES + 1)]
disturbance_vec[0:min(NUM_OF_SAMPLES + 1, len(disturbance))] = disturbance[0:]
err_integral = 0
err_derivative = 0

start_time = time.time()
for i in range(0, NUM_OF_SAMPLES + 1):
    theta_err[i] = (theta_d - theta_curr)
    err_integral = err_integral + theta_err[i]
    if i > 0:
        err_derivative = theta_err[i] - theta_err[i - 1]
    ctrl_torque = K_p * theta_err[i] + K_i * (err_integral * dt) + K_d * (err_derivative / dt)
    ctrl_torque = np.sign(ctrl_torque) * min(abs(ctrl_torque), MAX_TORQUE)
    theta_next = theta_t_plus_dt(
        dt,
        theta_curr,
        theta_prev,
        ctrl_torque +
        disturbance_vec[i],
        damping_k)
    omega_vec[i] = (theta_next - theta_curr) / dt
    theta_prev = theta_curr
    theta_curr = theta_next
    theta_vec[i] = theta_curr
    torque_vec[i] = ctrl_torque

print("Computation time =", time.time() - start_time)
print("Final error % =", (theta_d - theta_vec[-1]) / theta_d * 100)

theta_plot_vec = [theta_vec[i] / np.pi for i in range(0, NUM_OF_SAMPLES + 1)]
err_plot_vec = [theta_err[i] / theta_d * 100 for i in range(0, NUM_OF_SAMPLES + 1)]
myPlotter(t, theta_plot_vec, 't [s]', r"$\frac{\theta}{\pi}$ [rad]", 'angle vs time')
myPlotter(t, err_plot_vec, 't [s]', r"$(\theta_d - \theta)/\theta_d * 100 $ [%]", 'error vs time')
myPlotter(t, torque_vec, 't [s]', 'applied torque [Nm]', 'Applied torque')

# Energy / state space
myPlotter(theta_vec, omega_vec, 'theta [rad]', 'omega [rad/s]', 'State Space')
KE = [(omega_vec[i] * l)**2 * m / 2 for i in range(0, NUM_OF_SAMPLES + 1)]
U = [m * g * l * (1 - np.cos(theta_vec[i])) for i in range(0, NUM_OF_SAMPLES + 1)]
E = [KE[i] + U[i] for i in range(0, NUM_OF_SAMPLES + 1)]
myPlotter(KE, U, 'kinetic [J]', 'potential [J]', 'Energy')
myPlotter(t, E, 'time [s]', 'total energy [J]', 'Energy vs time')

if (animation):
    fig, ax = plt.subplots()
    ax.set_xlim(-l * 1.1, l * 1.1)
    ax.set_ylim(-l * 1.1, l * 1.1)
    ax.set_aspect('equal', adjustable='box')
    position_plot, = ax.plot(0, 0, marker='*')
    refresh_rate = 50  # Hz
    undersampling_rate = round(1 / dt / refresh_rate)
    for i in range(0, NUM_OF_SAMPLES + 1, undersampling_rate):
        position_plot.set_data([0, l * np.sin(theta_vec[i])], [0, - l * np.cos(theta_vec[i])])
        plt.pause(1 / refresh_rate)

plt.show()
