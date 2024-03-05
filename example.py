import numpy as np
import matplotlib.pyplot as plt

def derivs(t, state):
    """
    This function defines the differential equations for the double pendulum.
    t: time variable (not used here, but needed for odeint)
    state: array containing the angles and angular velocities of the two pendulum rods
    """
    # Unpack the state vector
    theta1, theta2, omega1, omega2 = state
    
    # Parameters
    L1 = 1.0  # length of pendulum 1
    L2 = 1.0  # length of pendulum 2
    m1 = 1.0  # mass of pendulum 1
    m2 = 1.0  # mass of pendulum 2
    g = 9.8   # acceleration due to gravity
    
    # Equations of motion
    dydt = [omega1,
            omega2,
            (-g * (2 * m1 + m2) * np.sin(theta1) - m2 * g * np.sin(theta1 - 2 * theta2) - 2 * np.sin(theta1 - theta2) * m2 * (omega2 ** 2 * L2 + omega1 ** 2 * L1 * np.cos(theta1 - theta2))) /
            (L1 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2))),
            (2 * np.sin(theta1 - theta2) * (omega1 ** 2 * L1 * (m1 + m2) + g * (m1 + m2) * np.cos(theta1) + omega2 ** 2 * L2 * m2 * np.cos(theta1 - theta2))) /
            (L2 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))]
    return dydt

def runge_kutta_step(f, t, y, h):
    """
    This function performs one step of the fourth-order Runge-Kutta method.
    f: function representing the system of differential equations
    t: current time
    y: current state
    h: step size
    """
    k1 = h * np.array(f(t, y))
    k2 = h * np.array(f(t + h / 2, y + k1 / 2))
    k3 = h * np.array(f(t + h / 2, y + k2 / 2))
    k4 = h * np.array(f(t + h, y + k3))
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

def simulate_double_pendulum(theta1_init, theta2_init, omega1_init, omega2_init, simulation_time, time_step):
    """
    Simulate the double pendulum motion using the fourth-order Runge-Kutta method.
    theta1_init: initial angle of pendulum 1
    theta2_init: initial angle of pendulum 2
    omega1_init: initial angular velocity of pendulum 1
    omega2_init: initial angular velocity of pendulum 2
    simulation_time: total time to simulate
    time_step: time step for the simulation
    """
    # Initialize arrays to store results
    times = np.arange(0, simulation_time, time_step)
    theta1_values = np.zeros(len(times))
    theta2_values = np.zeros(len(times))
    omega1_values = np.zeros(len(times))
    omega2_values = np.zeros(len(times))
    
    # Set initial conditions
    state = np.array([theta1_init, theta2_init, omega1_init, omega2_init])
    
    # Simulate the motion
    for i, t in enumerate(times):
        theta1_values[i] = state[0]
        theta2_values[i] = state[1]
        omega1_values[i] = state[2]
        omega2_values[i] = state[3]
        state = runge_kutta_step(derivs, t, state, time_step)
    
    # Return the results
    return times, theta1_values, theta2_values, omega1_values, omega2_values

# Initial conditions
theta1_init = np.pi / 4
theta2_init = np.pi / 2
omega1_init = 0
omega2_init = 0

# Simulation parameters
simulation_time = 10
time_step = 0.01

# Simulate the double pendulum motion
times, theta1_values, theta2_values, omega1_values, omega2_values = simulate_double_pendulum(theta1_init, theta2_init, omega1_init, omega2_init, simulation_time, time_step)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(times, theta1_values, label='Theta 1')
plt.plot(times, theta2_values, label='Theta 2')
plt.title('Double Pendulum Simulation')
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.legend()
plt.grid(True)
plt.show()

# Print the final conditions
final_index = -1  # index of the final time step
final_conditions = (theta1_values[final_index], theta2_values[final_index], omega1_values[final_index], omega2_values[final_index])
print("Final conditions: ", final_conditions)