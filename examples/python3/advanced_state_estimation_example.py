#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Advanced State Estimation Example

This example demonstrates advanced state estimation techniques for Cosserat beams
using Lie groups. It includes:

1. Extended Kalman Filter (EKF) for pose estimation
2. Uncertainty propagation using Lie algebra
3. Multi-modal sensor fusion
4. Real-time state estimation for beam control
"""

import numpy as np
import math
from Cosserat.LieGroups import SO3, SE3, PoseVel
import matplotlib.pyplot as plt

def print_section(title):
    """Helper to print section titles with formatting"""
    print("\n" + "="*80)
    print(f" {title} ".center(80, "-"))
    print("="*80)

class LieGroupEKF:
    """Extended Kalman Filter using Lie groups for beam state estimation"""

    def __init__(self, initial_pose, initial_covariance):
        """Initialize EKF with Lie group state representation"""
        self.state = initial_pose  # SE3 pose
        self.covariance = initial_covariance  # 6x6 covariance matrix
        self.process_noise = np.eye(6) * 0.01  # Process noise
        self.measurement_noise = np.eye(6) * 0.1  # Measurement noise

    def predict(self, dt, control_input=None):
        """Prediction step using Lie group integration"""
        # For simplicity, assume constant velocity model
        if control_input is not None:
            # Control input is in Lie algebra (twist)
            delta_pose = SE3.exp(control_input * dt)
            self.state = self.state * delta_pose
        else:
            # No control, assume constant pose
            pass

        # Propagate covariance using Lie group Jacobians
        # Simplified: assume identity Jacobian for demonstration
        F = np.eye(6)  # State transition Jacobian
        self.covariance = F @ self.covariance @ F.T + self.process_noise

    def update(self, measurement):
        """Update step with measurement in Lie algebra"""
        # Compute innovation (error in Lie algebra)
        innovation = (measurement * self.state.inverse()).log()

        # Measurement Jacobian (simplified)
        H = np.eye(6)

        # Kalman gain
        S = H @ self.covariance @ H.T + self.measurement_noise
        K = self.covariance @ H.T @ np.linalg.inv(S)

        # Update state using Lie group composition
        correction = SE3.exp(K @ innovation.log())
        self.state = correction * self.state

        # Update covariance
        I = np.eye(6)
        self.covariance = (I - K @ H) @ self.covariance

    def get_state(self):
        return self.state

    def get_covariance(self):
        return self.covariance

class BeamStateEstimator:
    """Advanced beam state estimator using multiple sensors"""

    def __init__(self, beam_geometry):
        self.beam_geometry = beam_geometry
        self.ekf = None
        self.initialize_estimator()

    def initialize_estimator(self):
        """Initialize the state estimator"""
        # Initial pose (straight beam)
        initial_pose = SE3()

        # Initial uncertainty
        initial_covariance = np.eye(6) * 0.1

        self.ekf = LieGroupEKF(initial_pose, initial_covariance)

    def estimate_state(self, sensor_measurements, dt, control_input=None):
        """Perform state estimation using multiple sensors"""
        # Prediction step
        self.ekf.predict(dt, control_input)

        # Fuse measurements from multiple sensors
        fused_measurement = self.fuse_measurements(sensor_measurements)

        # Update step
        if fused_measurement is not None:
            self.ekf.update(fused_measurement)

        return self.ekf.get_state(), self.ekf.get_covariance()

    def fuse_measurements(self, measurements):
        """Fuse measurements from multiple sensors"""
        if not measurements:
            return None

        # Simple averaging for demonstration
        # In practice, would use more sophisticated fusion
        poses = [m['pose'] for m in measurements if 'pose' in m]

        if not poses:
            return None

        # Average poses (simplified)
        avg_translation = np.mean([p.translation for p in poses], axis=0)
        # For rotations, would need proper averaging on SO(3)

        return SE3(SO3(), avg_translation)  # Simplified

class MultiSensorSimulator:
    """Simulate multiple sensors for beam state estimation"""

    def __init__(self, true_trajectory):
        self.true_trajectory = true_trajectory
        self.sensors = []

    def add_sensor(self, sensor_type, noise_std, bias=None):
        """Add a sensor to the simulation"""
        sensor = {
            'type': sensor_type,
            'noise_std': noise_std,
            'bias': bias or np.zeros(6)
        }
        self.sensors.append(sensor)

    def get_measurements(self, time_index):
        """Get noisy measurements from all sensors at given time"""
        true_pose = self.true_trajectory[time_index]
        measurements = []

        for sensor in self.sensors:
            # Add noise to true pose
            noise = np.random.normal(0, sensor['noise_std'], 6)
            noisy_twist = true_pose.log() + noise + sensor['bias']

            # Convert back to SE3
            measured_pose = SE3.exp(noisy_twist)

            measurements.append({
                'sensor_type': sensor['type'],
                'pose': measured_pose,
                'timestamp': time_index
            })

        return measurements

def simulate_beam_estimation():
    """Simulate beam state estimation with multiple sensors"""

    print_section("Beam State Estimation Simulation")

    # Create true beam trajectory (simulated deformation)
    num_steps = 100
    dt = 0.01

    true_trajectory = []
    for i in range(num_steps):
        # Simulate oscillating beam deformation
        angle = 0.2 * math.sin(2 * math.pi * i * dt)
        translation = np.array([20.0, 0.0, 0.0])

        pose = SE3(SO3(angle, np.array([0.0, 1.0, 0.0])), translation)
        true_trajectory.append(pose)

    # Set up multi-sensor simulation
    sensor_sim = MultiSensorSimulator(true_trajectory)

    # Add different types of sensors
    sensor_sim.add_sensor('imu', noise_std=0.05)  # Inertial measurement unit
    sensor_sim.add_sensor('camera', noise_std=0.1)  # Vision system
    sensor_sim.add_sensor('strain_gauge', noise_std=0.02)  # Strain sensors

    # Initialize state estimator
    estimator = BeamStateEstimator(None)  # Simplified geometry

    # Run estimation
    estimated_trajectory = []
    covariances = []

    for i in range(num_steps):
        # Get sensor measurements
        measurements = sensor_sim.get_measurements(i)

        # Estimate state
        estimated_pose, covariance = estimator.estimate_state(measurements, dt)

        estimated_trajectory.append(estimated_pose)
        covariances.append(covariance)

    # Compute estimation errors
    errors = []
    for i in range(num_steps):
        error = (true_trajectory[i] * estimated_trajectory[i].inverse()).log()
        errors.append(np.linalg.norm(error))

    print(f"Estimation completed over {num_steps} time steps")
    print(".4f"
    print(".4f"
    print(".4f"

    # Plot results
    plot_estimation_results(true_trajectory, estimated_trajectory, errors, covariances)

def plot_estimation_results(true_trajectory, estimated_trajectory, errors, covariances):
    """Plot state estimation results"""

    # Extract data for plotting
    true_positions = [p.translation for p in true_trajectory]
    est_positions = [p.translation for p in estimated_trajectory]

    true_angles = [p.rotation.angle for p in true_trajectory]
    est_angles = [p.rotation.angle for p in estimated_trajectory]

    # Create plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # Position tracking
    ax1.plot([p[0] for p in true_positions], label='True X')
    ax1.plot([p[0] for p in est_positions], '--', label='Estimated X')
    ax1.set_xlabel('Time step')
    ax1.set_ylabel('X Position')
    ax1.legend()
    ax1.set_title('Position Tracking (X)')

    # Orientation tracking
    ax2.plot(true_angles, label='True Angle')
    ax2.plot(est_angles, '--', label='Estimated Angle')
    ax2.set_xlabel('Time step')
    ax2.set_ylabel('Angle (rad)')
    ax2.legend()
    ax2.set_title('Orientation Tracking')

    # Estimation error
    ax3.plot(errors)
    ax3.set_xlabel('Time step')
    ax3.set_ylabel('Estimation Error')
    ax3.set_title('Estimation Error Over Time')
    ax3.grid(True)

    # Uncertainty evolution
    uncertainties = [np.trace(cov) for cov in covariances]
    ax4.plot(uncertainties)
    ax4.set_xlabel('Time step')
    ax4.set_ylabel('Trace of Covariance')
    ax4.set_title('Estimation Uncertainty')
    ax4.grid(True)

    plt.tight_layout()
    plt.savefig('beam_state_estimation.png', dpi=150, bbox_inches='tight')
    print("Results saved to 'beam_state_estimation.png'")

def demonstrate_uncertainty_propagation():
    """Demonstrate uncertainty propagation using Lie groups"""

    print_section("Uncertainty Propagation with Lie Groups")

    # Initial pose and uncertainty
    initial_pose = SE3()
    initial_covariance = np.diag([0.01, 0.01, 0.01, 0.001, 0.001, 0.001])

    print("Initial pose:", initial_pose.translation)
    print("Initial uncertainty (diagonal):", np.diag(initial_covariance))

    # Simulate motion with control input
    dt = 0.1
    control_input = np.array([0.1, 0.0, 0.0, 0.0, 0.0, 0.1])  # Linear + angular velocity

    # Propagate uncertainty through multiple time steps
    current_pose = initial_pose
    current_covariance = initial_covariance

    poses = [current_pose]
    covariances = [current_covariance]

    for i in range(10):
        # Apply motion
        delta_pose = SE3.exp(control_input * dt)
        current_pose = current_pose * delta_pose

        # Propagate covariance (simplified - would use proper Lie group Jacobians)
        # In practice, this requires computing the state transition Jacobian
        F = np.eye(6)  # Simplified
        process_noise = np.eye(6) * 0.001
        current_covariance = F @ current_covariance @ F.T + process_noise

        poses.append(current_pose)
        covariances.append(current_covariance)

    print("Final pose after 10 steps:", current_pose.translation)
    print("Final uncertainty (diagonal):", np.diag(current_covariance))
    print(".3f"

def main():
    """Main function demonstrating advanced state estimation"""

    print("Advanced State Estimation for Cosserat Beams")
    print("="*80)

    # Demonstrate uncertainty propagation
    demonstrate_uncertainty_propagation()

    # Run full state estimation simulation
    simulate_beam_estimation()

    print("\n" + "="*80)
    print("Advanced state estimation demonstration completed!")
    print("Features demonstrated:")
    print("- Extended Kalman Filter using Lie groups")
    print("- Multi-sensor fusion")
    print("- Uncertainty propagation")
    print("- Real-time state estimation")
    print("="*80)

if __name__ == "__main__":
    main()</content>
</xai:function_call">Now let me create an example for multi-section beam support with section-specific properties. This will show how to handle beams with different properties in each section. 

<xai:function_call name="write">
<parameter name="filePath">examples/python3/multi_section_beam_example.py