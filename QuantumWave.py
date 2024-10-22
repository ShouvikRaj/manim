from manim import *
import numpy as np

class QuantumWaveFunction(Scene):
    def construct(self):
        # Parameters for the wave function
        sigma = 0.5 # Width of the Gaussian packet
        k0 = 2.0 # Initial momentum
        m = 1.0 # Mass
        hbar = 1.0 # Reduced Planck's constant
        
        # Create axes
        axes = Axes(
            x_range=[-5, 5, 1],
            y_range=[-1, 1, 0.5],
            axis_config={"color": BLUE},
            x_length=10,
            y_length=6
        )
        
        # Add labels
        x_label = axes.get_x_axis_label(r"x")
        y_label = axes.get_y_axis_label(r"|\psi(x,t)|")
        labels = VGroup(x_label, y_label)
        
        # Function to calculate wave function at time t
        def wave_function(x, t):
            # Complex Gaussian wave packet
            factor = 1 / (sigma * np.sqrt(2 * np.pi))
            x_spread = sigma**2 + (1j * hbar * t) / (2 * m)
            psi = factor * np.exp(-(x**2) / (4 * x_spread)) * np.exp(1j * k0 * x)
            return np.real(psi), np.imag(psi)
        
        # Create the initial wave function
        x_vals = np.linspace(-5, 5, 200)
        real_part, imag_part = wave_function(x_vals, 0)
        
        # Create the wave function plot
        wave_plot = axes.plot(
            lambda x: np.exp(-(x**2)/(4*sigma**2)) * np.cos(k0*x),
            color=YELLOW
        )
        
        # Add everything to the scene
        self.play(
            Create(axes),
            Write(labels)
        )
        self.play(Create(wave_plot))
        
        # Animate the evolution
        def update_wave(mob, dt):
            t = self.time
            real_part, imag_part = wave_function(x_vals, t)
            new_wave = axes.plot(
                lambda x: np.exp(-(x**2)/(4*(sigma**2 + (1j*hbar*t)/(2*m)))) * 
                         np.cos(k0*x - (k0**2*hbar*t)/(2*m)),
                color=YELLOW
            )
            mob.become(new_wave)
        
        # Add updater and let it run
        wave_plot.add_updater(update_wave)
        self.wait(6)
        
        # Clean up
        wave_plot.clear_updaters()
        self.wait(0.5)

if __name__ == "__main__":
    config.pixel_height = 720
    config.pixel_width = 1280
    config.frame_rate = 30
    scene = QuantumWaveFunction()
    scene.render()

class QuantumWave3D(ThreeDScene):
    def construct(self):
        # Parameters
        sigma = 0.5 # Width of the Gaussian packet
        k0 = 2.0 # Initial momentum
        m = 1.0 # Mass
        hbar = 1.0 # Reduced Planck's constant
        omega = 2.0 # Angular frequency for spiral motion

        # Set up the 3D axes
        axes = ThreeDAxes(
            x_range=[-5, 5, 1],
            y_range=[-2, 2, 1],
            z_range=[-2, 2, 1],
            x_length=10,
            y_length=4,
            z_length=4,
        )

        # Add labels
        x_label = axes.get_x_axis_label(r"x")
        y_label = axes.get_y_axis_label(r"Re(\psi)")
        z_label = axes.get_z_axis_label(r"Im(\psi)")
        labels = VGroup(x_label, y_label, z_label)

        # Create the initial wave function data
        x_vals = np.linspace(-5, 5, 100)
        
        # Function to calculate wave function at time t
        def wave_function(x, t):
            factor = 1 / (sigma * np.sqrt(2 * np.pi))
            x_spread = sigma**2 + (1j * hbar * t) / (2 * m)
            psi = factor * np.exp(-(x**2) / (4 * x_spread)) * np.exp(1j * k0 * x)
            # Add spiral motion
            spiral_factor = np.exp(1j * omega * t)
            return psi * spiral_factor

        # Create initial 3D curve
        def get_wave_points(t):
            psi = wave_function(x_vals, t)
            return np.array([[x, psi.real, psi.imag] for x, psi in zip(x_vals, psi)])

        wave_curve = ParametricFunction(
            lambda x: np.array([
                x,
                np.exp(-(x**2)/(4*sigma**2)) * np.cos(k0*x),
                np.exp(-(x**2)/(4*sigma**2)) * np.sin(k0*x)
            ]),
            t_range=[-5, 5],
            color=YELLOW
        )

        # Set up the scene
        self.set_camera_orientation(phi=70 * DEGREES, theta=30 * DEGREES)
        self.begin_ambient_camera_rotation(rate=0.1)

        # Add everything to the scene
        self.play(
            Create(axes),
            Write(labels)
        )
        self.play(Create(wave_curve))

        # Create a value tracker for time
        time_tracker = ValueTracker(0)

        # Update function for the wave
        def update_wave(mob):
            t = time_tracker.get_value()
            new_points = get_wave_points(t)
            mob.set_points_smoothly([
                axes.coords_to_point(*point) for point in new_points
            ])

        # Add updater and animate
        wave_curve.add_updater(update_wave)
        self.play(
            time_tracker.animate.set_value(6),
            run_time=6,
            rate_func=linear
        )

        # Clean up
        wave_curve.clear_updaters()
        self.wait(0.5)

if __name__ == "__main__":
    config.pixel_height = 720
    config.pixel_width = 1280
    config.frame_rate = 30
    scene = QuantumWave3D()
    scene.render()
