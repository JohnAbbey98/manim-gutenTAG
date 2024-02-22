# load libraries
from manim import *
from manim import RED
import math
import numpy as np

# Function for normal distribution PDF
def PDF_normal(x, mu, sigma):
    '''
    General form of probability density function of univariate normal distribution
    '''
    pdf_value = math.exp(-((x-mu)**2)/(2*sigma**2))/(sigma*math.sqrt(2*math.pi))

    # Scale the PDF value
    scaled_pdf_value = pdf_value * 25000

    return scaled_pdf_value


class metapeakFull(MovingCameraScene):
    def construct(self):
        self.camera.frame.save_state()

        # Define the axes
        ax = Axes(
            x_range=[850, 925, 50],
            y_range=[0, 15000, 5000],
            axis_config={"color": BLUE, "include_numbers": True}
        )

        # Add labels for x and y axis
        x_label = MathTex("m/z").next_to(ax.x_axis, DOWN)
        y_label = MathTex("Counts").next_to(ax.y_axis, UP)

        # Gaussian function using NumPy
        def gaussian(x, mu=910, sigma=1.5, scale = 12500):
            return scale * np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))

        # Generate x and y values
        x_vals = np.linspace(850, 925, 1000)  # Adjust range and number of points as needed
        y_vals = gaussian(x_vals, scale=12500)

        # Create the Gaussian curve
        curve = ax.plot_line_graph(x_vals, y_vals, add_vertex_dots=False, line_color=YELLOW, stroke_width=2)

        curve_fill = curve.copy().set_fill(
            color=YELLOW_B,  
            opacity=0.2    
        )

        # curve 2
        x_vals2 = np.linspace(850, 925, 1000)  # Adjust range and number of points as needed
        y_vals2 = gaussian(x_vals, mu = 872, scale=11000)
        curve2 = ax.plot_line_graph(x_vals2, y_vals2, add_vertex_dots=False, line_color=YELLOW, stroke_width=2)
        curve_fill2 = curve2.copy().set_fill(
            color=YELLOW_B,  
            opacity=0.2    
        )

        # curve 3 
        x_vals3 = np.linspace(850, 925, 1000)  # Adjust range and number of points as needed
        y_vals3 = gaussian(x_vals, mu = 886, sigma = 2, scale=25000)
        curve3 = ax.plot_line_graph(x_vals3, y_vals3, add_vertex_dots=False, line_color=YELLOW, stroke_width=2)
        curve_fill3 = curve3.copy().set_fill(
            color=YELLOW_B,  
            opacity=0.2    
        )

        # curve 4
        x_vals4 = np.linspace(850, 925, 1000)  # Adjust range and number of points as needed
        y_vals4 = gaussian(x_vals, mu = 856, sigma = 2, scale=20000)
        curve4 = ax.plot_line_graph(x_vals4, y_vals4, add_vertex_dots=False, line_color=YELLOW, stroke_width=2)
        curve_fill4 = curve4.copy().set_fill(
            color=YELLOW_B,  
            opacity=0.2    
        )


        # Define the vertical lines at multiple x-values
        x_spectrum_1 = [871.5, 872.5, 873.5]
        y_spectrum_1 = [2000, 1500, 1000]

        # Make spectrum 1
        spectrum_1 = VGroup()  # Create a vector group for the lines
        for (x, y) in zip(x_spectrum_1, y_spectrum_1):
            start_point = ax.c2p(x, 0)  # bottom of the line (x-value, y_min)
            end_point = ax.c2p(x, y)   # top of the line (x-value, y_max)
            line = DashedLine(start_point, end_point, color=RED_C, dash_length=0.07, stroke_width=2)
            spectrum_1.add(line)

        x_spectrum_2 = [909.5, 910.5, 911.5]
        y_spectrum_2 = [2250, 1750, 750]

        spectrum_2 = VGroup()  # Create a vector group for the lines
        for (x, y) in zip(x_spectrum_2, y_spectrum_2):
            start_point = ax.c2p(x, 0)  # bottom of the line (x-value, y_min)
            end_point = ax.c2p(x, y)   # top of the line (x-value, y_max)
            line = DashedLine(start_point, end_point, color=RED_C, dash_length=0.07, stroke_width=2)
            spectrum_2.add(line)

        x_spectrum_3 = [884.5, 885.5, 886.5, 887.5]
        y_spectrum_3 = [700, 4000, 3000, 1000]

        spectrum_3 = VGroup()  # Create a vector group for the lines
        for (x, y) in zip(x_spectrum_3, y_spectrum_3):
            start_point = ax.c2p(x, 0)  # bottom of the line (x-value, y_min)
            end_point = ax.c2p(x, y)   # top of the line (x-value, y_max)
            line = DashedLine(start_point, end_point, color=RED_C, dash_length=0.07, stroke_width=2)
            spectrum_3.add(line)

        x_spectrum_4 = [854.5, 855.5, 856.5, 857.5]
        y_spectrum_4 = [100, 3000, 2000, 1500]

        spectrum_4 = VGroup()  # Create a vector group for the lines
        for (x, y) in zip(x_spectrum_4, y_spectrum_4):
            start_point = ax.c2p(x, 0)  # bottom of the line (x-value, y_min)
            end_point = ax.c2p(x, y)   # top of the line (x-value, y_max)
            line = DashedLine(start_point, end_point, color=RED_C, dash_length=0.07, stroke_width=2)
            spectrum_4.add(line)



        # target line
        x_value_target = 908.5
        y_value_target = 4000
        target_line = VGroup()

        start_point = ax.c2p(x_value_target, 0)  # bottom of the line (x-value, y_min)
        end_point = ax.c2p(x_value_target, y_value_target)   # top of the line (x-value, y_max)
        #line = DashedLine(start_point, end_point, color=GREEN, dash_length=0.05)
        line = Line(start_point, end_point, color=PURPLE, stroke_width=3)
        target_line.add(line)

        # center line
        x_value_center = 910
        y_value_center = 4000
        center_line = VGroup()

        start_point = ax.c2p(x_value_center, 0)  # bottom of the line (x-value, y_min)
        end_point = ax.c2p(x_value_center, y_value_center)   # top of the line (x-value, y_max)
        #line = DashedLine(start_point, end_point, color=GREEN, dash_length=0.05)
        line = Line(start_point, end_point, color=GREEN_E, stroke_width=3)
        center_line.add(line)

        # Shift everything to avoid margin issues
        everything_group = VGroup(ax, x_label, curve, curve_fill, curve2, curve_fill2, curve3, curve_fill3, curve4, curve_fill4, y_label, target_line, center_line, spectrum_1, spectrum_2, spectrum_3, spectrum_4)  # Include all other objects that are part of the scene
        everything_group.shift(RIGHT * 0.5)

        labels = VGroup(x_label)
        labels.shift(UP*0.3)

        self.add(ax, x_label, y_label, spectrum_1, spectrum_2, spectrum_3, spectrum_4)

        self.wait(1)

        self.play(self.camera.frame.animate.scale(0.4).move_to(spectrum_2))
        

        # Add target line to scene
        self.play(Create(target_line), run_time = 0.5)

        self.wait(1)

         # Apply curve
        self.play(Create(curve), run_time = 3)

        # Apply the curve fill
        self.play(Write(curve_fill, lag_ratio=50, run_time = 1))

        self.play(Create(curve2), run_time = 0.1)
        self.play(Write(curve_fill2, lag_ratio=50, run_time = 0.1))

        self.play(Create(curve3), run_time = 0.5)
        self.play(Write(curve_fill3, lag_ratio=50, run_time = 0.1))

        self.play(Create(curve4), run_time = 0.5)
        self.play(Write(curve_fill4, lag_ratio=50, run_time = 0.1))

        # Apply wave
        self.play(ApplyWave(curve, amplitude = 0.1), 
                  ApplyWave(curve_fill, amplitude = 0.05))
        
        self.wait(1)
        
        # Add center line to scene
        self.play(Create(center_line), run_time = 0.5)

        self.wait(2)

        self.play(Restore(self.camera.frame))
        self.wait()