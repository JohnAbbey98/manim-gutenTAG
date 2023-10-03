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
    scaled_pdf_value = pdf_value * 20000

    return scaled_pdf_value


# animation for generateMetapeaks
class generateMetapeaks(Scene):
    '''
    Scene to observe how adjustments to the standard deviation of a normal
    distribution influences the shape of its probability density function.
    '''

    def construct(self):
        ax = Axes(
            x_range=[900, 920, 5],
            y_range=[0, 20000, 4000], 

            x_axis_config={
                "color": WHITE,
                "include_numbers": True
            },
            y_axis_config={
                "color": WHITE,
                "include_numbers": True 
            },
            axis_config={"color": WHITE}
        )

        # Add labels for x and y axis
        x_label = MathTex("m/z").next_to(ax.x_axis, DOWN)

        y_label = MathTex("Counts").next_to(ax.y_axis, UP)

        self.add(ax, x_label, y_label)

        # Initialize sigma and mu ValueTrackers
        mu = ValueTracker(910)
        sigma = ValueTracker(1)

        # Text to display distribution standard deviation and mean
        mu_text = MathTex(r'\mu = ').next_to(ax, UP, buff=0.2).set_color(YELLOW)
        sigma_text = MathTex(r'\sigma = ').next_to(ax, UP, buff=-0.2).set_color(YELLOW)
        
        # Always redraw the decimal value for mu and sigma for each frame
        mu_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2)
            .set_value(mu.get_value())
            .next_to(mu_text, RIGHT, buff=0.2)
            .set_color(YELLOW)
        )

        sigma_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2)
            .set_value(sigma.get_value())
            .next_to(sigma_text, RIGHT, buff=0.2)
            .set_color(YELLOW)
        )

        # Define PDF curve, always redraw for each frame
        curve_to_fill = always_redraw(
            lambda: ax.plot(
                lambda x: PDF_normal(x, mu.get_value(), sigma.get_value()), color=YELLOW)
        )

        curve_fill = curve_to_fill.copy().set_fill(
            color=YELLOW,  
            opacity=0.5    
        )

        # Define the vertical lines at multiple x-values
        x_values = [909.5, 910.5, 911.5]
        y_values = [4000, 3000, 2000]

        lines = VGroup()  # Create a vector group for the lines
        for (x, y) in zip(x_values, y_values):
            start_point = ax.c2p(x, 0)  # bottom of the line (x-value, y_min)
            end_point = ax.c2p(x, y)   # top of the line (x-value, y_max)
            line = DashedLine(start_point, end_point, color=BLUE, dash_length=0.1)
            lines.add(line)

        # shift rightwards to avoid margin issues
        everything_group = VGroup(ax, curve_to_fill, curve_fill, lines, x_label, y_label)  # Include all other objects that are part of the scene
        everything_group.shift(RIGHT*0.5)
        everything_group.shift(UP*0.3)


        # Start animation
        #self.add(ax, mu_text, mu_value_text, sigma_text, sigma_value_text, lines)
        self.add(ax, lines)

        # Apply curve
        self.play(Create(curve_to_fill), run_time = 3)

        # Apply the curve fill
        self.play(Write(curve_fill, lag_ratio=50, run_time = 2))

        # Apply wave
        self.play(ApplyWave(curve_to_fill, amplitude = 0.1), ApplyWave(curve_fill, amplitude = 0.05))



   
