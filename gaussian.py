from manim import *
from manim import RED
import math
import numpy as np

RED = "#FF0000"

def PDF_normal(x, mu, sigma):
    '''
    General form of probability density function of univariate normal distribution
    '''
    return math.exp(-((x-mu)**2)/(2*sigma**2))/(sigma*math.sqrt(2*math.pi))


class AdjustMu(Scene):
    '''
    Scene to observe how adjustments to the mean of a normal distrubtion
    influences the shape of its probability density function
    '''

    def construct(self):
        ax = Axes(
            x_range = [700, 720, 5],
            y_range = [0, 0.5, 0.1],
            axis_config = {'include_numbers':True}
        )

        # Initialize mu (distribution mean) ValueTracker to 0
        mu = ValueTracker(710)

        # Text to display distrubtion mean
        mu_text = MathTex(r'\mu = ').next_to(ax, UP, buff=0.2).set_color(YELLOW)
        # Always redraw the decimal value for mu for each frame
        mu_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2)
            .set_value(mu.get_value())
            .next_to(mu_text, RIGHT, buff=0.2)
            .set_color(YELLOW)
        )

        # Define PDF curve, always redraw for each frame
        curve = always_redraw(
            lambda: ax.plot(
                lambda x: PDF_normal(x, mu.get_value(), 1), color=YELLOW)
                #lambda x: PDF_normal(x, mu.get_value(), 1), color=YELLOW)

        )

        # Start animation
        self.add(ax, mu_text, mu_value_text)

        self.play(Create(curve))
        self.play(
            mu.animate.set_value(710), run_time=2,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            mu.animate.set_value(705), run_time=1.5,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            mu.animate.set_value(715), run_time=2,
            rate_func=rate_functions.smooth
        )
        self.play(Uncreate(curve))


class AdjustSigma(Scene):
    '''
    Scene to observe how adjustments to the standard deviation of a normal
    distrubtion influences the shape of its probability density function
    '''

    def construct(self):
        ax = Axes(
            x_range = [-4, 4, 1],
            y_range = [0, 1, 0.2],
            axis_config = {'include_numbers':True}
        )

        # Initialize sigma (distribution standard deviation) ValueTracker to 1
        sigma = ValueTracker(1)

        # Text to display distrubtion standard deviation
        sigma_text = MathTex(r'\sigma = ').next_to(ax, UP, buff=0.2).set_color(YELLOW)
        # Always redraw the decimal value for sigma for each frame
        sigma_value_text = always_redraw(
            lambda: DecimalNumber(num_decimal_places=2)
            .set_value(sigma.get_value())
            .next_to(sigma_text, RIGHT, buff=0.2)
            .set_color(YELLOW)
        )

        # Define PDF curve, always redraw for each frame
        curve = always_redraw(
            lambda: ax.plot(
                lambda x: PDF_normal(x, 0, sigma.get_value()), color=YELLOW)
        )

        # Start animation
        self.add(ax, sigma_text, sigma_value_text)
        self.play(Create(curve))
        self.play(
            sigma.animate.set_value(1.5), run_time=1,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            sigma.animate.set_value(0.5), run_time=1.5,
            rate_func=rate_functions.smooth
        )
        self.wait()
        self.play(
            sigma.animate.set_value(1), run_time=1.25,
            rate_func=rate_functions.smooth
        )
        self.play(Uncreate(curve))



class AdjustMuSigma(Scene):
    '''
    Scene to observe how adjustments to the standard deviation of a normal
    distrubtion influences the shape of its probability density function
    '''

    def construct(self):
        ax = Axes(
            x_range = [700, 720, 5],
            y_range = [0, 1, 0.2],
            axis_config = {'include_numbers':True}
        )

        # Initialize sigma (distribution standard deviation) ValueTracker to 1
        # Initialize mu (distribution mean) ValueTracker to 0
        mu = ValueTracker(710)
        sigma = ValueTracker(1)

        # Text to display distrubtion standard deviation
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

        
        ######

                
        
        ######

        #Define PDF curve, always redraw for each frame
        curve = always_redraw(
            lambda: ax.plot(
                lambda x: PDF_normal(x, mu.get_value(), sigma.get_value()), color=YELLOW)
        )

        # Create a filled area under the curve using VGroup
        def get_filled_curve():
            filled_curve = ax.get_area(curve, color=YELLOW, opacity=0.3)
            return VGroup(curve, filled_curve)

        filled_curve_group = always_redraw(get_filled_curve)

        # Define the vertical lines at multiple x-values
        x_values = [709.5, 710.5, 711.5]
        y_values = [0.2, 0.13, 0.07]
        lines = VGroup()  # Create a vector group for the lines

        for (x, y) in zip(x_values, y_values):
            start_point = ax.c2p(x, 0)  # bottom of the line (x-value, y_min)
            end_point = ax.c2p(x, y)   # top of the line (x-value, y_max)
            line = DashedLine(start_point, end_point, color=BLUE, dash_length=0.1)
            lines.add(line)

        # Start animation
        self.add(ax, mu_text, mu_value_text, sigma_text, sigma_value_text, lines)
        self.play(Create(curve), run_time = 4)
        self.wait(0.5)  # Wait for 1 second
        # Apply the FadeInFromLeft effect to the fill.
        #self.play(FadeIn(filled_curve_group,  shift=LEFT))
        self.play(Transform(curve, filled_curve_group))  # Transform the curve into the filled version



class generateMetapeaks(Scene):
    '''
    Scene to observe how adjustments to the standard deviation of a normal
    distribution influences the shape of its probability density function.
    '''

    def construct(self):
        ax = Axes(
            x_range=[700, 720, 5],
            y_range=[0, 1, 0.2],
            axis_config={'include_numbers': True}
        )

        # Initialize sigma (distribution standard deviation) ValueTracker to 1
        # Initialize mu (distribution mean) ValueTracker to 0
        mu = ValueTracker(710)
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
            color=YELLOW,      # Starting color
            opacity=0.5      # Transparency
        )

        # Define the vertical lines at multiple x-values
        x_values = [709.5, 710.5, 711.5]
        y_values = [0.2, 0.13, 0.07]
        lines = VGroup()  # Create a vector group for the lines

        for (x, y) in zip(x_values, y_values):
            start_point = ax.c2p(x, 0)  # bottom of the line (x-value, y_min)
            end_point = ax.c2p(x, y)   # top of the line (x-value, y_max)
            line = DashedLine(start_point, end_point, color=BLUE, dash_length=0.1)
            lines.add(line)

        # Start animation
        self.add(ax, mu_text, mu_value_text, sigma_text, sigma_value_text, lines)
        self.play(Create(curve_to_fill), run_time = 3)

        # Wait for 1 second.
        #self.wait(0.5)

        # Apply the FadeInFromLeft effect to the fill.
        self.play(Write(curve_fill, lag_ratio=50, run_time = 2))
        #self.wait(0.1)
        self.play(ApplyWave(curve_to_fill, amplitude = 0.1), ApplyWave(curve_fill, amplitude = 0.05))



   
