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
                "color": BLUE,
                "include_numbers": True
            },
            y_axis_config={
                "color": BLUE,
                "include_numbers": True 
            },
            axis_config={"color": BLUE}
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
            line = DashedLine(start_point, end_point, color=RED_B, dash_length=0.1)
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



# animation for generateMetapeaks
class assignMetapeaks(Scene):
    '''
    Scene to observe how adjustments to the standard deviation of a normal
    distribution influences the shape of its probability density function.
    '''

    def construct(self):
        ax = Axes(
            x_range=[900, 920, 5],
            y_range=[0, 20000, 4000], 

            x_axis_config={
                "color": BLUE,
                "include_numbers": True
            },
            y_axis_config={
                "color": BLUE,
                "include_numbers": True 
            },
            axis_config={"color": BLUE}
        )

        # Add labels for x and y axis
        x_label = MathTex("m/z").next_to(ax.x_axis, DOWN)
        y_label = MathTex("Counts").next_to(ax.y_axis, UP)


        #self.add(ax, x_label, y_label)

        # Initialize sigma and mu ValueTrackers
        mu = ValueTracker(910)
        sigma = ValueTracker(1)

        # Define PDF curve, always redraw for each frame
        curve_to_fill = always_redraw(
            lambda: ax.plot(
                lambda x: PDF_normal(x, mu.get_value(), sigma.get_value()), color=YELLOW)
        )

        curve_fill = curve_to_fill.copy().set_fill(
            color=YELLOW_B,  
            opacity=0.2    
        )

        # Define vertical line for target m/z value
        x_value_target = 909
        y_value_target = 8000
        target_line = VGroup()

        start_point = ax.c2p(x_value_target, 0)  # bottom of the line (x-value, y_min)
        end_point = ax.c2p(x_value_target, y_value_target)   # top of the line (x-value, y_max)
        #line = DashedLine(start_point, end_point, color=GREEN, dash_length=0.05)
        line = Line(start_point, end_point, color=PURPLE, stroke_width=5)
        target_line.add(line)


        # Define the vertical lines at multiple x-values
        x_values = [909.5, 910.5, 911.5]
        y_values = [4000, 3000, 2000]

        lines = VGroup()  # Create a vector group for the lines
        for (x, y) in zip(x_values, y_values):
            start_point = ax.c2p(x, 0)  # bottom of the line (x-value, y_min)
            end_point = ax.c2p(x, y)   # top of the line (x-value, y_max)
            line = DashedLine(start_point, end_point, color=RED_C, dash_length=0.1)
            lines.add(line)

        # Define vertical line for metapeak center
        x_value_green = 910
        y_value_green = 8000
        green_line = VGroup()
        
        start_point = ax.c2p(x_value_green, 0)  # bottom of the line (x-value, y_min)
        end_point = ax.c2p(x_value_green, y_value_green)   # top of the line (x-value, y_max)
        #line = DashedLine(start_point, end_point, color=GREEN, dash_length=0.05)
        line = Line(start_point, end_point, color=GREEN_E, stroke_width=5)
        green_line.add(line)

        # shift rightwards to avoid margin issues
        everything_group = VGroup(ax, curve_to_fill, curve_fill, lines, x_label, y_label, green_line, target_line)  # Include all other objects that are part of the scene
        everything_group.shift(RIGHT*0.5)
        everything_group.shift(UP*0.3)


        # Start animation
        #self.add(ax, mu_text, mu_value_text, sigma_text, sigma_value_text, lines)
        self.add(ax, lines, x_label, y_label, target_line)

        # Apply curve
        self.play(Create(curve_to_fill), run_time = 3)

        # Apply the curve fill
        self.play(Write(curve_fill, lag_ratio=50, run_time = 1))

        # Apply wave
        self.play(ApplyWave(curve_to_fill, amplitude = 0.1), ApplyWave(curve_fill, amplitude = 0.05))

        # Add green line to scene
        self.play(Create(green_line), run_time = 0.5)




class ZoomTest(MovingCameraScene):
    def construct(self):
        self.camera.frame.save_state()

        # define stage (axes)
        ax = Axes(x_range=[850, 925, 50], y_range=[0, 25000, 5000],            
            axis_config={"color": BLUE, "include_numbers": True})
        
        # Add labels for x and y axis
        x_label = MathTex("m/z").next_to(ax.x_axis, DOWN)
        y_label = MathTex("Counts").next_to(ax.y_axis, UP)
        
        # Initialize sigma and mu ValueTrackers
        mu = ValueTracker(910)
        sigma = ValueTracker(1)

        # Define PDF curve, always redraw for each frame
        curve_to_fill = always_redraw(
            lambda: ax.plot(
                lambda x: PDF_normal(x, mu.get_value(), sigma.get_value()), color=YELLOW)
        )

        curve_fill = curve_to_fill.copy().set_fill(
            color=YELLOW_B,  
            opacity=0.2    
        )


        #graph = ax.plot(lambda x: np.sin(x), color=BLUE, x_range=[850, 950])
        #graph = ax.plot(lambda x: curve_to_fill, color=YELLOW, x_range=[800, 1000 * PI])

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

        x_spectrum_3 = [939.5, 940.5]
        y_spectrum_3 = [2000, 1500]

        spectrum_3 = VGroup()  # Create a vector group for the lines
        for (x, y) in zip(x_spectrum_3, y_spectrum_3):
            start_point = ax.c2p(x, 0)  # bottom of the line (x-value, y_min)
            end_point = ax.c2p(x, y)   # top of the line (x-value, y_max)
            line = DashedLine(start_point, end_point, color=RED_C, dash_length=0.07, stroke_width=2)
            spectrum_3.add(line)

        x_spectrum_4 = [881.5, 882.5, 883.5]
        y_spectrum_4 = [2000, 1000, 500]
        
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

        # shift rightwards to avoid margin issues
        everything_group = VGroup(ax, x_label, curve_to_fill, curve_fill, y_label, target_line, spectrum_1, spectrum_2, spectrum_3, spectrum_4)  # Include all other objects that are part of the scene
        everything_group.shift(RIGHT*0.5)
        labels = VGroup(x_label)
        labels.shift(UP*0.3)

        #dot_1 = Dot(ax.i2gp(graph.t_min, graph))
        #dot_2 = Dot(ax.i2gp(graph.t_max, graph))
        self.add(ax, x_label, y_label, spectrum_1, spectrum_2, spectrum_3, spectrum_4)

        self.wait(1)

        self.play(self.camera.frame.animate.scale(0.25).move_to(spectrum_2))
        
        #self.play(self.camera.frame.animate.move_to(dot_2))

        #self.wait(1)

        # Add green line to scene
        self.play(Create(target_line), run_time = 0.5)

               # Apply curve
        self.play(Create(curve_to_fill), run_time = 3)

        # Apply the curve fill
        self.play(Write(curve_fill, lag_ratio=50, run_time = 1))

        # Apply wave
        self.play(ApplyWave(curve_to_fill, amplitude = 0.1), 
                  ApplyWave(curve_fill, amplitude = 0.05))
        #self.play(Restore(self.camera.frame))
        #self.wait()



class testMetapeak(MovingCameraScene):
    
    def construct(self):

        # Define the axes
        ax = Axes(
            x_range=[850, 925, 50],
       
            y_range=[0, 15000, 5000],
            axis_config={"color": BLUE, "include_numbers": True}
        )

        # Add labels for x and y axis
        x_label = MathTex("m/z").next_to(ax.x_axis, DOWN)
        y_label = MathTex("Counts").next_to(ax.y_axis, UP)


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


        x_spectrum_2 = [875.5, 876.5, 877.5]
        y_spectrum_2 = [2000, 1500, 1000]
    
        # Make spectrum 1
        spectrum_2 = VGroup()  # Create a vector group for the lines
        for (x, y) in zip(x_spectrum_2, y_spectrum_2):
            start_point = ax.c2p(x, 0)  # bottom of the line (x-value, y_min)
            end_point = ax.c2p(x, y)   # top of the line (x-value, y_max)
            line = DashedLine(start_point, end_point, color=RED_C, dash_length=0.07, stroke_width=2)
            spectrum_1.add(line)
        
        self.add(ax, x_label, y_label, spectrum_1, spectrum_2)

        self.play(self.camera.frame.animate.scale(0.4).move_to(spectrum_1))


