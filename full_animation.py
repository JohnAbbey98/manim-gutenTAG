"""
Manim animation illustrating the gutenTAG metapeak algorithm for MALDI-MSI data.

Faithfully animates each step of the generateMetapeaks() pipeline:
  1. countPeaks — raw peak-count histogram
  2. smoothPeakCounts — Gaussian convolution
  3. generateSeedMz — local maxima detection
  4. segmentPeakCounts — watershed segmentation
  5. estimateMetapeaks — center, max, limits
  6. assignMetapeaks — panel target matching

Render: manim -pql full_animation.py MetapeakAlgorithm
"""

from manim import *
import numpy as np
from scipy.signal import find_peaks
from collections import deque

# Note: Tex() is used for all text (not Text()) due to manimpango 0.4.4 Cython type bugs.


# ---------------------------------------------------------------------------
# Colour constants (matching plotQC.R conventions, adapted for dark background)
# ---------------------------------------------------------------------------
RAW_COUNT_COLOR = WHITE
SMOOTH_COLOR = ORANGE
SEED_COLOR = TEAL
MAX_COLOR = GREEN
TARGET_COLOR = PURPLE
METAPEAK_FILL = GREY
TARGETED_FILL = RED


# ---------------------------------------------------------------------------
# Synthetic data generation
# ---------------------------------------------------------------------------

def generate_synthetic_maldi_data(
    mz_range=(850, 920),
    bin_size=0.05,
    seed=42,
):
    """
    Generate a realistic-looking MALDI peak-count histogram with isotopic envelopes.

    Species are arranged in dense neighbouring clusters so isotopic envelopes
    interleave, illustrating the challenge the metapeak algorithm solves.

    Returns
    -------
    mz_bins : np.ndarray
        m/z bin centers
    raw_counts : np.ndarray
        Simulated peak counts per bin
    panel_targets : np.ndarray
        Panel target m/z values (the "expected" feature masses)
    species_mz : list[float]
        Monoisotopic m/z of each placed species (for reference)
    """
    rng = np.random.default_rng(seed)
    mz_bins = np.arange(mz_range[0], mz_range[1], bin_size)
    raw_counts = np.zeros_like(mz_bins, dtype=float)

    # 5 clearly separated species, each with monoisotopic peak + 2-3 isotopic peaks
    species_mz_list = [858.0, 872.0, 886.0, 900.0, 912.0]

    # Isotopic envelopes: monoisotopic is strongest, then decreasing
    isotope_patterns = [
        [1.0, 0.65, 0.30, 0.10],   # species 1: M, M+1, M+2, M+3
        [1.0, 0.70, 0.35, 0.10],         # species 2: M, M+1, M+2, M+3
        [1.0, 0.75, 0.40, 0.12],   # species 3
        [1.0, 0.60, 0.25],         # species 4
        [1.0, 0.70, 0.35, 0.10],   # species 5
    ]

    # Moderate variation in abundance
    base_counts = [350, 280, 400, 220, 320]

    for i, (mono_mz, base_count) in enumerate(zip(species_mz_list, base_counts)):
        pattern = isotope_patterns[i % len(isotope_patterns)]
        for j, rel_intensity in enumerate(pattern):
            peak_mz = mono_mz + j * 1.003  # carbon-13 spacing
            count = base_count * rel_intensity * rng.uniform(0.85, 1.15)
            # Add as a narrow Gaussian (sigma ~ 0.08 Da) to the histogram
            sigma_peak = 0.08
            peak_contribution = count * np.exp(
                -0.5 * ((mz_bins - peak_mz) / sigma_peak) ** 2
            )
            raw_counts += peak_contribution

    # Constant low-level background noise
    raw_counts += rng.poisson(2, size=len(mz_bins)).astype(float)
    raw_counts = np.maximum(raw_counts, 0)

    # Panel targets: offset 0.2-0.5 Da from monoisotopic peak to show realistic mass shift
    panel_targets = np.array([
        858.0 + 0.3,   # +0.3 Da shift
        872.0 - 0.4,   # -0.4 Da shift
        886.0 + 0.2,   # +0.2 Da shift
        900.0 - 0.3,   # -0.3 Da shift
        912.0 + 0.5,   # +0.5 Da shift
    ])

    return mz_bins, raw_counts, panel_targets, species_mz_list


# ---------------------------------------------------------------------------
# Algorithm helpers (faithful to gutenTAG R implementations)
# ---------------------------------------------------------------------------

def mz_scaling_factor(mz_bins):
    """Matches .mzScalingFactor() in gutenTAG."""
    return (mz_bins[-1] - mz_bins[0]) / len(mz_bins)


def smooth_peak_counts(mz_bins, counts, smooth_factor=1.0):
    """
    Gaussian convolution matching smoothPeakCounts.R.
    Uses circular convolution via np.convolve (open mode for simplicity).
    """
    n = len(counts)
    sigma = smooth_factor / mz_scaling_factor(mz_bins)
    x = np.linspace(-n / 2, n / 2, n)
    kernel = np.exp(-x ** 2 / sigma ** 2)
    kernel /= kernel.sum()
    smoothed = np.convolve(counts, kernel, mode="same")
    return smoothed


def smooth_peak_counts_variable(mz_bins, counts, sigma_value):
    """Smooth with an explicit sigma value (for animation interpolation)."""
    if sigma_value <= 0:
        return counts.copy()
    n = len(counts)
    x = np.linspace(-n / 2, n / 2, n)
    kernel = np.exp(-x ** 2 / sigma_value ** 2)
    kernel /= kernel.sum()
    return np.convolve(counts, kernel, mode="same")


def generate_seed_mz(mz_bins, smoothed, detection_threshold=0, sparsity=3):
    """
    Find local maxima of smoothed histogram. Matches generateSeedMz.R.
    Uses scipy.signal.find_peaks instead of splus2R::peaks.
    """
    counts = smoothed.copy()
    counts[counts < detection_threshold] = 0
    min_distance = max(1, int(sparsity / mz_scaling_factor(mz_bins)))
    peaks_idx, _ = find_peaks(counts, distance=min_distance)
    # Filter out peaks at zero
    peaks_idx = peaks_idx[counts[peaks_idx] > detection_threshold]
    return mz_bins[peaks_idx], peaks_idx


def segment_peak_counts(mz_bins, smoothed, seed_indices, detection_threshold=0):
    """
    1D watershed segmentation matching segmentPeakCounts.R (EBImage::propagate).
    Each bin is assigned to its nearest seed; bins below threshold are label 0.
    """
    n = len(mz_bins)
    labels = np.zeros(n, dtype=int)
    counts = smoothed.copy()
    counts[counts < detection_threshold] = 0

    # BFS from each seed simultaneously (Voronoi-like propagation)
    queue = deque()
    visited = np.zeros(n, dtype=bool)
    for i, idx in enumerate(seed_indices):
        labels[idx] = i + 1
        visited[idx] = True
        queue.append(idx)

    while queue:
        pos = queue.popleft()
        for neighbor in [pos - 1, pos + 1]:
            if 0 <= neighbor < n and not visited[neighbor] and counts[neighbor] > 0:
                visited[neighbor] = True
                labels[neighbor] = labels[pos]
                queue.append(neighbor)

    return labels


def segment_peak_counts_stepwise(mz_bins, smoothed, seed_indices, detection_threshold=0):
    """
    BFS watershed segmentation that records each step for animation.
    Returns (final_labels, steps) where steps is a list of label snapshots.
    """
    n = len(mz_bins)
    labels = np.zeros(n, dtype=int)
    counts = smoothed.copy()
    counts[counts < detection_threshold] = 0

    queue = deque()
    visited = np.zeros(n, dtype=bool)
    for i, idx in enumerate(seed_indices):
        labels[idx] = i + 1
        visited[idx] = True
        queue.append(idx)

    steps = [labels.copy()]

    while queue:
        # Process one "wavefront" — all items currently in the queue
        frontier_size = len(queue)
        for _ in range(frontier_size):
            pos = queue.popleft()
            for neighbor in [pos - 1, pos + 1]:
                if 0 <= neighbor < n and not visited[neighbor] and counts[neighbor] > 0:
                    visited[neighbor] = True
                    labels[neighbor] = labels[pos]
                    queue.append(neighbor)
        steps.append(labels.copy())

    return labels, steps


def estimate_metapeaks(mz_bins, raw_counts, labels):
    """
    Compute metapeak parameters matching estimateMetapeaks.R.
    Returns list of dicts with keys: center, max_mz, lower_limit, upper_limit.
    """
    num_segments = labels.max()
    metapeaks = []
    for k in range(1, num_segments + 1):
        mask = labels == k
        if not mask.any():
            continue
        peak_mzs = mz_bins[mask]
        peak_counts = raw_counts[mask]

        center = np.average(peak_mzs, weights=peak_counts) if peak_counts.sum() > 0 else peak_mzs.mean()
        max_mz = peak_mzs[np.argmax(peak_counts)]

        # 1st and 99th percentile limits, re-centered around max_mz
        cumsum_freq = np.cumsum(peak_counts) / peak_counts.sum()
        lower_idx = np.searchsorted(cumsum_freq, 0.01)
        upper_idx = np.searchsorted(cumsum_freq, 0.99)
        lower_limit_raw = peak_mzs[min(lower_idx, len(peak_mzs) - 1)]
        upper_limit_raw = peak_mzs[min(upper_idx, len(peak_mzs) - 1)]
        half_width = (upper_limit_raw - lower_limit_raw) / 2
        lower_limit = max_mz - half_width
        upper_limit = max_mz + half_width

        metapeaks.append({
            "center": center,
            "max_mz": max_mz,
            "lower_limit": lower_limit,
            "upper_limit": upper_limit,
            "label": k,
        })
    return metapeaks


def assign_metapeaks_to_targets(metapeaks, panel_targets, mz_threshold=1.5):
    """
    Match panel targets to nearest metapeak (k-NN with threshold).
    Returns dict: target_mz -> metapeak dict (or None).
    """
    assignments = {}
    for target in panel_targets:
        best = None
        best_dist = mz_threshold
        for mp in metapeaks:
            dist = abs(mp["max_mz"] - target)
            if dist < best_dist:
                best_dist = dist
                best = mp
        assignments[target] = best
    return assignments


# ---------------------------------------------------------------------------
# Pre-compute all data for the animation
# ---------------------------------------------------------------------------
MZ_BINS, RAW_COUNTS, PANEL_TARGETS, SPECIES_MZ = generate_synthetic_maldi_data()
SMOOTH_FACTOR = 1.0
SIGMA_TARGET = SMOOTH_FACTOR / mz_scaling_factor(MZ_BINS)
SMOOTHED = smooth_peak_counts(MZ_BINS, RAW_COUNTS, smooth_factor=SMOOTH_FACTOR)
DETECTION_THRESHOLD = np.max(RAW_COUNTS) * 0.025
SEED_MZ, SEED_IDX = generate_seed_mz(MZ_BINS, SMOOTHED, detection_threshold=DETECTION_THRESHOLD, sparsity=3)
LABELS, BFS_STEPS = segment_peak_counts_stepwise(MZ_BINS, SMOOTHED, SEED_IDX, detection_threshold=DETECTION_THRESHOLD)
METAPEAKS = estimate_metapeaks(MZ_BINS, RAW_COUNTS, LABELS)
ASSIGNMENTS = assign_metapeaks_to_targets(METAPEAKS, PANEL_TARGETS)

# Pre-compute smoothed curves at various sigma values for progressive animation
N_SIGMA_STEPS = 90  # More steps for smooth animation (was 30)
SIGMA_STEPS = np.linspace(0, SIGMA_TARGET, N_SIGMA_STEPS)
SMOOTHED_CACHE = [smooth_peak_counts_variable(MZ_BINS, RAW_COUNTS, s) for s in SIGMA_STEPS]

# Segment colour palette (muted, for ~15 segments)
SEGMENT_COLORS = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
    "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac",
    "#86bcb6", "#8cd17d", "#b6992d", "#d4a6c8", "#499894",
    "#f1ce63", "#d37295", "#a0cbe8", "#ffbe7d", "#8b8b8b",
]


# ---------------------------------------------------------------------------
# Main animation scene
# ---------------------------------------------------------------------------

class MetapeakAlgorithm(MovingCameraScene):
    def construct(self):
        self.camera.frame.save_state()

        # Normalise y-values for manim coordinate system
        self.y_scale = 1.0 / RAW_COUNTS.max() * 12000
        self.raw_scaled = RAW_COUNTS * self.y_scale
        self.smooth_scaled = SMOOTHED * self.y_scale * 1.2  # 20% boost for smoothed peaks
        self.smoothed_cache_scaled = [s * self.y_scale * 1.2 for s in SMOOTHED_CACHE]

        # Step indicator (persistent, follows camera)
        #self.step_label = Tex("", font_size=28, color=YELLOW).to_corner(UL, buff=0.5)
        self.step_label = Tex("", font_size=28, color=YELLOW).to_edge(UP, buff=0.5)
        self.add(self.step_label)

        self.title_card()
        self.step1_count_peaks()
        self.step2_smooth()
        self.step3_seeds()
        self.step4_segment()
        self.step5_estimate()
        self.step6_assign()
        self.step7_summary()

    # -- Helpers --

    def update_step_label(self, text):
        #new_label = Tex(text, font_size=28, color=YELLOW).to_corner(UL, buff=0.15)
        new_label = Tex(text, font_size=28, color=YELLOW).to_edge(UP, buff=0.15)
        self.play(Transform(self.step_label, new_label), run_time=0.5)

    def make_line_graph(self, mz, counts, color, stroke_width=1.5, opacity=1.0):
        """Create a manim line graph from arrays using the scene axes."""
        graph = self.ax.plot_line_graph(
            mz, counts,
            add_vertex_dots=False,
            line_color=color,
            stroke_width=stroke_width,
        )
        if opacity < 1.0:
            graph.set_opacity(opacity)
        return graph

    def mz_to_scene(self, mz_val):
        """Convert an m/z value to a scene x-coordinate using the axes."""
        return self.ax.c2p(mz_val, 0)[0]

    # -- Steps --

    def title_card(self):
        title = Tex("gutenTAG: Metapeak Algorithm", font_size=44, color=WHITE)
        subtitle = Tex(
            "MALDI Mass Spectrometry Peak Segmentation",
            font_size=24, color=GREY_B,
        ).next_to(title, DOWN, buff=0.3)
        group = VGroup(title, subtitle)
        self.play(FadeIn(group), run_time=1.5)
        self.wait(1.5)
        self.play(FadeOut(group), run_time=1)

    def step1_count_peaks(self):
        self.update_step_label("Step 1: countPeaks")

        # Create axes
        y_max = float(self.raw_scaled.max() * 1.15)
        self.ax = Axes(
            x_range=[850, 920, 10],
            y_range=[0, y_max, y_max / 4],
            x_length=12,
            y_length=6,
            axis_config={"color": BLUE, "include_numbers": True},
            x_axis_config={"numbers_to_include": np.arange(860, 920, 10)},
            y_axis_config={"numbers_to_include": []},
        )
        x_label = Tex("m/z", font_size=22).next_to(self.ax.x_axis, DOWN, buff=0.3)
        y_label = Tex("Peak Count", font_size=22).next_to(self.ax.y_axis, UP, buff=0.4)

        self.play(FadeIn(self.ax), FadeIn(x_label), FadeIn(y_label), run_time=1)

        # Store for later use
        self.x_label = x_label
        self.y_label = y_label

        # Raw histogram line
        self.raw_graph = self.make_line_graph(
            MZ_BINS, self.raw_scaled, RAW_COUNT_COLOR, stroke_width=1.0
        )
        self.play(Create(self.raw_graph), run_time=4)
        self.wait(0.5)

        # Label isotope peaks for one species to educate viewers (no zoom)
        isotope_labels = VGroup()
        ref_mz = 872.0  # monoisotopic mass of species 2
        for j, label_text in enumerate(["M", "M+1", "M+2", "M+3"]):
            peak_mz = ref_mz + j * 1.003
            # Find nearest bin to get y position
            bin_idx = np.argmin(np.abs(MZ_BINS - peak_mz))
            peak_y = self.raw_scaled[bin_idx]
            label = Tex(label_text, font_size=15, color=YELLOW).move_to(
                self.ax.c2p(peak_mz, peak_y + self.raw_scaled.max() * 0.04)
            )
            isotope_labels.add(label)

        self.play(
            LaggedStart(*[FadeIn(l, shift=UP * 0.1) for l in isotope_labels], lag_ratio=0.15),
            run_time=1,
        )
        self.wait(1)
        self.play(FadeOut(isotope_labels), run_time=0.5)
        self.wait(0.5)

    def step2_smooth(self):
        self.update_step_label("Step 2: smoothPeakCounts")

        # Show the Gaussian kernel shape as a brief inset
        kernel_ax = Axes(
            x_range=[-3, 3, 1],
            y_range=[0, 1.1, 0.5],
            x_length=2.5,
            y_length=1.5,
            axis_config={"color": GREY_B, "stroke_width": 1, "include_numbers": False},
        ).to_corner(UR, buff=0.5)
        kernel_x = np.linspace(-3, 3, 100)
        kernel_y = np.exp(-0.5 * kernel_x ** 2)
        kernel_curve = kernel_ax.plot_line_graph(
            kernel_x, kernel_y,
            add_vertex_dots=False,
            line_color=SMOOTH_COLOR,
            stroke_width=2,
        )
        kernel_label = Tex(
            r"Gaussian kernel ($\sigma$ = " + f"{SMOOTH_FACTOR:.1f})",
            font_size=16, color=SMOOTH_COLOR,
        ).next_to(kernel_ax, DOWN, buff=0.1)
        kernel_group = VGroup(kernel_ax, kernel_curve, kernel_label)
        self.play(FadeIn(kernel_group), run_time=1)
        self.wait(1)

        # Static sigma label
        sigma_text = Tex(
            r"Gaussian convolution: $\sigma$ = " + f"{SMOOTH_FACTOR:.1f}",
            font_size=20, color=SMOOTH_COLOR,
        ).next_to(self.ax, DOWN, buff=0.6)
        self.play(FadeIn(sigma_text), run_time=0.5)

        # Progressive sigma animation using pre-computed cache
        sigma_tracker = ValueTracker(0)

        def get_smooth_graph():
            t = sigma_tracker.get_value()
            idx = int(t * (N_SIGMA_STEPS - 1))
            idx = min(idx, N_SIGMA_STEPS - 1)
            y_data = self.smoothed_cache_scaled[idx]
            return self.ax.plot_line_graph(
                MZ_BINS, y_data,
                add_vertex_dots=False,
                line_color=SMOOTH_COLOR,
                stroke_width=2.5,
            )

        smooth_graph = always_redraw(get_smooth_graph)
        self.add(smooth_graph)

        # Animate sigma from 0 to 1 (linear for smooth cache interpolation)
        self.play(sigma_tracker.animate.set_value(1.0), run_time=5, rate_func=linear)
        self.wait(0.5)

        # Fade out kernel inset and sigma text
        self.play(FadeOut(kernel_group), FadeOut(sigma_text), run_time=0.5)

        # Replace always_redraw with static final graph
        self.remove(smooth_graph)
        self.smooth_graph = self.make_line_graph(
            MZ_BINS, self.smooth_scaled, SMOOTH_COLOR, stroke_width=2.5
        )
        self.add(self.smooth_graph)

        # Zoom into dense cluster C (~870-875) to show merged isotope peaks
        zoom_center = self.ax.c2p(872, self.smooth_scaled.max() * 0.5)
        self.play(
            self.camera.frame.animate.scale(0.35).move_to(zoom_center),
            run_time=1.5,
        )
        merge_label = Tex(
            "Isotopic Peaks Merged", font_size=14, color=SMOOTH_COLOR,
        ).move_to(self.ax.c2p(872, self.smooth_scaled.max() * 0.85))
        self.play(FadeIn(merge_label), run_time=0.5)
        self.wait(2)
        self.play(FadeOut(merge_label), run_time=0.3)
        self.play(Restore(self.camera.frame), run_time=1.5)
        self.wait(0.5)

    def step3_seeds(self):
        self.update_step_label("Step 3: generateSeedMz")

        self.seed_dots = VGroup()
        self.seed_lines = VGroup()

        for mz_val, idx in zip(SEED_MZ, SEED_IDX):
            y_val = self.smooth_scaled[idx]
            # Dot on the smoothed curve
            dot_pos = self.ax.c2p(mz_val, y_val)
            dot = Dot(dot_pos, radius=0.06, color=SEED_COLOR)
            self.seed_dots.add(dot)

            # Dashed vertical line to x-axis
            start = self.ax.c2p(mz_val, 0)
            end = self.ax.c2p(mz_val, y_val)
            line = DashedLine(start, end, color=SEED_COLOR, dash_length=0.08, stroke_width=3.5)
            self.seed_lines.add(line)

        self.play(
            LaggedStart(*[FadeIn(d, scale=0.5) for d in self.seed_dots], lag_ratio=0.08),
            run_time=2,
        )
        self.play(
            LaggedStart(*[Create(l) for l in self.seed_lines], lag_ratio=0.08),
            run_time=2,
        )
        self.wait(1)

    def _build_segment_polygons(self, labels_snapshot):
        """Build filled polygons for segments present in a label snapshot."""
        fills = VGroup()
        num_segments = labels_snapshot.max()
        for k in range(1, num_segments + 1):
            mask = labels_snapshot == k
            if not mask.any():
                continue
            indices = np.where(mask)[0]
            seg_mz = MZ_BINS[indices]
            seg_smooth = self.smooth_scaled[indices]

            points = []
            points.append(self.ax.c2p(seg_mz[0], 0))
            for m, s in zip(seg_mz[::2], seg_smooth[::2]):
                points.append(self.ax.c2p(m, s))
            points.append(self.ax.c2p(seg_mz[-1], seg_smooth[-1]))
            points.append(self.ax.c2p(seg_mz[-1], 0))

            color = SEGMENT_COLORS[k % len(SEGMENT_COLORS)]
            polygon = Polygon(
                *points,
                fill_color=color,
                fill_opacity=0.25,
                stroke_width=0,
            )
            fills.add(polygon)
        return fills

    def step4_segment(self):
        self.update_step_label("Step 4: segmentPeakCounts")

        # Show detection threshold line
        threshold_y = DETECTION_THRESHOLD * self.y_scale
        threshold_line = DashedLine(
            self.ax.c2p(850, threshold_y),
            self.ax.c2p(920, threshold_y),
            color=GREY, dash_length=0.15, stroke_width=1,
        )
        threshold_text = Tex("detection threshold", font_size=16, color=GREY).next_to(
            threshold_line, RIGHT, buff=0.2
        )
        self.play(Create(threshold_line), FadeIn(threshold_text), run_time=1)
        self.wait(0.5)

        # Animated watershed flooding using pre-computed BFS steps
        flood_tracker = ValueTracker(0)
        n_steps = len(BFS_STEPS)

        def get_flood_fills():
            t = flood_tracker.get_value()
            step_idx = int(t * (n_steps - 1))
            step_idx = min(step_idx, n_steps - 1)
            return self._build_segment_polygons(BFS_STEPS[step_idx])

        flood_fills = always_redraw(get_flood_fills)
        self.add(flood_fills)

        # Slow pedagogical flooding (~10 seconds)
        self.play(
            flood_tracker.animate.set_value(1.0),
            run_time=10,
            rate_func=linear,
        )
        self.wait(0.5)

        # Replace always_redraw with static final fills
        self.remove(flood_fills)
        self.segment_fills = self._build_segment_polygons(LABELS)
        self.add(self.segment_fills)

        # Zoom to a boundary between two segments
        if len(SEED_MZ) >= 2:
            boundary_mz = (SEED_MZ[4] + SEED_MZ[5]) / 2 if len(SEED_MZ) > 5 else (SEED_MZ[0] + SEED_MZ[1]) / 2
            zoom_center = self.ax.c2p(boundary_mz, self.smooth_scaled.max() * 0.3)
            self.play(
                self.camera.frame.animate.scale(0.3).move_to(zoom_center),
                run_time=1.5,
            )
            boundary_label = Tex(
                "watershed boundary", font_size=14, color=WHITE,
            ).move_to(self.ax.c2p(boundary_mz, self.smooth_scaled.max() * 0.15))
            self.play(FadeIn(boundary_label), run_time=0.5)
            self.wait(2)
            self.play(FadeOut(boundary_label), run_time=0.3)
            self.play(Restore(self.camera.frame), run_time=1.5)

        # Clean up threshold line
        self.play(FadeOut(threshold_line), FadeOut(threshold_text), run_time=0.5)
        self.wait(0.5)

    def step5_estimate(self):
        self.update_step_label("Step 5: estimateMetapeaks")

        self.max_lines = VGroup()
        self.center_lines = VGroup()
        self.metapeak_rects = VGroup()
        y_top = float(self.raw_scaled.max() * 1.1)

        for mp in METAPEAKS:
            # Max line (green dashed)
            max_line = DashedLine(
                self.ax.c2p(mp["max_mz"], 0),
                self.ax.c2p(mp["max_mz"], y_top * 0.6),
                color=MAX_COLOR, dash_length=0.08, stroke_width=3.5,
            )
            self.max_lines.add(max_line)

            # Center line (green solid, thinner)
            center_line = Line(
                self.ax.c2p(mp["center"], 0),
                self.ax.c2p(mp["center"], y_top * 0.5),
                color=MAX_COLOR, stroke_width=2.5, stroke_opacity=0.8,
            )
            self.center_lines.add(center_line)

            # Metapeak range rectangle (grey)
            rect = Rectangle(
                width=abs(self.mz_to_scene(mp["upper_limit"]) - self.mz_to_scene(mp["lower_limit"])),
                height=self.ax.c2p(0, y_top)[1] - self.ax.c2p(0, 0)[1],
                fill_color=METAPEAK_FILL,
                fill_opacity=0.15,
                stroke_width=0.5,
                stroke_color=GREY,
            )
            rect_center_x = (self.mz_to_scene(mp["lower_limit"]) + self.mz_to_scene(mp["upper_limit"])) / 2
            rect_center_y = (self.ax.c2p(0, 0)[1] + self.ax.c2p(0, y_top)[1]) / 2
            rect.move_to([rect_center_x, rect_center_y, 0])
            self.metapeak_rects.add(rect)

        # Zoom into one segment to show detail
        if METAPEAKS:
            focus_mp = METAPEAKS[len(METAPEAKS) // 2]
            zoom_center = self.ax.c2p(focus_mp["max_mz"], self.smooth_scaled.max() * 0.4)
            self.play(
                self.camera.frame.animate.scale(0.3).move_to(zoom_center),
                run_time=1.5,
            )

            # Show max, center, limits for this one metapeak
            focus_idx = len(METAPEAKS) // 2
            self.play(Create(self.max_lines[focus_idx]), run_time=0.8)

            max_label = Tex("max", font_size=14, color=MAX_COLOR).next_to(
                self.max_lines[focus_idx].get_top(), RIGHT, buff=0.05
            )
            self.play(FadeIn(max_label), run_time=0.3)

            self.play(Create(self.center_lines[focus_idx]), run_time=0.8)
            center_label = Tex("center", font_size=14, color=MAX_COLOR).next_to(
                self.center_lines[focus_idx].get_top(), LEFT, buff=0.05
            )
            self.play(FadeIn(center_label), run_time=0.3)

            self.play(FadeIn(self.metapeak_rects[focus_idx]), run_time=0.8)
            self.wait(1.5)

            # Clean up labels and restore
            self.play(FadeOut(max_label), FadeOut(center_label), run_time=0.3)
            self.play(Restore(self.camera.frame), run_time=1.5)

        # Show all metapeak annotations at once
        remaining_max = VGroup(*[l for i, l in enumerate(self.max_lines) if i != len(METAPEAKS) // 2])
        remaining_center = VGroup(*[l for i, l in enumerate(self.center_lines) if i != len(METAPEAKS) // 2])
        remaining_rects = VGroup(*[r for i, r in enumerate(self.metapeak_rects) if i != len(METAPEAKS) // 2])

        self.play(
            LaggedStart(*[Create(l) for l in remaining_max], lag_ratio=0.03),
            LaggedStart(*[Create(l) for l in remaining_center], lag_ratio=0.03),
            LaggedStart(*[FadeIn(r) for r in remaining_rects], lag_ratio=0.03),
            run_time=2,
        )
        self.wait(1)

    def step6_assign(self):
        self.update_step_label("Step 6: assignMetapeaks")

        y_top = float(self.raw_scaled.max() * 1.1)

        # Purple target lines
        self.target_lines = VGroup()
        for target_mz in PANEL_TARGETS:
            line = Line(
                self.ax.c2p(target_mz, 0),
                self.ax.c2p(target_mz, y_top * 0.7),
                color=TARGET_COLOR, stroke_width=3.5,
            )
            self.target_lines.add(line)

        self.play(
            LaggedStart(*[FadeIn(l) for l in self.target_lines], lag_ratio=0.1),
            run_time=1.5,
        )
        self.wait(0.5)

        # Matching arrows and colour transitions
        arrows = VGroup()
        matched_rect_indices = []

        for target_mz, mp in ASSIGNMENTS.items():
            if mp is None:
                continue
            # Find the rect index for this metapeak
            for i, m in enumerate(METAPEAKS):
                if m["label"] == mp["label"]:
                    matched_rect_indices.append(i)
                    break

            # Curved arrow from target to metapeak max
            start = self.ax.c2p(target_mz, y_top * 0.65)
            end = self.ax.c2p(mp["max_mz"], y_top * 0.55)
            if abs(target_mz - mp["max_mz"]) > 0.01:
                arrow = CurvedArrow(start, end, color=YELLOW_C, stroke_width=2.5, tip_length=0.15)
            else:
                arrow = Arrow(start, end, color=YELLOW_C, stroke_width=2.5, tip_length=0.15, buff=0.05)
            arrows.add(arrow)

        self.play(
            LaggedStart(*[Create(a) for a in arrows], lag_ratio=0.15),
            run_time=2,
        )

        # Transition matched rectangles grey -> red
        animations = []
        for idx in matched_rect_indices:
            if idx < len(self.metapeak_rects):
                animations.append(
                    self.metapeak_rects[idx].animate.set_fill(TARGETED_FILL, opacity=0.2)
                )
        if animations:
            self.play(*animations, run_time=1.5)

        self.wait(0.5)

        # Zoom into one matched pair
        matched_targets = [t for t, m in ASSIGNMENTS.items() if m is not None]
        if matched_targets:
            focus_target = matched_targets[len(matched_targets) // 2]
            focus_mp = ASSIGNMENTS[focus_target]
            zoom_center = self.ax.c2p(
                (focus_target + focus_mp["max_mz"]) / 2,
                self.smooth_scaled.max() * 0.4,
            )
            self.play(
                self.camera.frame.animate.scale(0.3).move_to(zoom_center),
                run_time=1.5,
            )
            self.wait(2)
            self.play(Restore(self.camera.frame), run_time=1.5)

        # Clean up arrows
        self.play(FadeOut(arrows), run_time=0.5)
        self.wait(0.5)

    def step7_summary(self):
        self.update_step_label("")

        self.wait(0.5)

        # Build legend in top-right corner
        legend_items = [
            (Line(ORIGIN, RIGHT * 0.4, color=RAW_COUNT_COLOR, stroke_width=1.5), "Raw counts"),
            (Line(ORIGIN, RIGHT * 0.4, color=SMOOTH_COLOR, stroke_width=2.5), "Smoothed"),
            (DashedLine(ORIGIN, RIGHT * 0.4, color=SEED_COLOR, stroke_width=3.5, dash_length=0.06), "Seeds"),
            (Rectangle(width=0.4, height=0.2, fill_color=SEGMENT_COLORS[0], fill_opacity=0.25, stroke_width=0.5, stroke_color=SEGMENT_COLORS[0]), "Watershed segment"),
            (DashedLine(ORIGIN, RIGHT * 0.4, color=MAX_COLOR, stroke_width=3.5, dash_length=0.06), "Metapeak max"),
            (Line(ORIGIN, RIGHT * 0.4, color=TARGET_COLOR, stroke_width=3.5), "Panel target"),
            (Rectangle(width=0.4, height=0.2, fill_color=METAPEAK_FILL, fill_opacity=0.2, stroke_width=0.5), "Metapeak range"),
            (Rectangle(width=0.4, height=0.2, fill_color=TARGETED_FILL, fill_opacity=0.2, stroke_width=0.5), "Targeted range"),
        ]

        legend = VGroup()
        for i, (icon, label_text) in enumerate(legend_items):
            label = Tex(label_text, font_size=16, color=WHITE)
            row = VGroup(icon, label).arrange(RIGHT, buff=0.15)
            legend.add(row)

        legend.arrange(DOWN, buff=0.12, aligned_edge=LEFT)
        legend_bg = SurroundingRectangle(
            legend, color=WHITE, fill_color=BLACK, fill_opacity=0.7,
            buff=0.15, corner_radius=0.1, stroke_width=0.5,
        )
        legend_group = VGroup(legend_bg, legend).to_corner(UR, buff=0.3)

        self.play(FadeIn(legend_group), run_time=1.5)
        self.wait(3)
        self.play(FadeOut(legend_group), run_time=1)
        self.wait(1)
