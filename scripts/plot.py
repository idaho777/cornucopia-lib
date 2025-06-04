import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import shutil

from pathlib import Path
from tqdm.contrib.concurrent import process_map
from datetime import datetime
from matplotlib.axes import Axes

from generate_index import generate_index
from make_runs_manifest import make_runs_manifest


def read_csv_data(
    filename: str,
) -> tuple[
    list[tuple[float, float]],
    tuple[float, float, float, float, float, float, float, float],
    list[tuple[float, float]],
    tuple[float, float, float, float, float, float, float, float],
    list[tuple[float, float]],
]:
    """
    Reads a file with three blocks:
      Sampled points
      Fitted Anticlothoid
      Actual Anticlothoid

    each block has a header line, then data lines.
    For the fitted/actual blocks, the first data line is the 4-value transform,
    then a header "x,y", then (x,y) point lines.
    """
    # load & strip
    lines = [L.strip() for L in open(filename, "r") if L.strip()]

    # find where each block starts
    i_s = lines.index("Sampled points")
    i_f = lines.index("Fitted Anticlothoid")
    i_a = lines.index("Actual Anticlothoid")

    # --- sampled points: skip header + its "x,y" line ---
    sampled_lines = lines[i_s + 2 : i_f]
    sampled_pts = [tuple(map(float, L.split(","))) for L in sampled_lines]

    # --- fitted: the very next line is transform, then skip its "x,y" header ---
    fitted_data = tuple(map(float, lines[i_f + 2].split(",")))
    fitted_pts_lines = lines[i_f + 4 : i_a]
    fitted_pts = [tuple(map(float, L.split(","))) for L in fitted_pts_lines]

    # --- actual: same pattern as fitted ---
    actual_data = tuple(map(float, lines[i_a + 2].split(",")))
    actual_pts_lines = lines[i_a + 4 :]
    actual_pts = [tuple(map(float, L.split(","))) for L in actual_pts_lines]

    return sampled_pts, fitted_data, fitted_pts, actual_data, actual_pts


def draw_circle(
    ax: Axes,
    data: tuple[float, ...],
    curve_color: str,
    arrow_scale: float = 0.5,
) -> None:
    tx, ty, a, theta = data[:4]
    t = np.linspace(0.0, 2.0 * np.pi, 200)
    _ = ax.plot(tx + a * np.cos(t), ty + a * np.sin(t), "--", color=curve_color)
    _ = ax.plot(tx, ty, "o", color=curve_color)

    cos_t, sin_t = math.cos(theta), math.sin(theta)
    ux = np.array([cos_t, sin_t]) * a * arrow_scale
    uy = np.array([-sin_t, cos_t]) * a * arrow_scale

    _ = ax.arrow(tx, ty, ux[0], ux[1], color="red", width=0.02, head_width=0.1)
    _ = ax.arrow(tx, ty, uy[0], uy[1], color="blue", width=0.02, head_width=0.1)


def format_circle_info(prefix: str, data: tuple, pre=True) -> str:
    tx, ty, a, theta, start_ang, end_ang, length, start_rad, end_rad, error = data
    if pre:
        return f"{prefix}:\ncenter=({tx:.2f},{ty:.2f})\nr={a:.2f}\nθ={theta:.2f},\n\nstart_angle={start_ang:.2f}°\nend_angle={end_ang:.2f}°\nlength={length:.2f}\nstart_radius={start_rad:.2f}\nend_radius={end_rad:.2f}\nerror={error:.2f}"
    return f"{prefix}:\n({tx:.2f},{ty:.2f})=center\n{a:.2f}=r\n{theta:.2f}=θ,\n\n{start_ang:.2f}°=start_angle\n{end_ang:.2f}°=end_angle\n{length:.2f}=length\n{start_rad:.2f}=start_radius\n{end_rad:.2f}=end_radius\n{error:.2f}=error"


def plot_scene(csv_path: Path, output_dir: Path) -> None:
    sampled, fitted_t, fitted, actual_t, actual = read_csv_data(csv_path)

    fitted_color = "limegreen"
    actual_color = "darkorange"

    # 2×2 grid
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    # make room on the right for a single legend
    fig.subplots_adjust(right=0.75)

    # --- subplot 1: Sampled vs Fitted ---
    ax = axes[0, 0]
    if sampled:
        xs, ys = zip(*sampled)
        ax.plot(xs, ys, "k.-", label="Sampled")
    if fitted:
        xs, ys = zip(*fitted)
        ax.plot(xs, ys, "-", color=fitted_color, label="Fitted")
    if fitted_t:
        draw_circle(ax, fitted_t, fitted_color)
    ax.set_title("Sampled vs Fitted")
    ax.set_aspect("equal")

    # --- subplot 2: Sampled vs Actual ---
    ax = axes[0, 1]
    if sampled:
        xs, ys = zip(*sampled)
        ax.plot(xs, ys, "k.-", label="Sampled")
    if actual:
        xs, ys = zip(*actual)
        ax.plot(xs, ys, "-", color=actual_color, label="Actual")
    if actual_t:
        draw_circle(ax, actual_t, actual_color)
    ax.set_title("Sampled vs Actual")
    ax.set_aspect("equal")

    # --- subplot 3: Fitted vs Actual ---
    ax = axes[1, 0]
    if fitted:
        xs, ys = zip(*fitted)
        ax.plot(xs, ys, "-", color=fitted_color, label="Fitted")
    if actual:
        xs, ys = zip(*actual)
        ax.plot(xs, ys, "-", color=actual_color, label="Actual")
    if fitted_t:
        draw_circle(ax, fitted_t, fitted_color)
    if actual_t:
        draw_circle(ax, actual_t, actual_color)
    ax.set_title("Fitted vs Actual")
    ax.set_aspect("equal")

    # --- subplot 4: All Curves ---
    ax = axes[1, 1]
    if sampled:
        xs, ys = zip(*sampled)
        ax.plot(xs, ys, "k.-", label="Sampled")
    if fitted:
        xs, ys = zip(*fitted)
        ax.plot(xs, ys, "-", color=fitted_color, label="Fitted")
    if actual:
        xs, ys = zip(*actual)
        ax.plot(xs, ys, "-", color=actual_color, label="Actual")
    if fitted_t:
        draw_circle(ax, fitted_t, fitted_color)
    if actual_t:
        draw_circle(ax, actual_t, actual_color)
    ax.set_title("All Curves")
    ax.set_aspect("equal")

    # grid on all
    for a in axes.flat:
        a.grid(True)

    # --- single legend on the right ---
    handles, labels = axes[1, 1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(0.9, 0.4),
        frameon=False,
        fontsize=9,
    )

    # --- info boxes to the right ---
    if fitted_t:
        fig.text(
            0.9,
            0.6,
            format_circle_info("Fitted", fitted_t, pre=True),
            color=fitted_color,
            va="center",
            ha="right",
            fontsize=9,
            # bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=fitted_color),
        )
    if actual_t:
        fig.text(
            0.9,
            0.6,
            format_circle_info("Actual", actual_t, pre=False),
            color=actual_color,
            va="center",
            ha="left",
            fontsize=9,
            # bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=actual_color),
        )

    stem = csv_path.stem
    outpath = output_dir / f"{stem}.png"
    fig.savefig(outpath, dpi=70, bbox_inches="tight")
    plt.close(fig)


def _worker(args: tuple[Path, Path]) -> None:
    csv_path, img_dir = args
    plot_scene(csv_path, img_dir)


if __name__ == "__main__":
    runs_dir = Path("./runs")
    yy_mm_dd = datetime.now().strftime("%y%m%d")
    f_generate_imgs = True
    f_generate_index = True

    if len(sys.argv) > 1:
        if sys.argv[1] == "--skip":
            print("Skipping plotting...")
            f_generate_imgs = False
            f_generate_index = False
        else:
            yy_mm_dd = sys.argv[1]

    fit_dir = Path(f"./runs/fits-{yy_mm_dd}")
    csv_dir = fit_dir / "csv"
    img_dir = fit_dir / "img"

    if f_generate_imgs:
        fit_dir.mkdir(parents=True, exist_ok=True)
        if img_dir.exists():
            shutil.rmtree(img_dir)
        img_dir.mkdir()

        csv_files = list(csv_dir.glob("*.csv"))

        # parallel map with progress bar
        process_map(
            _worker,
            [(p, img_dir) for p in csv_files],
            desc="Plotting CSVs",
            chunksize=10,
            max_workers=6,
        )

    if f_generate_index:
        generate_index(str(img_dir), str(fit_dir / "index.json"))

    make_runs_manifest(runs_dir, runs_dir / "runs.json")
