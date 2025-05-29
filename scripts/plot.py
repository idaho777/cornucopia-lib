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
    tuple[float, float, float, float] | None,
    list[tuple[float, float]],
    tuple[float, float, float, float] | None,
    list[tuple[float, float]],
]:
    sampled_pts: list[tuple[float, float]] = []
    fitted_transform: tuple[float, float, float, float] | None = None
    fitted_pts: list[tuple[float, float]] = []
    actual_transform: tuple[float, float, float, float] | None = None
    actual_pts: list[tuple[float, float]] = []

    section = None
    saw_header = False

    with open(filename, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line == "Sampled points":
                section = "sampled"
                saw_header = False
                continue
            if line == "Fitted Anticlothoid":
                section = "fitted"
                saw_header = False
                continue
            if line == "Actual Anticlothoid":
                section = "actual"
                saw_header = False
                continue

            if section == "sampled":
                if not saw_header:
                    saw_header = True
                    continue
                parts = line.split(",")
                if len(parts) >= 2:
                    try:
                        x, y = float(parts[0]), float(parts[1])
                        sampled_pts.append((x, y))
                    except ValueError:
                        pass

            elif section == "fitted":
                if fitted_transform is None:
                    if line.startswith("startx"):
                        continue
                    parts = line.split(",")
                    if len(parts) == 4:
                        try:
                            fitted_transform = tuple(map(float, parts))
                        except ValueError:
                            fitted_transform = None
                elif not saw_header:
                    if line.startswith("x"):
                        saw_header = True
                        continue
                else:
                    parts = line.split(",")
                    if len(parts) >= 2:
                        try:
                            x, y = float(parts[0]), float(parts[1])
                            fitted_pts.append((x, y))
                        except ValueError:
                            pass

            elif section == "actual":
                if actual_transform is None:
                    if line.startswith("startx"):
                        continue
                    parts = line.split(",")
                    if len(parts) == 4:
                        try:
                            actual_transform = tuple(map(float, parts))
                        except ValueError:
                            actual_transform = None
                elif not saw_header:
                    if line.startswith("x"):
                        saw_header = True
                        continue
                else:
                    parts = line.split(",")
                    if len(parts) >= 2:
                        try:
                            x, y = float(parts[0]), float(parts[1])
                            actual_pts.append((x, y))
                        except ValueError:
                            pass

    return sampled_pts, fitted_transform, fitted_pts, actual_transform, actual_pts


def draw_circle(
    ax: Axes,
    transform: tuple[float, float, float, float],
    curve_color: str,
    arrow_scale: float = 0.5,
) -> None:
    tx, ty, a, theta = transform
    t = np.linspace(0.0, 2.0 * np.pi, 200)
    ax.plot(tx + a * np.cos(t), ty + a * np.sin(t), "--", color=curve_color)
    ax.plot(tx, ty, "o", color=curve_color)

    cos_t, sin_t = math.cos(theta), math.sin(theta)
    ux = np.array([cos_t, sin_t]) * a * arrow_scale
    uy = np.array([-sin_t, cos_t]) * a * arrow_scale

    ax.arrow(tx, ty, ux[0], ux[1], color="red", width=0.02, head_width=0.1)
    ax.arrow(tx, ty, uy[0], uy[1], color="blue", width=0.02, head_width=0.1)


def format_circle_info(
    prefix: str, transform: tuple[float, float, float, float]
) -> str:
    tx, ty, a, theta = transform
    return f"{prefix}:\ncenter=({tx:.2f},{ty:.2f})\nr={a:.2f}\nθ={theta:.2f}"


def plot_scene(csv_path: Path, output_dir: Path) -> None:
    sampled, fitted_t, fitted, actual_t, actual = read_csv_data(csv_path)

    fitted_color = "darkorange"
    actual_color = "dimgrey"

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.subplots_adjust(right=0.8)

    # 1) Sampled vs Fitted
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
    ax.legend(loc="upper left")
    ax.set_aspect("equal")

    # 2) Sampled vs Actual
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
    ax.legend(loc="upper left")
    ax.set_aspect("equal")

    # 3) Fitted vs Actual
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
    ax.legend(loc="upper left")
    ax.set_aspect("equal")

    # 4) All three
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
    ax.legend(loc="upper left")
    ax.set_aspect("equal")

    for a in axes.flat:
        a.grid(True)

    # Two info boxes on the right
    if fitted_t:
        fig.text(
            0.82,
            0.65,
            format_circle_info("Fitted", fitted_t),
            color=fitted_color,
            va="center",
            ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=fitted_color),
        )
    if actual_t:
        fig.text(
            0.82,
            0.35,
            format_circle_info("Actual", actual_t),
            color=actual_color,
            va="center",
            ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=actual_color),
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
