#!/usr/bin/env bash
#set -e

# --------------------------------------------------------------------------------
# By default we do *nothing* until you ask.  At the end, if you asked for no
# steps at all, we'll turn them all on.
# --------------------------------------------------------------------------------
DO_CMAKE=false
DO_MAKE=false
DO_TEST=false
DO_PLOT=false
DO_GIFS=false

CMAKE_ARGS=""
MAKE_ARGS="-j$(sysctl -n hw.ncpu)"
EXEC_ARGS=""
PLOT_ARGS=""
GIF_ARGS=""

STEPS_SPECIFIED=false

# --------------------------------------------------------------------------------
# Usage / help
# --------------------------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $0 [cmake] [make] [run] [plot] [gifs] [--cmake ARGS] [--make ARGS] [--run ARGS] [--plot CSV] [--gifs DATE]

If you give no step-names (cmake/make/run/plot/gifs), all five will run by default.
If you list one or more of cmake/make/run/plot/gifs, *only* those will run (in the
order shown below), plus any of the --* argument-flags you supply.

Commands:
  cmake         Configure project with cmake
  make          Build with make
  run           Run the test executable
  plot          Run plot script
  gifs          Run GIF-maker script

Options:
  --cmake ARGS  Passed verbatim to cmake
  --make ARGS   Passed verbatim to make
  --run ARGS    Passed verbatim to the test binary
  --plot CSV    Passed as argument to scripts/plot.py
  --gifs DATE   Passed as argument to scripts/make-gifs.sh

EOF
}

# --------------------------------------------------------------------------------
# Parse positional commands + options
# --------------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
  cmake)
    DO_CMAKE=true
    STEPS_SPECIFIED=true
    shift
    ;;
  make)
    DO_MAKE=true
    STEPS_SPECIFIED=true
    shift
    ;;
  run)
    DO_TEST=true
    STEPS_SPECIFIED=true
    shift
    ;;
  plot)
    DO_PLOT=true
    STEPS_SPECIFIED=true
    shift
    ;;
  gifs)
    DO_GIFS=true
    STEPS_SPECIFIED=true
    shift
    ;;
  --cmake)
    CMAKE_ARGS="$2"
    shift 2
    ;;
  --make)
    MAKE_ARGS="$2"
    shift 2
    ;;
  --run)
    EXEC_ARGS="$2"
    shift 2
    ;;
  --plot)
    PLOT_ARGS="$2"
    shift 2
    ;;
  --gifs)
    GIF_ARGS="$2"
    shift 2
    ;;
  help | --help | -h)
    usage
    exit 0
    ;;
  *)
    echo "❌ Unknown argument: $1"
    usage
    exit 1
    ;;
  esac
done

# --------------------------------------------------------------------------------
# If the user never specified ANY of the five steps, turn them ALL on
# --------------------------------------------------------------------------------
if ! $STEPS_SPECIFIED; then
  DO_CMAKE=true
  DO_MAKE=true
  DO_TEST=true
  DO_PLOT=true
  DO_GIFS=true
fi

# --------------------------------------------------------------------------------
# Run steps in canonical order
# --------------------------------------------------------------------------------

if $DO_CMAKE; then
  echo "🛠️  Running cmake..."
  cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug -S . $CMAKE_ARGS
fi

if $DO_MAKE; then
  echo "⚙️  Running make..."
  make -C build $MAKE_ARGS
fi

if $DO_TEST; then
  echo "🧲  Running tests..."
  ./build/Test/Test $EXEC_ARGS
fi

if $DO_PLOT; then
  echo "📊  Plotting..."
  if [[ -n "$PLOT_ARGS" ]]; then
    python scripts/plot.py "$PLOT_ARGS"
  else
    python scripts/plot.py
  fi
fi

if $DO_GIFS; then
  echo "🎥  Generating GIFs... disabled"
  # if [[ -n "$GIF_ARGS" ]]; then
  #   ./scripts/make-gifs.sh "$GIF_ARGS"
  # else
  #   ./scripts/make-gifs.sh
  # fi
fi
