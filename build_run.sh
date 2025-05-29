#!/bin/bash
# set -e

# By default, run cmake, make, and run
DO_CMAKE=true
DO_MAKE=true
DO_TEST=true
DO_PLOT=true
DO_GIFS=true

CMAKE_ARGS=""
MAKE_ARGS="-j$(sysctl -n hw.ncpu)"
EXEC_ARGS=""
PLOT_ARGS=""
GIF_ARGS=""

# Parse arguments
if [[ "$1" == "help" || "$1" == "--help" || "$1" == "-h" ]]; then
  echo "Usage: $0 [cmake|make|run|plot [CSV]|gifs] [--cmake ARGS] [--make ARGS] [--run ARGS] [--plot CSV] [--gifs DATE]"
  echo ""
  echo "Commands:"
  echo "  cmake         Run only cmake step"
  echo "  make          Run only make step"
  echo "  run           Run only test executable"
  echo "  plot [CSV]    Run plot script (optionally with CSV file)"
  echo "  gifs [DATE]   Run GIF creation script (optionally with DATE)"
  echo ""
  echo "Options:"
  echo "  --cmake ARGS  Pass ARGS to cmake"
  echo "  --make ARGS   Pass ARGS to make"
  echo "  --run ARGS    Pass ARGS to test executable"
  echo "  --plot CSV    Pass CSV file to plot script"
  echo "  --gifs DATE   Pass DATE to GIF script"
  echo ""
  echo "If no command is given, runs cmake, make, test, plot, and gifs in order."
  exit 0
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
  cmake)
    DO_CMAKE=true
    DO_MAKE=false
    DO_TEST=false
    DO_PLOT=false
    DO_GIFS=false
    shift
    ;;
  make)
    DO_CMAKE=false
    DO_MAKE=true
    DO_TEST=false
    DO_PLOT=false
    DO_GIFS=false
    shift
    ;;
  run)
    DO_CMAKE=false
    DO_MAKE=false
    DO_TEST=true
    DO_PLOT=false
    DO_GIFS=false
    shift
    ;;
  plot)
    DO_CMAKE=false
    DO_MAKE=false
    DO_TEST=false
    DO_PLOT=true
    DO_GIFS=false
    shift
    ;;
  gifs)
    DO_CMAKE=false
    DO_MAKE=false
    DO_TEST=false
    DO_PLOT=false
    DO_GIFS=true
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
  *)
    echo "❌ Unknown argument: $1"
    exit 1
    ;;
  esac
done

# Run selected steps
$DO_CMAKE && {
  echo "🛠️Running CMake..."
  cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug -S . "$CMAKE_ARGS"
}

$DO_MAKE && {
  echo "⚙️ Running Make..."
  make -C build "$MAKE_ARGS"
}

$DO_TEST && {
  echo "🧲 Running ./build/Test/Test $EXEC_ARGS"
  ./build/Test/Test "$EXEC_ARGS"
}

$DO_PLOT && {
  echo "📊 Plotting..."
  if [[ -n "$PLOT_ARGS" ]]; then
    echo "📄 plot.py: $PLOT_ARGS"
    python scripts/plot.py "$PLOT_ARGS"
  else
    echo "📂 No CSV file specified. Plot new data."
    python scripts/plot.py
  fi
}

$DO_GIFS && {
  echo "🎥 Making gifs..."
  if [[ -n "$GIF_ARGS" ]]; then
    echo " Using DATE: $GIF_ARGS"
    ./scripts/make-gifs.sh "$GIF_ARGS"
  else
    ./scripts/make-gifs.sh
  fi
}
