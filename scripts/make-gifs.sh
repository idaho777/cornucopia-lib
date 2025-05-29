#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 [YYMMDD]"
  echo "  YYMMDD   date subfolder under runs/fits- (defaults to today)"
  exit 1
}
# accept zero or one argument
if [[ $# -gt 1 ]]; then
  usage
fi

# 1) Figure out where the script lives (project root)
# 1) Figure out where the script lives (project root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# 2) Date: either passed in or today
if [[ $# -eq 1 ]]; then
  YYMMDD="$1"
else
  YYMMDD=$(date +%y%m%d)
fi

BASE="${PROJECT_DIR}/runs/fits-${YYMMDD}"
IMGDIR="${BASE}/img"

if [[ ! -d "$IMGDIR" ]]; then
  echo "❌ Cannot find image directory: $IMGDIR" >&2
  exit 1
fi

# 4) Find unique prefixes: strip off everything from "_T" onward
#    We'll collect them in a variable, one per line.
prefixes=$(
  find "$IMGDIR" -maxdepth 1 -type f -name 'FittingTest_s\[*\]_T\[*\]*' \
    -exec basename {} \; |
    sed 's/_a.*//' |
    sort -u
)

if [[ -z "$prefixes" ]]; then
  echo "⚠️  No matching files in $IMGDIR"
  exit 0
fi

echo "Found prefixes:"
while IFS= read -r p; do
  echo "  $p"
done <<EOF
$prefixes
EOF

# 5) Build GIF for each prefix
while IFS= read -r prefix; do
  echo "⏳  Building GIF for $prefix …"

  # escape literal [ and ] so the shell glob works
  safe="${prefix//[/\\[}"
  safe="${safe//]/\\]}"

  # collect all frames
  frames=("$IMGDIR"/${safe}*)

  echo "    found ${#frames[@]} frames"
  if ((${#frames[@]})); then
    out="${BASE}/${prefix}.gif"
    magick -delay 1 -loop 0 "${frames[@]}" "$out"
    echo "✅  Saved $out"
  else
    echo "⚠️   No frames found for $prefix"
  fi
done <<EOF
$prefixes
EOF

echo "🎉  All GIFs built!"
