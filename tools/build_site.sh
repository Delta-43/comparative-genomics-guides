#!/usr/bin/env bash
# Build the website into ./_site from a staging copy of ONLY the site sources.
# Why: Quarto scans the whole project directory; a clean copy keeps large local data folders
# (FASTQ, BAM, pipeline results) out of the scan and out of the published site.
#   tools/build_site.sh            # build once  -> _site/
#   tools/build_site.sh preview    # live preview on http://localhost:4200
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
stage="$(mktemp -d "${TMPDIR:-/tmp}/guides-site.XXXXXX")"
trap 'rm -rf "$stage"' EXIT
QUARTO="${QUARTO:-quarto}"

# Keep this list in sync with project.render in _quarto.yml
mkdir -p "$stage/pipelines"
cp -r "$root/_quarto.yml" "$root/index.qmd" "$root/about.qmd" "$root/assets" "$root/guides" "$stage/"
cp "$root/pipelines/README.md" "$stage/pipelines/"

if [[ "${1:-}" == "preview" ]]; then
  (cd "$stage" && "$QUARTO" preview --port 4200 --no-browser)
else
  (cd "$stage" && "$QUARTO" render)
  rm -rf "$root/_site" && cp -r "$stage/_site" "$root/_site"
  echo "Site built: $root/_site/index.html"
fi
