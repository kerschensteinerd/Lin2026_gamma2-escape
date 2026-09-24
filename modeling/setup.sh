#!/usr/bin/env bash
# setup.sh - check Python/NEURON/deps, install what pip can, point to docs otherwise.
# Usage: ./setup.sh
set -e

echo "Checking Python..."
python3 --version || { echo "Python 3.8+ required: https://www.python.org/downloads/"; exit 1; }

echo "Checking NEURON..."
if python3 -c "import neuron" 2>/dev/null; then
    python3 -c "from neuron import h; print('  NEURON OK, version', h.nrnversion())"
elif [[ "$(uname -s)" == MINGW* || "$(uname -s)" == MSYS* ]]; then
    # Native Windows (Git Bash) - no pip wheel; NEURON ships a binary installer only.
    echo "  NEURON not found. On native Windows, install via the official .exe:"
    echo "  -> https://github.com/neuronsimulator/nrn/releases"
    echo "  Step-by-step guide: https://neuron.yale.edu/ftp/neuron/nrn_mswin_install.pdf"
    echo "  (Python must already be installed separately.) Re-run this script after."
    echo "  Alternative: install NEURON inside WSL instead (pip install works there)."
    exit 1
else
    # macOS / Linux / WSL - pip install is the documented method.
    echo "  NEURON not found - installing via pip..."
    if ! pip3 install neuron; then
        echo "  pip install failed. Full install docs (binary/source alternatives):"
        echo "  -> https://nrn.readthedocs.io/en/latest/install/install_instructions.html"
        exit 1
    fi
fi

echo "Installing remaining Python dependencies..."
pip3 install -r requirements.txt

echo "Setup complete. Try: python3 scripts/example/run_example.py"
