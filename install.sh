#!/usr/bin/env bash
# Create the btc-ccs conda environment and check that it works.
# Usage: ./install.sh [env-name]   (default env name: btc-ccs)
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
env_name="${1:-btc-ccs}"

# prefer mamba (much faster solves), fall back to conda
if command -v mamba >/dev/null 2>&1; then
  conda_cmd=mamba
elif command -v conda >/dev/null 2>&1; then
  conda_cmd=conda
else
  echo "Neither conda nor mamba was found. Install Miniforge first: https://github.com/conda-forge/miniforge" >&2
  exit 1
fi

# bioconda has no Apple Silicon build of CNTools, so build an Intel env (runs under Rosetta 2)
if [ "$(uname -s)" = "Darwin" ] && [ "$(uname -m)" = "arm64" ]; then
  echo "Apple Silicon detected: creating an osx-64 environment (requires Rosetta 2)"
  if ! arch -x86_64 /usr/bin/true 2>/dev/null; then
    echo "Rosetta 2 is not installed. Install it first with:" >&2
    echo "  softwareupdate --install-rosetta --agree-to-license" >&2
    exit 1
  fi
  export CONDA_SUBDIR=osx-64
fi

writable() {  # true if the directory is writable, or doesn't exist yet but its parent is
  if [ -d "$1" ]; then [ -w "$1" ]; else [ -w "$(dirname "$1")" ]; fi
}

# a machine-wide install (e.g. /opt/miniforge3 owned by root) can't hold user envs;
# fall back to ~/.conda, which conda also searches when activating by name
conda_base="$(conda info --base)"
if writable "$conda_base/envs"; then
  env_prefix="$conda_base/envs/$env_name"
else
  env_prefix="$HOME/.conda/envs/$env_name"
  echo "$conda_base/envs is not writable; installing to $env_prefix instead"
  mkdir -p "$HOME/.conda/envs"
  if ! conda info --json | grep -q "\"$HOME/.conda/envs\""; then
    echo "Adding $HOME/.conda/envs to envs_dirs in ~/.condarc so the env can be activated by name"
    conda config --append envs_dirs "$HOME/.conda/envs"
  fi
fi
if ! writable "$conda_base/pkgs"; then
  export CONDA_PKGS_DIRS="$HOME/.conda/pkgs"
fi

if [ -e "$env_prefix" ]; then
  echo "$env_prefix already exists. Remove it first with:  conda env remove -p $env_prefix" >&2
  exit 1
fi

echo "Creating conda environment '$env_name' with $conda_cmd"
"$conda_cmd" env create -p "$env_prefix" -f "$repo_dir/environment.yml"

# keep later installs into this env on the same platform
if [ "${CONDA_SUBDIR:-}" = "osx-64" ]; then
  conda run -p "$env_prefix" conda config --env --set subdir osx-64
fi

echo
echo "Testing the installation"
conda run --no-capture-output -p "$env_prefix" bash "$repo_dir/test_install.sh"

echo
echo "Done. Activate the environment with:  conda activate $env_name"
