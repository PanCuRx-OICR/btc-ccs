#!/usr/bin/env bash
# Check that the eCCS and gCCS classifiers run and reproduce the bundled example outputs.
# Usage: ./test_install.sh   (inside the btc-ccs conda env, or any R with the dependencies)

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/btc-ccs-test.XXXXXX")" || { echo "Could not create a temporary directory" >&2; exit 1; }
trap 'rm -rf "$tmp_dir"' EXIT

failed=0
pass() { printf '  \033[32mPASS\033[0m  %s\n' "$1"; }
fail() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; failed=1; }

echo "Checking R and packages"
if ! command -v Rscript >/dev/null 2>&1; then
  fail "Rscript not found on PATH (did you 'conda activate btc-ccs'?)"
  exit 1
fi
if ! r_version="$(Rscript --version 2>&1)"; then
  fail "Rscript found but could not run:"
  echo "$r_version" | sed 's/^/          /'
  if [[ "$r_version" == *"Bad CPU type"* ]]; then
    echo "          This R is an Intel build; install Rosetta 2 with:"
    echo "          softwareupdate --install-rosetta --agree-to-license"
  fi
  exit 1
fi
pass "$(echo "$r_version" | head -1)"

for pkg in data.table optparse stringr dplyr tibble CNTools; do
  if Rscript -e "suppressMessages(library($pkg))" >/dev/null 2>&1; then
    pass "R package $pkg"
  else
    fail "R package $pkg could not be loaded"
  fi
done

run_classifier() {  # name, expected output, script + args...
  local name=$1 expected=$2; shift 2
  local out="$tmp_dir/$name.txt" log="$tmp_dir/$name.log" want="$tmp_dir/$name.expected"
  tr -d '\r' <"$expected" >"$want"
  if ! Rscript "$@" -o "$out" >"$log" 2>&1; then
    fail "$name failed to run; last lines of its log:"
    tail -5 "$log" | sed 's/^/          /'
  elif diff -q "$want" "$out" >/dev/null; then
    pass "$name matches $(basename "$expected"): $(tail -1 "$out" | tr '\t' ' ')"
  else
    fail "$name output differs from $(basename "$expected"):"
    diff "$want" "$out" | sed 's/^/          /'
  fi
}

echo "Running classifiers on the example data"
run_classifier eCCS "$repo_dir/data/eCCS.txt" \
  "$repo_dir/call_btc_eCCS.R" -i "$repo_dir/data/tpm.txt"
run_classifier gCCS "$repo_dir/data/gCCS.txt" \
  "$repo_dir/call_btc_gCCS.R" -s "$repo_dir/data/example.seg" -v "$repo_dir/data/mean_vaf.txt"

echo
if [ "$failed" -eq 0 ]; then
  echo "All checks passed. btc-ccs is ready to use."
else
  echo "Some checks failed (see above)."
  exit 1
fi
