#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# Build the documentation from the MAIN_BRANCH or DEV_BRANCH
# and push it to TARGET_BRANCH.
MAIN_BRANCH="main"
DEV_BRANCH="development"
TARGET_BRANCH="gh-pages"

mkdir out

# if on the dev branch, use the dev_layout.html template to get the 
# links correct
if [ "$GITHUB_BRANCH" = "$DEV_BRANCH" ]; then
    mv sphinx_docs/source/_templates/dev_layout.html sphinx_docs/source/_templates/layout.html 
fi

# Build the Sphinx documentation
cd sphinx_docs
make html
cd ../

mkdir -p out/docs/
if [ "$GITHUB_BRANCH" = "$MAIN_BRANCH" ]; then
    mv sphinx_docs/build/html/* out/docs
else 
    mkdir -p out/docs/dev/
    mv sphinx_docs/build/html/* out/docs/dev
fi
