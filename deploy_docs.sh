#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# Build the documentation from the MASTER_BRANCH or DEV_BRANCH
# and push it to TARGET_BRANCH.
MASTER_BRANCH="master"
DEV_BRANCH="development"
TARGET_BRANCH="gh-pages"

# Pull requests and commits to other branches shouldn't try to deploy
if [ "$TRAVIS_PULL_REQUEST" != "false" ] || { [ "$TRAVIS_BRANCH" != "$MASTER_BRANCH" ] && [ "$TRAVIS_BRANCH" != "$DEV_BRANCH" ]; }; then
    echo "Skipping deploy on $TRAVIS_BRANCH. We only deploy docs automatically from development and master."
    exit 0
fi

# Save some useful information
REPO=`git config remote.origin.url`
SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
SHA=`git rev-parse --verify HEAD`

# Add rsa keys to the ssh agent to push to GitHub
gpg --output ../id_maestro_rsa --batch --passphrase $DECRYPT_GITHUB_AUTH --decrypt id_maestro_rsa.enc
chmod 600 ../id_maestro_rsa
eval `ssh-agent -s`
ssh-add ../id_maestro_rsa
#ls ../id_rsa_travis

# Clone the existing gh-pages for this repo into out/
# Create a new empty branch if gh-pages doesn't exist yet (should only happen on first deply)
git clone $REPO out
cd out

# # Regenerate the API documentation
# git checkout $SOURCE_BRANCH || git checkout --orphan $SOURCE_BRANCH
# sphinx-apidoc -f -M -o docs/source/ .

# # Add new docs to the repo
# git add docs/source --all

# # Commit and push to SOURCE_BRANCH
# git commit -m "Regenerate API Documentation: ${SHA}" || true
# git push $SSH_REPO $SOURCE_BRANCH || true

git checkout $TARGET_BRANCH || git checkout --orphan $TARGET_BRANCH
cd ..

# Clean out existing contents. If it's the master branch, then 
# move dev docs to a temporary directory then move back again after 
# everything is deleted
if [ "$TRAVIS_BRANCH" = "$MASTER_BRANCH" ]; then
    mkdir tmp
    mv -r out/docs/dev tmp || true
    rm -rf out/docs/**/* || exit 0
    mv tmp/dev out/docs || true
    rmdir tmp
else 
    rm -rf out/docs/dev/**/* || exit 0
fi

# Pull from SOURCE_BRANCH again
git pull || true

# if on the dev branch, use the dev_layout.html template to get the 
# links correct
if [ "$TRAVIS_BRANCH" = "$DEV_BRANCH" ]; then
    mv sphinx_docs/source/_templates/dev_layout.html sphinx_docs/source/_templates/layout.html 
fi

# Build the Sphinx documentation
cd sphinx_docs
make html
cd ../

mkdir -p out/docs/
if [ "$TRAVIS_BRANCH" = "$MASTER_BRANCH" ]; then
    mv sphinx_docs/build/html/* out/docs
else 
    mkdir -p out/docs/dev/
    mv sphinx_docs/build/html/* out/docs/dev
fi
touch out/.nojekyll

# Now let's go have some fun with the cloned repo
cd out
git config user.name "Travis CI"
git config user.email "$COMMIT_AUTHOR_EMAIL"

echo "doing git add/commit/push"

# Commit the "changes", i.e. the new version.
# The delta will show diffs between new and old versions.
git add --all

# Exit if there are no docs changes
if git diff --staged --quiet; then
   echo "exiting with no docs changes"
   exit 0
fi

# Otherwise, commit and push
git commit -m "Deploy to GitHub Pages: ${SHA}"
git push $SSH_REPO $TARGET_BRANCH
cd ..

# Kill the ssh-agent
ssh-agent -k
