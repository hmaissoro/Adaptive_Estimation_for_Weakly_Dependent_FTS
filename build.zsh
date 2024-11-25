#!/bin/zsh

# Define the package directory and output directory
PACKAGE_DIR="adaptiveFTS"
OUTPUT_DIR="out"

# Get the architecture and operating system
ARCH=$(uname -m)
OS=$(uname -s)

# Navigate to the package directory
cd $PACKAGE_DIR

# Check the package
echo "Checking the package..."
R CMD check .

# Build the package and specify the output directory
echo "Building the package..."
R CMD build .

# Compile the package
echo "Compiling the package..."
R CMD INSTALL .

# Build the documentation
# make sure the last line of the DESCRIPTION file is empty
if [ "$(tail -n 1 DESCRIPTION)" != "" ]; then
    echo "" >> DESCRIPTION
fi
echo "Building the documentation..."
Rscript -e "devtools::document()"
Rscript -e "devtools::build_manual()"
# output should be in ../adaptiveFTS_[VER].pdf

PACKAGE_NAME=$(ls | grep .tar.gz)
cd ..
DOCS_NAME=$(ls | grep .pdf)

mkdir -p $OUTPUT_DIR
mv "$PACKAGE_DIR/$PACKAGE_NAME" "$OUTPUT_DIR/${PACKAGE_NAME%.tar.gz}_${OS}_${ARCH}.tar.gz"

mkdir -p docs
mv $DOCS_NAME docs

echo "Done."