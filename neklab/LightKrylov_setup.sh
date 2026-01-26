#!/bin/bash

# LightKrylov Setup Script – custom version
echo "LightKrylov Setup Script"
echo "-----------------------"
echo "This script will:"
echo "1. Install fpm with the desired package manager (if desired)."
echo "2. Clone the LightKrylov repository (if desired)."
echo "3. Build and optionally test LightKrylov."
echo "-----------------------"

# Install dependencies if the user agrees
read -p "Do you need to install fpm ? (y/n): " confirm
    if [[ "$confirm" =~ ^[yY]$ ]]; then
    # Ask the user for a package manager to install fpm
    echo "Select the package manager to install fpm:"
    echo "  1) conda"
    echo "  2) brew"
    echo "  3) pip"
    echo "  4) spack"
    read -p "Enter the number of your choice: " pm_choice

    case "$pm_choice" in
      1) PM="conda" ;;
      2) PM="brew" ;;
      3) PM="pip" ;;
      4) PM="spack" ;;
      *) echo "Invalid choice. Exiting."; exit 1 ;;
    esac

    echo "You selected $PM to install fpm."

    echo "Installing dependencies with $PM..."
    case "$PM" in
      conda)
        conda config --add channels conda-forge
        conda install -y fpm
        ;;
      brew)
        brew tap fortran-lang/fortran
        brew update
        brew install fpm
        ;;
      pip)
        # fpm is available on PyPI.
        pip install --upgrade fpm
        ;;
      spack)
        # Assume spack is already installed and configured
        spack install fpm
        ;;
    esac
else
    echo "Skipping dependencies installation."
fi

# Check for existing LightKrylov directory
if [ -d "LightKrylov" ]; then
    read -p "Directory LightKrylov exists. Remove it? (y/n): " confirm
    if [[ "$confirm" =~ ^[yY]$ ]]; then
        echo "Removing LightKrylov directory..."
        rm -rf LightKrylov
        clone_needed="yes"
    else
        echo "Retaining existing LightKrylov directory."
        clone_needed="no"
    fi
else
    clone_needed="yes"
fi

# Clone LightKrylov if needed
if [ "$clone_needed" = "yes" ]; then
    read -p "Clone the LightKrylov repository? (y/n): " confirm
    if [[ "$confirm" =~ ^[yY]$ ]]; then
        echo "Cloning LightKrylov..."
        git clone https://github.com/nekStab/LightKrylov.git
        cd LightKrylov
        git checkout -b pinned-version dcd85b3e027e2cda016e1f48d8a28f146c9eb5b3
    else
        echo "Skipping cloning."
        echo "Ensure LightKrylov directory is present for subsequent steps."
    fi
else
    cd LightKrylov
fi

# Clean any previous builds
fpm clean --all

# Build and install LightKrylov
echo "Installing LightKrylov..."
fpm install --profile release

# Run unit tests if desired
read -p "Run LightKrylov unit tests? (y/n): " confirm
if [[ "$confirm" =~ ^[yY]$ ]]; then
    fpm test --profile release
fi

echo "LightKrylov setup complete."

