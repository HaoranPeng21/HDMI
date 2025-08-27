#!/bin/bash

# HDMI Installation Script
# HGT Detection from MAGs in Individual

echo "=== HDMI Installation Script ==="
echo "HGT Detection from MAGs in Individual"
echo ""

# Check if conda is available
if command -v conda &> /dev/null; then
    echo "✓ Conda found"
    
    # Create conda environment
    echo "Creating conda environment 'HDMI'..."
    conda env create -f environment.yml
    
    if [ $? -eq 0 ]; then
        echo "✓ Conda environment created successfully"
        echo ""
        echo "To activate the environment and use HDMI:"
        echo "  conda activate HDMI"
        echo "  HDMI --help"
        echo ""
        echo "Or install HDMI globally:"
        echo "  pip install -e ."
        echo "  HDMI --help"
    else
        echo "✗ Failed to create conda environment"
        exit 1
    fi
else
    echo "Conda not found. Installing with pip..."
    
    # Install with pip
    pip install -e .
    
    if [ $? -eq 0 ]; then
        echo "✓ HDMI installed successfully"
        echo ""
        echo "To use HDMI:"
        echo "  HDMI --help"
    else
        echo "✗ Failed to install HDMI"
        exit 1
    fi
fi

echo ""
echo "=== Installation Complete ==="
echo "For more information, see README.md"
