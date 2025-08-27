#!/usr/bin/env python3
"""
HDMI: HGT Detection from MAGs in Individual
Entry point for pip installation
"""

import sys
import os

# Get the directory containing this script
script_dir = os.path.dirname(os.path.abspath(__file__))
bin_dir = os.path.join(script_dir, 'bin')

# Add the bin directory to the path
if bin_dir not in sys.path:
    sys.path.insert(0, bin_dir)

# Import and run the main function from bin/HDMI
import HDMI as hdmi_module
hdmi_module.main()
