#!/usr/bin/env python3
"""
HDMI: HGT Detection from MAGs in Individual
Entry point for pip installation
"""

import sys
import os

# Add the bin directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'bin'))

# Import and run the main function from bin/HDMI
from HDMI import main

if __name__ == "__main__":
    main()
