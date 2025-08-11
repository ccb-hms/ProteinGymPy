#!/usr/bin/env python3
"""
Basic test from the original test directory.
"""

import sys
import unittest


class TestBasic(unittest.TestCase):
    """Basic test class."""
    
    def test_hello_world(self):
        """Test basic functionality."""
        # This test is a placeholder from the original test.py
        result = "Hello World"
        self.assertEqual(result, "Hello World")
    
    def test_python_version(self):
        """Test that we're running on Python 3.13+."""
        self.assertGreaterEqual(sys.version_info.major, 3)
        self.assertGreaterEqual(sys.version_info.minor, 13)


def main():
    """Main function for direct execution."""
    print("Hello World")


if __name__ == "__main__":
    unittest.main()