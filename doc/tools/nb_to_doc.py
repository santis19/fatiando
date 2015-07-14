#! /usr/bin/env python
"""
Convert empty IPython notebook to a sphinx doc page.

"""
import sys
from subprocess import check_call as sh


def convert_nb(nbname):

    # Clear notebook output
    sh(["ipython", "nbconvert", "--to", "notebook", "--inplace",
        "--ClearOutputPreprocessor.enabled=True", nbname])

    # Execute the notebook
    sh(["ipython", "nbconvert", "--to", "notebook",
        "--execute", "--inplace", nbname])

    # Convert to .rst for Sphinx
    sh(["ipython", "nbconvert", "--to", "rst", nbname])

    # Clear notebook output
    sh(["ipython", "nbconvert", "--to", "notebook", "--inplace",
        "--ClearOutputPreprocessor.enabled=True", nbname])


if __name__ == "__main__":

    for nbname in sys.argv[1:]:
        convert_nb(nbname)
