"""PSI Command Line Interface (CLI).

$ psi               # Run PSI in console mode
$ psi file.py       # Run a script in console
$ psi -i file.py    # Run file and enter interactive mode
"""

import argparse
import inspect
import sys
import os

from tqdm import tqdm

from psi import VERSION
from psi.app import App


def main():
    parser = argparse.ArgumentParser(
                            prog="PSI",
                            formatter_class=argparse.RawTextHelpFormatter,
                            description=__doc__,
                            # epilog=LICENSE,
                            usage=argparse.SUPPRESS
                            )

    parser.add_argument("-V",
                        "--version",
                        action="version",
                        version="%(prog)s " + VERSION,
                        help="show version number and exit"
                        )

    parser.add_argument("-i",
                        "--interactive",
                        dest="is_interactive",
                        action="store_true",
                        help="enter interactive mode after running script"
                        )

    parser.add_argument("file",
                        metavar="FILE",
                        nargs="?",
                        help="run program from script file"
                        )

    # parse arguments
    args = parser.parse_args()
    if args.file:
        # print("running script file in batch mode")
        app = App()

        tqdm.write("PSI Design and Analysis")
        tqdm.write("Version: %s" % VERSION)
        tqdm.write("Design Codes: All Codes")
        tqdm.write("")
        tqdm.write("Input File: %s" % args.file)
        tqdm.write("")

        # print("starting the gui application")
        num_lines = sum(1 for line in open(args.file, "r"))
        with open(args.file, "r") as fp:
            for line in tqdm(fp, total=num_lines):
                app.interp.push(line)

        if args.is_interactive:
            app.run()

        app.quit()

    else:
        # print("starting the shell")
        App().run()


if __name__ == "__main__":
    App().run()
