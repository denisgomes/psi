"""OpenPIPE Command Line Interface (OPCLI)

$ openpipe                  # Run OpenPIPE in console
$ openpipe file.py          # Run a script in console
$ openpipe -i file.py       # Run file and start the GUI
$ openpipe -i               # Start the GUI
"""

import argparse

from openpipe import VERSION
from openpipe.app import App


def main():
    parser = argparse.ArgumentParser(
                            prog="OpenPIPE",
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
                        help="run program in interactive mode"
                        )

    parser.add_argument("file",
                        metavar="FILE",
                        nargs="?",
                        help="run program from script file"
                        )

    # parse arguments
    args = parser.parse_args()
    if args.file:
        # print "running script file in batch mode"
        if args.is_interactive:
            # print "starting the gui application"
            App().run()
    elif args.is_interactive:
        print "starting the gui application"
    else:
        # print "starting the shell"
        App().run()


if __name__ == "__main__":
    App().run()
