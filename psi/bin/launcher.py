"""Pipe Stress Infinity Command Line Interface (CLI).

$ psi               # Run PSI in console mode
$ psi file.inp      # Run a script in console
$ psi -i file.inp   # Run file and enter interactive mode
"""

import argparse
import sys
import time
import os
import logging
from contextlib import redirect_stdout

from tqdm import tqdm

from psi import VERSION
from psi.app import App


class TqdmLoggingHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)

            # force writing to default __stdout__
            # which should be the terminal screen
            with redirect_stdout(sys.__stdout__):
                tqdm.write(msg)

            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


def setup_logger(outfile, errfile):
    # logger setup
    tqdmlogger = logging.getLogger("tqdm")
    tqdmlogger.setLevel(logging.INFO)
    tqdmlogger.addHandler(TqdmLoggingHandler())

    stdoutlogger = logging.getLogger("stdout")
    stdoutlogger.setLevel(logging.INFO)
    stdoutlogger.addHandler(logging.FileHandler(outfile))

    stderrlogger = logging.getLogger("stderr")
    stderrlogger.setLevel(logging.INFO)
    stderrlogger.addHandler(logging.FileHandler(errfile))


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

    cwd = os.path.abspath(os.curdir)
    if args.file:
        basename, ext = os.path.splitext(os.path.basename(args.file))
        if ext == ".inp":
            outfile = os.path.join(cwd, basename+".out")
            errfile = os.path.join(cwd, basename+".err")
        else:
            raise IOError("Invalid file type or extension")

        # setup loggers
        setup_logger(outfile, errfile)

        app = App()

        # print("starting the gui application")
        num_lines = sum(1 for line in open(args.file, "r"))
        with open(args.file, "r") as fp:
            null = open(os.devnull, "w")
            bar = tqdm(fp, total=num_lines)

            header = ("PSI Design and Analysis\n"
                      "Version: %s \n"
                      "Design Codes: All Codes\n\n"
                      "Input File: %s \n" %
                      (VERSION, args.file))

            bar.write(header)

            bar.set_description("Processing...")
            for lno, line in enumerate(bar, 1):
                # suppress all superfluous output by interp
                with redirect_stdout(null):

                    try:
                        app.interp.push(line)

                    # catch all errors
                    except:
                        app.interp.showtraceback()
                        break

                if lno == num_lines:
                    bar.set_description("Done!")
                else:
                    time.sleep(0.01)

        if args.is_interactive:
            app.run()

    else:
        # print("starting the shell")
        outfile = os.path.join(cwd, "file.out")
        errfile = os.path.join(cwd, "file.err")

        setup_logger(outfile, errfile)

        App().run()


if __name__ == "__main__":
    App().run()
