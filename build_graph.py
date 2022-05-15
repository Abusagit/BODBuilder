import os
import sys
import argparse
import logging
import shutil

package_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(package_dir, 'graph_dir'))

import graph_structure as gs


class TqdmHandler(logging.Handler):
    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)  # , file=sys.stderr)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


def main():
    parser = argparse.ArgumentParser(description="De-Bruijn graph builder, symplyfier and drawer")
    parser.add_argument("-i", "--input", nargs="+", type=str, help="Path to your file with data for graph")
    parser.add_argument("-k", "--kmer-size", type=int, help="Size of a kmer to be used for graph building")
    parser.add_argument("-b", "--bad_cov", type=float, default=100, help="Threshold for clipping edges defined \
                                                                                    as Edge_cov / mean(Edge_cov)")
    parser.add_argument("--draw", action="store_true", help="Stay this option to get drawn final graph")
    parser.add_argument("--force", action="store_true", help="Force override dir]")
    parser.add_argument("-o", "--outdir")

    args = parser.parse_args()

    force_message = ""
    if os.path.isdir(args.outdir) and not args.force:
        raise FileExistsError(f"Directory {args.outdir} already exists! Specify another one or use '--force' to override")
    elif os.path.isdir(args.outdir) and args.force:
        shutil.rmtree(args.outdir, ignore_errors=True)

        force_message = f"Force overriding directory {args.outdir}"

    os.mkdir(args.outdir)

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(os.path.join(args.outdir, "graph_building.log"), mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    #
    # tqdm_handler = TqdmHandler()
    # tqdm_handler.setLevel(logging.DEBUG)

    root.addHandler(console_handler)
    root.addHandler(file_handler)
    # root.addHandler(tqdm_handler)
    if force_message:
        root.warning(force_message)

    input_files = '\n'.join(args.input)
    root.info(f"Started building graph on: \n{input_files}")
    root.debug(f"Parsed arguments are: {args}")
    root.info(f"Results will be saved in {args.outdir}/")

    graph = gs.DBGraph(args.input, args.kmer_size, outdir=args.outdir, ratio=args.bad_cov)

    graph.complete_graph_building()


if __name__ == '__main__':
    main()
