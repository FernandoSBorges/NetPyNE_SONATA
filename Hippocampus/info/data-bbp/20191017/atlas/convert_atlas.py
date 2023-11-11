#!/usr/bin/env python2

import os
import sys
import argparse

import numpy as np

from voxcell.nexus.voxelbrain import Atlas


def main(args):
    atlas = Atlas.open(args.url, cache_dir=args.atlas_cache)

    layers = args.layers.split(",")

    layer_ratio = np.array(args.layer_ratio.split(","), dtype=np.float32)
    layer_ratio /= np.sum(layer_ratio)
    y01 = np.insert(np.cumsum(layer_ratio), 0, 0)

    height = atlas.load_data('height')
    for layer, y0, y1 in zip(layers, y01[:-1], y01[1:]):
        print >>sys.stderr, "PH[%s]..." % layer
        raw = np.zeros(height.raw.shape + (2,), dtype=np.float32)
        raw[..., 0] = y0 * height.raw
        raw[..., 1] = y1 * height.raw
        height.with_data(raw).save_nrrd(os.path.join(args.output_dir, "[PH]%s.nrrd" % layer))

    print >>sys.stderr, "PH[y]..."
    distance = atlas.load_data('distance')
    distance.save_nrrd(os.path.join(args.output_dir, "[PH]y.nrrd"))

    print >>sys.stderr, "Done!"


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "url",
        help="Atlas folder / URL"
    )
    parser.add_argument(
        "-o", "--output-dir",
        help="Output folder"
    )
    parser.add_argument(
        "--atlas-cache",
        help="Atlas cache dir",
        default=None
    )
    parser.add_argument(
        "--layers", help="Layer names ('bottom' to 'top', comma-separated)", required=True
    )
    parser.add_argument(
        "--layer-ratio", help="Layer thickness ratio (comma-separated)", required=True
    )
    main(parser.parse_args())