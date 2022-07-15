#!/usr/bin/env python
"""
meta-sweeper - for performing parametric sweeps of simulated
metagenomic sequencing experiments.
Copyright (C) 2016 "Matthew Z DeMaere"

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from collections import defaultdict, OrderedDict


def filter_shared(d, shared):
    unshared = set(d.keys()) - shared
    for oi in unshared:
        del d[oi]


def inv_dict(d):
    cl = set()
    for v in d.values():
        cl.update(v)
    inv_d = {ci: set() for ci in cl}
    for oi, cl in d.items():
        for ci in cl:
            inv_d[ci].add(oi)
    return inv_d


def build_overlap_set_dict(obj2cl):
    cl2obj = inv_dict(obj2cl)
    ovl_dict = defaultdict(set)
    for neighbors in cl2obj.values():
        for oi in neighbors:
            ovl_dict[oi].update(neighbors)
    return ovl_dict


def build_clsize_dict(cl2obj, ovl_dict):
    out = {}
    for oi in cl2obj:
        out[oi] = {oj: len(cl2obj[oi] & cl2obj[oj]) for oj in ovl_dict[oi]}
    return out


def bcubed_recall(ovl_dict, c2c_size, t2t_size, td, weight=None):
    recall = 0.0
    if weight:
        for oi in td:
            ovl = ovl_dict[oi]
            ri = sum(weight[oj] * min(t2t_size[oi][oj], c2c_size[oi].get(oj, 0)) / t2t_size[oi][oj] for oj in ovl)
            wi = sum(weight[oj] for oj in ovl)
            recall += ri / wi
    else:
        for oi in td:
            ovl = ovl_dict[oi]
            ri = sum(min(t2t_size[oi][oj], c2c_size[oi].get(oj, 0)) / t2t_size[oi][oj] for oj in ovl)
            recall += ri / len(ovl)

    return recall / len(td)


def bcubed_precision(ovl_dict, c2c_size, t2t_size, cd, weight=None):
    precision = 0.0
    if weight:
        for oi in cd:
            ovl = ovl_dict[oi]
            pi = sum(weight[oj] * min(t2t_size[oi].get(oj, 0), c2c_size[oi][oj]) / c2c_size[oi][oj] for oj in ovl)
            wi = sum(weight[oj] for oj in ovl)
            precision += pi / wi
    else:
        for oi in cd:
            ovl = ovl_dict[oi]
            pi = sum(min(t2t_size[oi].get(oj, 0), c2c_size[oi][oj]) / c2c_size[oi][oj] for oj in ovl)
            precision += pi / len(ovl)
    return precision / len(cd)


def bcubed_fscore(return_queue, result, td, cd, weight=None):
    cd = cd.copy()
    td = td.copy()
    shared = set(td.keys()) & set(cd.keys())
    assert len(shared) > 0, 'There was no overlap between objects in truth table and clustering'

    # ratio of shared to all objects in analysis
    shared_ovlp = len(shared) / len(set(cd.keys()) | set(td.keys()))
    # print(f'Shared overlap: {shared_ovlp}')
    shared2tr = len(shared) / len(td)
    # print(f'Shared fraction relative to truth: {shared2tr}')
    shared2cl = len(shared) / len(cd)
    # print(f'Shared fraction relative to clustering: {shared2cl}')

    filter_shared(cd, shared)
    filter_shared(td, shared)

    # print('Building overlap set lookup')
    cd_ovl = build_overlap_set_dict(cd)
    td_ovl = build_overlap_set_dict(td)

    # print('Bulding size lookup')
    c2c_size = build_clsize_dict(cd, cd_ovl)
    t2t_size = build_clsize_dict(td, td_ovl)

    # print('Calculating precision')
    pre = bcubed_precision(cd_ovl, c2c_size, t2t_size, cd, weight)
    # print('Calculating recall')
    rec = bcubed_recall(td_ovl, c2c_size, t2t_size, td, weight)

    result['precision'] = pre
    result['recall'] = rec
    result['f_score'] = 2.0*pre*rec / (pre+rec)
    result['truth_count'] = len(td)
    result['pred_count'] = len(cd)
    result['shared_count'] = len(shared)
    result['shared_overlap'] = shared_ovlp
    result['shared2clustering'] = shared2cl
    result['shared2truth'] = shared2tr
    return_queue.put(result)


if __name__ == '__main__':

    import truthtable as tt
    import argparse
    import pandas as pd
    import sys
    import os
    import glob
    from multiprocessing import Pool, Manager

    def write_msg(stream, msg):
        stream.write(msg + '\n')

    parser = argparse.ArgumentParser(description='Calculate extended bcubed metric')
    parser.add_argument('truth', metavar='TRUTH', help='Truth table (json format)')
    parser.add_argument('pred', metavar='PREDICTION', help='Prediction MCL (file or directory)')
    parser.add_argument('output', metavar='OUTPUT', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout, help='Output file (stdout)')
    parser.add_argument('-N', '--ncpu', type=int, help='Number of CPU to use for multi-input processing')
    parser.add_argument('--tfmt', choices=['json', 'yaml'], default='json',
                        help='Data format of truth table [json]')
    parser.add_argument('-w', '--weighted', default=False, action='store_true', help='Use truth object weights')
    parser.add_argument('--hard', action='store_true', default=False, help='Extract hard truth prior to analysis')
    args = parser.parse_args()

    try:
        # read truth and convert to basic soft table
        truth = tt.read_truth(args.truth, args.tfmt)
        assert len(truth) > 0, f'Truth table contains no assignments: {args.truth}'

        # collect object weights if requested
        weights = truth.get_weights() if args.weighted else None

        # convert to a plain dict representation, either soft (1:*) or hard (1:1)
        truth = truth.hard(True, use_set=True) if args.hard else truth.soft(True)

        # read clustering and convert to basic soft table
        if not os.path.exists(args.pred):
            raise IOError(f'Path does not exist: {args.pred}')
        elif os.path.isfile(args.pred):
            print('Processing as a single file input')
            pred_files = [args.pred]
        else:
            print('Processing as directory input')
            pred_files = [fn for fn in glob.glob(os.path.join(args.pred, '*.mcl'))]

        clusterings = []
        for fn in pred_files:
            print(f'Read clustering from: {fn}')
            cl = tt.read_mcl(fn)
            assert len(cl) > 0, f'Clustering contains no assignments: {fn}'
            clusterings.append([fn, cl.soft(True)])

    except Exception as e:
        print(e)
        sys.exit(1)

    manger = Manager()
    shared_queue = manger.Queue()
    input_args = []
    for cl in clusterings:
        info = OrderedDict({
            'truth_table': args.truth,
            'prediction': cl[0],
            'use_weights': args.weighted,
            'use_hard_truth': args.hard})
        input_args.append([shared_queue, info, truth, cl[1], weights])

    with Pool(args.ncpu) as pool:
        pool.starmap(bcubed_fscore, input_args)

    result_table = []
    while not shared_queue.empty():
        result_table.append(shared_queue.get())

    # write out as a single line csv file
    pd.DataFrame(result_table).to_csv(args.output, index=False, float_format='%.8f')
