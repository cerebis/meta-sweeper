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
from Bio import SeqIO
import truthtable as ttable
import argparse
import os


parser = argparse.ArgumentParser(description='Extract clustered sequences')
parser.add_argument('-v', dest='verbose', default=False, action='store_true', help='Verbose standard output')
parser.add_argument('-f', '--force', default=False, action='store_true', help='Force overwriting of files')
parser.add_argument('--cid-list', nargs='*', help='Specific cluster IDs')
parser.add_argument('--clustering', required=True, metavar='CLUSTERING', help='MCL format clustering file')
parser.add_argument('--fasta', required=True, metavar='FASTA', help='Fasta sequences supporting clustering')
parser.add_argument('-o', '--out-dir', help='Output directory')
args = parser.parse_args()

if args.out_dir:
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    elif not os.path.isdir(args.out_dir):
        raise IOError('Output path exists and is not a directory')

if args.verbose:
    print 'Getting sequence index'
seq_index = SeqIO.index(args.fasta, 'fasta')

if args.verbose:
    print 'Reading clustering'
tt = ttable.read_mcl(args.clustering)
cl2seq = tt.invert()


if args.cid_list:
    cid_list = args.cid_list
else:
    cid_list = [ci for ci in cl2seq]

for ci in cid_list:

    if args.verbose:
        print 'Collecting sequences for cluster {0}'.format(ci)

    seqs = []
    for si in cl2seq[int(ci)]:
        seqs.append(seq_index[si])

    if args.verbose:
        print 'Writing {0} sequences for cluster {1}'.format(len(seqs), ci)

    opath = os.path.join(args.out_dir, 'cl{0}.fasta'.format(ci))
    if not args.force and os.path.exists(opath):
        raise IOError('Path {0} already exists. Use --force to overwrite destination files'.format(opath))
    SeqIO.write(seqs, opath, 'fasta')
