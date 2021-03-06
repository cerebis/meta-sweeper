#!/usr/bin/env nextflow
/**
 * Aggregation of results from CCC (HiC/3C) clustering
 *
 * After chrcc-sweep.nf and chrcc-cluster.nf, this workflow to analyse the performance of
 * each point in the iteration and produce a single file compiling the results. The
 * results are formatted as an associative collection (Java map, Python dict) in YAML format.
 *
 * Usage: chrcc-aggregate.nf [--debug]
 */
/*
 * meta-sweeper - for performing parametric sweeps of simulated
 * metagenomic sequencing experiments.
 * Copyright (C) 2016 "Matthew Z DeMaere"
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
import MetaSweeper

MetaSweeper ms = MetaSweeper.fromFile(new File('chrcc.yaml'))

sweep = MetaSweeper.createSweep()

a_stats = ms.keyedFrom(sweep, file("$ms.options.output/*asmstat"))

// collect grpahs and their statistics
g_stats = ms.keyedFrom(sweep, file("$ms.options.output/*.graphml"))
        .flatCross(ms.keyedFrom(sweep, file("$ms.options.output/*.gstat")))
        .flatCross(ms.keyedFrom(sweep, file("$ms.options.output/*.geigh")))

// collect clusterings and their statistics
cl_stats = ms.keyedFrom(sweep, file("$ms.options.output/*.cl"))
        .flatCross(ms.keyedFrom(sweep, file("$ms.options.output/*.bc")))

// join the three channels together at the appropriate sweep depths
stat_sweep = sweep.joinChannels(a_stats, g_stats, 3)
stat_sweep = sweep.joinChannels(stat_sweep, cl_stats, 4)

process Aggregate {

    input:
    set key, file('asmstat'), file('graph'), file('gstat'), file('geigh'), file('cl'), file('bc') from stat_sweep

    output:
    set key, stdout into all_stats

    script:
    if (params.debug) {
        """
        echo $key
        """
    }
    else {
        """
        aggregate.py --fmt yaml asmstat gstat geigh bc
        """
    }
}

if (params.debug) {
    println "Debug does not attempt to serialize"
}
else {
    parser = ms.getYamlParser()
    all_stats = all_stats.map { [params: it.getKey().map()] + parser.load(it.dropKey()) }
    fout = file("$ms.options.output/all_stats.yaml")
    fout.write(parser.dump(all_stats.toList().get()))
    fout << '\n'
}
