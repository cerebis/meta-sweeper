#!/usr/bin/env nextflow
/**
 * Time series based deconvolution workflow
 *
 * Test metagenomic deconvolution over a range of simulated strain evolution parameters
 *
 * Usage: timeseries-deconvolute.nf [--debug]
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

MetaSweeper ms = MetaSweeper.fromFile(new File('timeseries.yaml'))

//
// Generate phylogenetic trees for each clade within each community
//

def sweep = MetaSweeper.createSweep()
        .withVariable('seed', ms.variables.seed)
        .withVariable('clade', ms.variables.community.clades)
        .describe('Tree Generation')

// channel composed of the permutation of variables
gen_in = sweep.permutedChannel()

process TreeGen {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, seed, clade from gen_in

    output:
    set key, file("${key}.nwk"), seed, clade into tree_out

    script:
    if (params.debug) {
        """
        echo "$key ${clade.tree}" > "${key}.nwk"
        """
    }
    else {
        if (clade.isDefined()) {
            """
            echo "${clade.getDefined()}" | tree_scaler.py --max-height 0.1 - ${key}.nwk
            """
        }
        else {
            // presently we always assume birth-death
            assert clade.isSupportedAlgorithm() : 'Only birth_death is currently supported'
            """
            tree_generator.py --seed $seed --prefix ${clade.prefix} --suppress-rooting --mode random \
                --max-height 0.1 --birth-rate ${clade.tree.birth_rate} --death-rate ${clade.tree.death_rate} \
                --format newick --num-taxa ${clade.ntaxa} ${key}.nwk
            """
        }
    }
}


//
// Generate evolved sequences for each clade from each community
//

(tree_out, evo_in) = tree_out.into(2)

// add variation on alpha
sweep.withVariable('alpha', ms.variables.alpha)
        .describe('Evolve Clades')

// extend the channel to include new parameter
evo_in = sweep.extendChannel(evo_in, 'alpha')

process Evolve {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, tree_file, seed, clade, alpha from evo_in

    output:
    set key, file("${key}.evo.fa"), seed, clade into evo_out

    script:
    if (params.debug) {
        """
        echo $key > "${key}.evo.fa"
        """
    }
    else {
        """
        scale_tree.py -a $alpha $tree_file scaled_tree
        \$EXT_BIN/sgevolver/sgEvolver --indel-freq=${ms.options.evo.indel_freq} \
            --small-ht-freq=${ms.options.evo.small_ht_freq} \
            --large-ht-freq=${ms.options.evo.large_ht_freq} \
            --inversion-freq=${ms.options.evo.inversion_freq} \
            --random-seed=$seed scaled_tree \
             $clade.ancestor $clade.donor "${key}.evo.aln" "${key}.evo.fa"
        strip_semis.sh "${key}.evo.fa"
        """
    }

}


//
// Merge evolved sequences from the clades into whole communities
//

(evo_out, merge_seq_in) = evo_out.into(2)

// group by a reduced key that is only the random seed and alpha
merge_seq_in = merge_seq_in.groupBy { it.getKey().selectedKey('seed', 'alpha') }
// convert the resulting map of sweep point results into table format and sort by file name
        .flatMap { it.collect { k, v -> [k, v.collect { vi -> vi[1] }.toSorted { a, b -> a.name <=> b.name }] } }

process MergeClades {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('clade_seq') from merge_seq_in

    output:
    set key, file("${key}.community.fa") into merge_seq_out

    script:
    if (params.debug) {
        """
        echo $key > "${key}.community.fa"
        """
    }
    else {
        """
        cat clade_seq* >> ${key}.community.fa
        """
    }
}

//
// Generate abundance profiles for each clade within each community
//
(merge_seq_out, prof_in) = merge_seq_out.into(2)

community = Channel.value(ms.variables['community'])

process ProfileGen {
    publishDir ms.options.output, mode: 'copy', overwrite: 'true'

    input:
    set key, file('community.fa') from prof_in
    val community

    output:
    set key, file("${key}.prf") into prof_out

    script:
    def mu = community.profile.mu
    def sigma = community.profile.sigma
    if (params.debug) {
        """
        echo "$key $mu $sigma" > "${key}.prf"
        """
    }
    else {
        """
        profile_generator.py --seed ${key['seed']} --dist lognormal --lognorm-mu $mu \
            --lognorm-sigma $sigma community.fa ${key}.prf
        """
    }
}


//
// Merge and pair community sequences and profiles
//

(merge_seq_out, seq_prof) = merge_seq_out.into(2)
(merge_prof_out, tmp) = prof_out.into(2)

// select just the community sequences
seq_prof = seq_prof.map { it.pick(1) }
// combine with their respective profiles, then flatten and simplify the rows
        .phase(tmp).map { it = it.flatten(); it.pick(1, 3) }

//
// Generate shotgun sequencing reads for for each whole community
//
(seq_prof, wgs_in) = seq_prof.into(2)

// Add wgs coverage to sweep
sweep.withVariable('xfold', ms.variables.xfold)
        .describe('WGS Read Generation')

// extend the channel
wgs_in = sweep.extendChannel(wgs_in, 'xfold')

process WGS_Reads {
    publishDir ms.options.output, mode: 'copy', overwrite: false

    input:
    set key, file(comm_seq), file(comm_prof), xfold from wgs_in
    val community

    output:
    set key, file("${key}.wgs.*.r1.fq.gz"), file("${key}.wgs.*.r2.fq.gz"), file("${key}.cov") into wgs_out

    script:
    def mu = community.profile.mu
    def sigma = community.profile.sigma
    if (params.debug) {
        """
        for ((n=1; n<=$ms.options.num_samples; n++))
        do
            echo "metaART.py -C gzip --profile $comm_prof -z $ms.options.num_samples -M $xfold -S ${key['seed']} \
                    -s $ms.options.wgs.ins_std -m $ms.options.wgs.ins_len -l $ms.options.wgs.read_len
                    --coverage-out ${key}.cov -n ${key}.wgs $comm_seq ." > ${key}.wgs.\$n.r1.fq.gz

            echo "metaART.py -C gzip --profile $comm_prof -z $ms.options.num_samples -M $xfold -S ${key['seed']} \
                    -s $ms.options.wgs.ins_std -m $ms.options.wgs.ins_len -l $ms.options.wgs.read_len \
                    --coverage-out ${key}.cov -n ${key}.wgs $comm_seq ." > ${key}.wgs.\$n.r2.fq.gz
        done

        touch ${key}.cov
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/art:\$PATH
        metaART.py -C gzip -z $ms.options.num_samples -M $xfold -S ${key['seed']} --dist lognormal --lognorm-mu $mu --lognorm-sigma $sigma \
                -s $ms.options.wgs.ins_std -m $ms.options.wgs.ins_len -l $ms.options.wgs.read_len \
                --coverage-out ${key}.cov -n ${key}.wgs $comm_seq .
        """
    }
}

//
// Map WGS reads to reference sequences
//

// ancestral sequence for community
ancestor_in = Channel.value(file(ms.variables.community.clades[0].ancestor))

(wgs_out, map_in) = wgs_out.into(2)

        // ref, R1, R2
map_in = map_in.map{ it.pick(1, 2) }
        // pair the R1s and R2s -- might move this to a method in MetaSweeper
        // then flatten nested lists to become one row per read-pair
        .flatMap{ GroovyCollections.transpose(it[1].sort(), it[2].sort())
                    // extend the key to include sample number, extracting it from read-pair file names
                    .collect{ pair ->
                        def nsamp = pair.collect { ri -> (ri =~ /wgs\.(\d+)\.r[12].fq.*$/)[0][1] }
                        assert nsamp.size() == 2 && nsamp[0] == nsamp[1] : 'Error: read-pairs do not share the same sample index'
                        [sweep.extendKey(it.getKey(), 'nsamp', nsamp[0]), *pair] }
            }

process WGSMap {
    publishDir ms.options.output, mode: 'copy', overwrite: false

    input:
    set key, file(reads1), file(reads2) from map_in
    file(anc) from ancestor_in

    output:
    set key, file("${key}.wgs2ref.bam"), file(reads1), file(reads2) into map_out

    script:
    if (params.debug) {
        """
        if [ ! -e $anc ]; then echo "no ancestral sequence found"; exit 1; fi
        echo $key > ${key}.wgs2ref.bam
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/a5/bin:\$PATH
        if [ ! -e "${anc}.bwt" ]
        then
            bwa index $anc
        fi
        bwa mem -t 1 $anc $reads1 $reads2 | samtools view -bS - | samtools sort -l 9 - ${key}.wgs2ref
        """
    }
}

//
// Deconvolve the SNVs into strain genotypes
//

// ancestral sequence for community
ancestor_in = Channel.value(file(ms.variables.community.clades[0].ancestor))

(map_out, deconv_in) = map_out.into(2)
            // just the bam files
deconv_in = deconv_in.map { it.pick(1) }
                // remove the sample number from key
                .map { [it.getKey().popLevels(1), *it.dropKey()] }
                // group each samples time-series bams on the new reduced key and sort by base filename
                .groupTuple(sort: {it.name})

process Deconvolve {
    publishDir ms.options.output, mode: 'copy', overwrite: false

    input:
    set key, file('tp*.bam') from deconv_in
    file(anc) from ancestor_in

    output:
    set key, file("${key}.decon.csv"), file("${key}.snv_file.data.R"), file("${key}.strains.tre") into deconv_out

    script:
    if (params.debug) {
        """
        if [ ! -e $anc]; then exit 1; fi
        echo $key > ${key}.decon.csv
        echo $key > ${key}.snv_file.data.R
        echo $key > ${key}.strains.tre
        """
    }
    else {
        """
        export PATH=\$EXT_BIN/lofreq_star:\$PATH
        snvbpnmft.py . $anc *.bam
        mv decon.csv ${key}.decon.csv
        mv elbos.csv ${key}.elbos.csv
        mv snv_file.data.R ${key}.snv_file.data.R
        java -Xmx1000m -jar \$EXT_BIN/beast/beast.jar beast.xml
        java -jar \$EXT_BIN/treeanno/treeannotator.jar -burnin 1000 -heights mean aln.trees ${key}.strains.tre
        """
    }
}

//
// Record the true strain genotypes
//

(seq_prof, truth_in) = seq_prof.into(2)

// ancestral sequence for community
ancestor_in = Channel.value(file(ms.variables.community.clades[0].ancestor))

// just the community reference sequence
truth_in = truth_in.map{ it.pick(1) }

process Truth {
    publishDir ms.options.output, mode: 'copy', overwrite: false

    input:
    set key, file(ref) from truth_in
    file(anc) from ancestor_in

    output:
    set key, file("${key}.truth.tsv") into truth_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.truth.tsv
        """
    }
    else {
        """
        strain_truth.py --mauve-path=\$MAUVEPATH -o ${key}.truth.tsv $ref $anc
        """
    }
}

//
// Measure accuracy of strain genotypes
//

// join truth and deconv outputs at sweep depth of 2.
accuracy_in = sweep.joinChannels(truth_out, deconv_out, 2)

process Accuracy {
    publishDir ms.options.output, mode: 'copy', overwrite: false

    input:
    set key, file(truthfile), file(snvbpnmf), file(snv_file), file(tree_file) from accuracy_in

    output:
    file("${key}.truth.report.txt") into accuracy_out

    script:
    if (params.debug) {
        """
        echo $key > ${key}.truth.report.txt
        """
    }
    else {
        """
        measure_accuracy.py --bpnmf=${snvbpnmf} --truth=${truthfile} --sites=${snv_file} > ${key}.truth.report.txt
        """
    }
}
