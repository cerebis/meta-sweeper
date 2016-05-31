from collections import OrderedDict, Iterable, Counter
import pandas as pd
import numpy as np
import copy
import numpy
import yaml

YAML_WIDTH = 1000


"""
 Dump OrderedDict like a regular dict. Will not reserialize ordered in this case.
"""
def order_rep(dumper, data):
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', data.items(), flow_style=False)

yaml.add_representer(OrderedDict, order_rep)


class Assignment:

    def __init__(self, mapping, weight=1):
        self.mapping = mapping
        self.weight = weight

    def get_classes(self):
        return self.mapping.keys()

    def get_primary_class(self):
        _v = sorted(self.mapping, key=self.mapping.get, reverse=True)
        return _v[0]

    def mean_proportion(self):
        return numpy.mean(self.mapping.values())

    def num_classes(self):
        return len(self.get_classes())

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'weight={0} mapping={1}'.format(self.weight, str(self.mapping))

"""
    Dump Assignment objects as a mapping of member variables
"""
def assignment_rep(dumper, data):
    memb_map = {'weight': data.weight, 'mapping': data.mapping}
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', memb_map, flow_style=True)

yaml.add_representer(Assignment, assignment_rep)


class TruthTable(object):
    """
    Class which represents a truth table in the general sense.
    This can be either a prediction or the ground truth.

    The object stores object:class assignments along with a value
    representing the support/strength of that assignment. The intention
    is that this support can be used to convert multi-class assignments
    to a single most significant assignment. It is up to the user
    supply to choose a useful support value.
    """
    def __init__(self):
        self.asgn_dict = {}
        self.label_map = {}
        self.label_count = Counter()

    def __len__(self):
        """
        Length of truth table is equal to the number of assignment classes.
        :return: number of assignment classes
        """
        return len(self.asgn_dict.keys())

    def num_symbols(self):
        return len(self.label_count)

    def num_assignments(self):
        return sum(self.label_count.values())

    def num_objects(self):
        return len(self.asgn_dict.keys())

    def degeneracy(self, lengths=None):
        nobj =  self.num_objects()
        if nobj == 0:
            return None
        if lengths:
            s = 0
            l = 0
            for k, v in self.asgn_dict.iteritems():
                s += len(v) * lengths[k]
                l += lengths[k]
            return s/float(l)
        else:
            return float(self.num_assignments()) / nobj

    def invert(self):
        cl = {}
        for oi, clist in self.asgn_dict.iteritems():
            for ci in clist.mapping:
                if ci not in cl:
                    cl[ci] = set()
                cl[ci].add(oi)
        return cl

    def mean_overlap(self, lengths=None):
        cl = self.invert()
        ckeys = cl.keys()
        nkeys = len(ckeys)
        if nkeys == 0:
            return None

        ovl = 0.0
        if lengths:
            for i in xrange(nkeys):
                for j in xrange(i+1, nkeys):
                    int_cls = cl[ckeys[i]] & cl[ckeys[j]]
                    sint = 0
                    for ic in int_cls:
                        sint += lengths[ic]
                    uni_cls = cl[ckeys[i]] | cl[ckeys[j]]
                    sovl = 0
                    for ic in uni_cls:
                        sovl += lengths[ic]
                    ovl += sint / float(sovl)
            return ovl / (2*nkeys)
        else:
            for i in xrange(nkeys):
                for j in xrange(i+1, nkeys):
                    int_cls = len(cl[ckeys[i]] & cl[ckeys[j]])
                    uni_cls = len(cl[ckeys[i]] | cl[ckeys[j]])
                    ovl += int_cls / float(uni_cls)
            return ovl / (2*nkeys)

    def overlaps(self, lengths=None):
        cl = self.invert()
        ckeys = cl.keys()
        nkeys = len(ckeys)
        if nkeys == 0:
            return None

        print nkeys
        ovl = np.zeros((nkeys,nkeys))

        if lengths:
            for i in xrange(nkeys):
                for j in xrange(nkeys):
                    int_cls = cl[ckeys[i]] & cl[ckeys[j]]
                    sint = 0
                    for ic in int_cls:
                        sint += lengths[ic]
                    uni_cls = cl[ckeys[i]] | cl[ckeys[j]]
                    sovl = 0
                    for ic in uni_cls:
                        sovl += lengths[ic]
                    ovl[i,j] = sint / float(sovl)
        else:
            for i in xrange(nkeys):
                for j in xrange(nkeys):
                    int_cls = len(cl[ckeys[i]] & cl[ckeys[j]])
                    uni_cls = len(cl[ckeys[i]] | cl[ckeys[j]])
                    ovl[i,j] = int_cls / float(uni_cls)

        return pd.DataFrame(ovl, index=ckeys, columns=ckeys)

    def print_tally(self):
        n_symbol = self.num_symbols()
        n_assignments = self.num_assignments()
        n_objects = self.num_objects()
        degen_ratio = 100.0 * self.degeneracy()

        print '{0} symbols in table, {1:.0f} assignments of {2} objects ({3:.1f}% degeneracy)'.format(
            n_symbol, n_assignments, n_objects, degen_ratio)

        print 'ext_symb\tint_symb\tcount\tpercentage'
        for ci in sorted(self.label_count, key=self.label_count.get, reverse=True):
            print '{0}\t{1}\t{2}\t{3:5.3f}'.format(ci, self.label_map[ci],
                self.label_count[ci], self.label_count[ci] / float(n_assignments))

    def refresh_counter(self):
        self.label_count = Counter()
        for k, v in self.asgn_dict.iteritems():
            self.label_count.update(v.mapping)

    def cluster_extents(self, obj_weights):
        cl_map = self.invert()
        extents = {}
        for ci in cl_map:
            extents[ci] = np.sum([obj_weights[oi] for oi in cl_map[ci]])
        return extents

    def cluster_N50(self, obj_weights):
        cl_map = self.invert()
        n50 = {}
        for ci in cl_map:
            desc_len = sorted([obj_weights[oi] for oi in cl_map[ci]], reverse=True)
            sum_len = np.sum(desc_len)
            sum_oi = 0
            for i in xrange(len(desc_len)):
                sum_oi += desc_len[i]
                if sum_oi > 0.5*sum_len:
                    break
            n50[ci] = desc_len[i]
        return n50

    def filter_extent(self, min_proportion, obj_weights):
        #print '##filter_started_with {0}'.format(len(self.label_count.keys()))
        # make a inverted mapping, to build the deletion collection
        cl_map = self.invert()
        cl_keys = cl_map.keys()
        sum_weight = float(np.sum(obj_weights.values()))

        for ci in cl_keys:
            rw = np.sum([obj_weights[oi] for oi in cl_map[ci]])/sum_weight
            #print ci, cl_map[ci],rw
            if rw < min_proportion:
                for oi in self.asgn_dict.keys():
                    if ci in self.asgn_dict[oi].mapping:
                        del self.asgn_dict[oi].mapping[ci]
                        if len(self.asgn_dict[oi].mapping) == 0:
                            del self.asgn_dict[oi]
                del self.label_map[ci]

        self.refresh_counter()
        #print '##filter_finished_with {0}'.format(len(self.label_count.keys()))

    def filter_class(self, min_proportion):
        """
        Remove classes which represent less than a threshold proportion of all
        objects in the table. This can be used to address problems of scale wrt
        algorithm performance.
        """
        print '##filter_started_with {0}'.format(len(self.label_count.keys()))
        n_obj = float(sum(self.label_count.values()))
        #counter = Counter()
        for ci in sorted(self.label_count, key=self.label_count.get, reverse=True):
            if self.label_count[ci] / n_obj < min_proportion:
                # remove assignments to class
                for k in self.asgn_dict.keys():
                    if ci in self.asgn_dict[k].mapping:
                        del self.asgn_dict[k].mapping[ci]
                        if len(self.asgn_dict[k].mapping) == 0:
                            del self.asgn_dict[k]
                # remove class
                del self.label_map[ci]

        #for k, v in self.asgn_dict.iteritems():
        #    counter.update(v)
        #self.label_count = counter
        self.refresh_counter()

        print '##filter_finished_with {0}'.format(len(self.label_count.keys()))

    def get_weights(self):
        _w = {}
        for k, asgn in self.asgn_dict.iteritems():
            _w[k] = asgn.weight
        return _w

    def soft(self, universal=False):
        """
        :return plain dictionary with degenerate classification
        """
        _s = OrderedDict()
        _keys = sorted(self.asgn_dict.keys())
        for k in _keys:
            clz = self.asgn_dict[k].mapping.keys()
            if universal:
                # relabel with universal symbols if requested
                clz = [self.label_map[ci] for ci in clz]
            _s[k] = clz
        return _s

    def hard(self, universal=False):
        """
        Convert TT to a plain dictionary with the single most significant classification only.
        In the case of a tie, no effort is made to be uniformly random in the case of a tie and
        dependent on the behaviour of sort.
        :return plain dict with only one class->cluster mapping.
        """
        _s = OrderedDict()
        _keys = sorted(self.asgn_dict.keys())
        for _k in _keys:
            pc = self.asgn_dict[_k].get_primary_class()
            if universal:
                pc = self.label_map[pc]
            _s[_k] = [pc]
        return _s

    def get(self, key):
        return self.asgn_dict.get(key)

    def put(self, key, value, weight=None):
        self.asgn_dict[key] = Assignment(value)
        if weight:
            self.asgn_dict[key].weight = weight

    def update_from_yaml(self, yd):
        """
        Initialise a TruthTable from a generic dict of dicts object stored
        in YAML format. We chose to avoid a custom class here for inter-codebase
        portability.

        :param yd: generic yaml object, dict of dicts
        :return:
        """
        for k, v in yd.iteritems():
            self.asgn_dict[k] = Assignment(v['mapping'], v['weight'])
            self.label_count.update(v['mapping'].keys())
        labels = sorted(self.label_count.keys())
        self.label_map = dict((l, n) for n, l in enumerate(labels, 1))

    def update(self, dt, weights=None, min_score=0):
        """
        Initialise the assignment dictionary and also generate a mapping of
        class symbol to the positive integers. We can use this as a universal
        symbol basis.

        :param dt: the dictionary to initialise from
        """
        all_asgn = 0
        filt_obj = 0
        filt_asgn = 0
        all_obj = len(dt)

        for k, v in dt.iteritems():
            v_filt = dict((kv, vv) for kv, vv in v.iteritems() if int(vv) >= min_score)
            filt_asgn += len(v) - len(v_filt)
            all_asgn += len(v)
            if len(v_filt) == 0:
                filt_obj += 1
                continue
            self.asgn_dict[str(k)] = Assignment(v_filt)
            if weights:
                self.asgn_dict[str(k)].weight = weights[k]
            self.label_count.update(v_filt.keys())

        if filt_asgn > 0:
            print 'Filtered {0}/{1} assignments and {2}/{3} objects below minimum score {4}'.format(
                filt_asgn, all_asgn, filt_obj, all_obj, min_score)

        labels = sorted(self.label_count.keys())
        self.label_map = dict((l, n) for n, l in enumerate(labels, 1))

    def to_vector(self):
        vec = {}
        keys = sorted(self.asgn_dict.keys())
        hd = self.hard()
        for k in keys:
            vec[k] = hd[k]
        return vec

    @staticmethod
    def _write_dict(dt, pathname):
        with open(pathname, 'w') as h_out:
            for obj, clz in dt.iteritems():
                h_out.write('{0} {1}\n'.format(obj, ' '.join(sorted(clz))))

    def write(self, pathname):
        """
        Write the full table in a simple (no custom classes) YAML format"
        :pathname the output path
        """
        with open(pathname, 'w') as h_out:
            yaml.dump(self.asgn_dict, h_out, default_flow_style=False, width=YAML_WIDTH)

    def write_soft(self, pathname):
        """
        Write a plain dictionary representation in YAML format.
        :pathname the output path
        """
        TruthTable._write_dict(self.soft(), pathname)

    def write_hard(self, pathname):
        """
        Write a plain dictionary representation of only the most significant
        object:class assignments.
        :pathname the output path
        """
        TruthTable._write_dict(self.hard(), pathname)



def read_truth(pathname):
    """
    Read a TruthTable in YAML format
    :param pathname: path to truth table
    :param min_score: ignore objects with scores/weights/lengths below the given threshold.
    :return: truth table
    """
    with open(pathname, 'r') as h_in:
        tt = TruthTable()
        tt.update_from_yaml(yaml.load(h_in))
        return tt


def read_mcl(pathname):
    """
    Read a MCL solution file converting this to a TruthTable
    :param pathname: mcl output file
    :return: truth table
    """

    with open(pathname, 'r') as h_in:
        # read the MCL file, which lists all members of a class on a single line
        # the class ids are implicit, therefore we use line number.
        mcl = {}
        for ci, line in enumerate(h_in, start=1):
            objects = line.rstrip().split()
            for oi in objects:
                if oi not in mcl:
                    mcl[oi] = {}
                mcl[oi][ci] = 1.0  # there are no weights, therefore give them all 1
        # initialise the table
        tt = TruthTable()
        tt.update(mcl)
    return tt


def unique_labels(dt):
    """
    Extract the unique set of class labels used in the truth table
    :param dt: dictionary representation of truth table to analyze
    :return: sorted set of unique labels
    """
    labels = set()
    for v in dt.values():
        if isinstance(v, Iterable):
            labels.update(v)
        else:
            labels.add(v)
    return sorted(labels)


def crosstab(dt1, dt2):
    """
    Cross-tabulate two truth tables on hard clusters.

    :param dt1: first dictionary rep of truth table
    :param dt2: second dictionary rep of truth table
    :return: pandas dataframe
    """
    joined_keys = sorted(set(dt1.keys() + dt2.keys()))

    rows = unique_labels(dt1)
    cols = unique_labels(dt2)
    ctab = pd.DataFrame(0, index=rows, columns=cols)

    for k in joined_keys:
        if k in dt1 and k in dt2:
            i1 = dt1[k]
            i2 = dt2[k]
            ctab.loc[i1, i2] += 1

    return ctab


def simulate_error(tt, p_mut, p_indel, extra_symb=[]):
    """
    Simple method for introducing error in a truth table. This is useful when
    testing clustering metrics (Fm, Vm, Bcubed, etc). By default, the list of possible
    symbols is taken from those already assigned, but a user may provide additional
    symbols. These can provide a useful source of novelty, when for instance
    an object is already assigned to all existing class symbols.

    :param tt: the truth table to add error
    :param p_mut: the probably of a class mutation
    :param p_indel: the probability of deletion or insertion of a class to an object
    :param extra_symb: extra class symbols for inserting
    :return: truth table mutatant
    """
    symbols = list(set(tt.label_count.keys() + extra_symb))
    print symbols
    mut_dict = copy.deepcopy(tt.asgn_dict)

    for o_i in mut_dict.keys():
        others = list(set(symbols) - set(mut_dict[o_i].mapping))
        if numpy.random.uniform() < p_mut:
            if len(others) > 0:
                c_mut = numpy.random.choice(others, 1)[0]
                c_old = numpy.random.choice(mut_dict[o_i].mapping.keys(), 1)[0]

                # retain the weighting from the original to mutated
                # we do this pedantically so its easy to read
                weight = mut_dict[o_i].mapping[c_old]
                mut_dict[o_i].mapping[c_mut] = weight
                del mut_dict[o_i].mapping[c_old]

        if numpy.random.uniform() < p_indel:
            # flip a coin, delete or insert
            if numpy.random.uniform() < 0.5:
                # delete
                c_del = numpy.random.choice(mut_dict[o_i].mapping.keys(), 1)[0]
                del mut_dict[o_i].mapping[c_del]
            elif len(others) > 0:
                # insert from 'others'
                c_add = numpy.random.choice(others, 1)[0]
                num_cl = mut_dict[o_i].num_classes()
                adj_fac = float(num_cl / (num_cl+1.))
                ins_prop = mut_dict[o_i].mean_proportion() * adj_fac
                for k in mut_dict[o_i].mapping:
                    mut_dict[o_i].mapping[k] *= adj_fac
                mut_dict[o_i].mapping[c_add] = ins_prop

    mut_tt = TruthTable()
    mut_tt.update_from_yaml(mut_dict)
    return mut_tt