'''
DedupUMI.py - Deduplicate reads that are coded with a UMI
=========================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

The perpose of this script is to deduplicate BAM files based
on the first mapping co-ordinate and the UMI attached to the read.
It is assumed that the FASTQ files were processed with extract_UMI.py
and thus the UMI is the last word of the read name.

Two reads that have the same start co-ordinate and the same UMI and the same
splicing are assumed to be PCR duplicates.

The following criteria are applied to select the read that will be retained:

1. The read with the lowest nubmer of hits
2. The read with the highest mapping quality

Otherwise a read is chosen at random.

The input file must be sorted.

.. note::
   You need to redirect either logging information or output to a file or turn
   it off logging via -v 0 in order to get a valid sam/bam file.

Options
-------

The script assumes that input and output are in BAM format. Use that
-i/--in-sam and -o/--out-sam options to specifiy that the input/output
are SAM format.

The option --subset can be used to randomly select subsets of the reads
for futher processing.

Usage
-----

    python dedup_umi.py -I infile.bam -O deduped.bam

or

    cat infile.bam | python dedup_umi.py > deduped.bam 2> deduped.log

Command line options
--------------------

'''

import sys
import pysam
import CGAT.Experiment as E
import random
import collections


def breadth_first_search(node, graph):
    '''Returns all nodes reachable from node in graph
    where graph is an adjecency list implemented as a dictionary'''

    found = set()
    queue = list()
    queue.append(node)
    while len(queue) != 0:
        current, queue = queue[-1], queue[:-1]
        found.add(current)
        for child in graph[current]:
            if child not in found:
                queue.append(child)
    return found


def connected_components(graph):
    '''Takes a graph as a dictionary based adjencency list
    and returns a list of sets, where each set is repesents
    a connected component in the graph'''

    found = list()
    components = list()
    for node in graph:
        if node not in found:
            component = breadth_first_search(node, graph)
            found.extend(component)
            components.append(component)
    return components


def edit_dist(first, second):
    ''' returns the edit distance/hamming distances between
    its two arguements '''

    dist = sum([not a == b for a, b in zip(first, second)])
    return dist


def get_adj_list(umis, threshold=1):
    '''Returns dictionary where each entry represents the nodes
    adjecent to the dictionary key'''

    return {umi: [umi2 for umi2 in umis
                  if 0 < edit_dist(umi, umi2) <= threshold]
            for umi in umis}


def remove_umis(adj_list, cluster, nodes):
    '''removes the specified nodes from the cluster and returns
    the remaining nodes '''

    # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
    nodes_to_remove = set([node
                           for x in nodes
                           for node in adj_list[x]] + nodes)

    return cluster - nodes_to_remove


def get_best(cluster, adj_list, counts):
    '''Finds the nodes that best explain the cluster by successively
    removing more and more nodes until the cluster is explained.
    Nodes are removed in order of their count. '''

    if len(cluster) == 1:
        return list(cluster)

    sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                          reverse=True)

    for i in range(len(sorted_nodes) - 1):
        if len(remove_umis(adj_list, cluster, sorted_nodes[:i+1])) == 0:
            return sorted_nodes[:i+1]


def cluster_and_reduce(bundle, threshold=1):
    ''' Recieves a bundle of reads, clusters them on edit distance
    and selects connected clusters using the edit distance threshold.
    Within each cluster the UMIs with the highest reads are progressively
    selected, until the whole cluster is explained. Yeilds the selected
    reads '''

    umis = bundle.keys()
    counts = {umi: bundle[umi]["count"] for umi in umis}

    adj_list = get_adj_list(umis, threshold)
    clusters = connected_components(adj_list)

    for cluster in clusters:
        for umi in get_best(cluster, adj_list, counts):
            yield bundle[umi]["read"]


def get_bundles(insam, ignore_umi=False, subset=None, paired=False,
                chrom=None, spliced=False, soft_clip_threshold=0):
    ''' Returns a dictionary of dictionaries, representing the unique reads at
    a position/spliced/strand combination. The key to the top level dictionary
    is a umi. Each dictionary contains a "read" entry with the best read, and a
    count entry with the number of reads with that position/spliced/strand/umi
    combination'''
    
    last_pos = 0
    last_chr = ""
    reads_dict = collections.defaultdict(
        lambda: collections.defaultdict(
            lambda: collections.defaultdict(dict)))
    read_counts = collections.defaultdict(
        lambda: collections.defaultdict(dict))
    
    if chrom:
        inreads = insam.fetch(reference=chrom)
    else:
        inreads = insam.fetch()
        
    for read in inreads:

        if subset:
            if random.random() >= subset:
                continue

        if read.is_unmapped:
            continue

        if read.mate_is_unmapped and paired:
            continue

        if read.is_read2:
            continue

        is_spliced = False

        if read.is_reverse:
            pos = read.aend
            if read.cigar[-1][0] == 4:
                pos = pos + read.cigar[-1][1]
            start = read.pos

            if ('N' in read.cigarstring or
                (read.cigar[0][0] == 4 and
                 read.cigar[0][1] > soft_clip_threshold)):
                is_spliced = True
        else:
            pos = read.pos
            if read.cigar[0][0] == 4:
                pos = pos - read.cigar[0][1]
            start = pos

            if ('N' in read.cigarstring or
                (read.cigar[-1][0] == 4 and
                 read.cigar[-1][1] > soft_clip_threshold)):
                is_spliced = True

        if start > (last_pos+1000) and not read.tid == last_chr:

            out_keys = [x for x in reads_dict.keys() if x <= start-1000]
            
            for p in out_keys:
                for bundle in reads_dict[p].itervalues():
                    yield bundle
                del reads_dict[p]
                del read_counts[p]
            
            last_pos = start
            last_chr = read.tid

        if ignore_umi:
            umi = ""
        else:
            umi = read.qname.split("_")[-1]

        key = (read.is_reverse, spliced & is_spliced, paired*read.tlen)

        try:
            reads_dict[pos][key][umi]["count"] += 1
        except KeyError:
            reads_dict[pos][key][umi]["read"] = read
            reads_dict[pos][key][umi]["count"] = 1
            read_counts[pos][key][umi] = 0
        else:
            if reads_dict[pos][key][umi]["read"].mapq > read.mapq:
                continue
                
            if reads_dict[pos][key][umi]["read"].mapq < read.mapq:
                reads_dict[pos][key][umi]["read"] = read
                read_counts[pos][key][umi] = 0
                continue

            if reads_dict[pos][key][umi]["read"].opt("NH") < read.opt("NH"):
                continue
            elif reads_dict[pos][key][umi]["read"].opt("NH") > read.opt("NH"):
                reads_dict[pos][key][umi]["read"] = read
                read_counts[pos][key][umi] = 0

            read_counts[pos][key][umi] += 1
            prob = 1.0/read_counts[pos][key][umi]

            if random.random() < prob:
                reads_dict[pos][key][umi]["read"] = read

    #yeild remaining bundles
#    for p in reads_dict:
#                for bundle in reads_dict[p].itervalues():
#                    yield bundle
               
@profile
def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-i", "--in-sam", dest="in_sam", action="store_true",
                      help="Input file is in sam format", default=False)
    parser.add_option("-o", "--out-sam", dest="out_sam", action="store_true",
                      help="Output alignments in sam format", default=False)
    parser.add_option("--ignore-umi", dest="ignore_umi", action="store_true",
                      help="Ignore UMI and dedup only on position",
                      default=False)
    parser.add_option("--subset", dest="subset", type="string",
                      help="Use only a fraction of reads, specified by subset",
                      default=1.1)
    parser.add_option("--spliced-is-unique", dest="spliced",
                      action="store_true",
                      help="Treat a spliced read as different to an unspliced"
                           " one",
                      default=False)
    parser.add_option("--soft-clip-threshold", dest="soft",
                      type="float",
                      help="number of bases clipped from 5' end before"
                           "read os counted as spliced",
                      default=4)
    parser.add_option("--cluster-umis", dest="cluster_umis",
                      action="store_true",
                      help="Select best reads by clustering umis and removing"
                           "potential sequencing errors", default=False)
    parser.add_option("--edit-distance-theshold", dest="threshold",
                      type="int",
                      help="Edit distance theshold at which to join two UMIs"
                           "when clustering", default=1)
    parser.add_option("--chrom", dest="chrom", type="string",
                      help="Restrict to one chromosome",
                      default=None)
    parser.add_option("--paired", dest="paired", action="store_true",
                      default=False,
                      help="Use second-in-pair position when deduping")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        in_name = "-"

    if options.stdout != sys.stdout:
        out_name = options.stdout.name
        options.stdout.close()
    else:
        out_name = "-"

    if options.in_sam:
        in_mode = "r"
    else:
        in_mode = "rb"

    if options.out_sam:
        out_mode = "w"
    else:
        out_mode = "wb"

    infile = pysam.Samfile(in_name, in_mode)
    outfile = pysam.Samfile(out_name, out_mode,
                            template=infile)

    nInput, nOutput = 0, 0

    for bundle in get_bundles(infile,
                              ignore_umi=options.ignore_umi,
                              subset=float(options.subset),
                              paired=options.paired,
                              chrom=options.chrom,
                              spliced=options.spliced,
                              soft_clip_threshold=options.soft):

      
        nOutput +=1
        nInput += sum([bundle[umi]["count"] for umi in bundle])

        if nOutput % 10000 == 0:
            E.debug("Outputted %i" % nOutput)

        if nInput % 1000000 == 0:
            E.debug("Read %i input reads" % nInput)

        if options.ignore_umi or not options.cluster_umis:
            for umi in bundle:
         
                outfile.write(bundle[umi]["read"])
        else:
            for read in cluster_and_reduce(bundle, options.threshold):
                outfile.write(read)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
