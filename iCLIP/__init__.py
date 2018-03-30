''' This a modeule that holds functions and classes useful for analysing iCLIP data '''

from counting import count_intervals, count_transcript
from utils import spread, rand_apply, randomiseSites, TranscriptCoordInterconverter
from meta import meta_gene, processing_index, get_binding_matrix, quantile_row_norm,  sum_row_norm, compress_matrix
from kmers import pentamer_enrichment, pentamer_frequency
from distance import calcAverageDistance, findMinDistance, corr_profile
from clusters import Ph, fdr, get_crosslink_fdr_by_randomisation
from getters import make_getter, countChr

import random
