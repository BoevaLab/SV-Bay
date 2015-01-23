import logging
import os
import pysam
import resource
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import bisect
import yaml

import fragment
import utils

logger = logging.getLogger('main_logger')

class ChrFragments(object):

    # Attributes
    # config
    # lists of fragments: ff, rf, fr, rr
    # chr
    # calc_hist - True for first chrom from config, False otherwise
    # collected stats: flag_direction, smallest_normal, biggest_normal, ...
    # mean, median, R, epsilon, mu, sigma, a, b, p
    # lists of abnormal fragments: ff_abn, rf_abn, fr_abn, rr_abn
    # num_abnormal
    # normal_fragments

    def __init__(self, sam_file, config):
        self.config = config
        self.translocations = []
        self.ff = []
        self.fr = []
        self.rf = []
        self.rr = []
        self.__LoadSam(sam_file)
        # Chromosome name (same for all fragments except for translocation, pick from one)
        if len(self.fr) > len(self.rf):
            self.chr = self.fr[0].first_read_chr
        else:
            self.chr = self.rf[0].first_read_chr
        self.calc_hist = False
        if self.chr == config['chromosomes'][0]:
            self.calc_hist = True

    def SaveTranslocationsToTmp(self):
        trans_temp_fname = '/tmp/trans_' + self.chr + '.tmp'
        trans_temp_file = open(trans_temp_fname, 'w')
        for tr in self.translocations:
            trans_temp_file.write(tr.to_string())
            trans_temp_file.write('\n')
        trans_temp_file.close()

    # Collect different stats for this chr
    def CollectStats(self):
        # Build confidence interval. Take fragment lengths for most popular type (it may be fr or rf)
        if len(self.fr) > len(self.rf):
            fragment_lengths = [frag.length for frag in self.fr]
            self.flag_direction='fr'
        else:
            fragment_lengths = [frag.length for frag in self.rf]
            self.flag_direction='rf'
        
        fragment_lengths = sorted(fragment_lengths)
        normal_lengths = fragment_lengths[0:bisect.bisect_left(fragment_lengths,10000)]
        h = int(round(self.config['alpha'] * len(fragment_lengths))) # Number of segments that should be excluded from each side
        normal_lengths = normal_lengths[int(round(h))+1:-int(round(h)) -1]
        #normal_lengths =[i for i in fragment_lengths if fragment_lengths[0]<=i<= fragment_lengths[-1]]
        #normal_lengths =[i for i in fragment_lengths if 900<=i<= 4000] 
        abnormal_lengths = fragment_lengths[0:int(round(h))+1] + fragment_lengths[int(round(h)) -1:-1]
        self.smallest_normal = normal_lengths[0]
        self.biggest_normal = normal_lengths[-1]
        
        # Mean, standart deviation and Radius
        self.mean = sum(normal_lengths) / len(normal_lengths)
        q_sum = 0
        for i in range(0, len(normal_lengths) - 1):
            q_sum += (normal_lengths[i] - self.mean)**2
        sigma = (q_sum / len(normal_lengths))**0.5
        normal_lengths.sort()
        if len(normal_lengths) % 2 == 0:
            n = len(normal_lengths)
            self.median = (normal_lengths[n/2-1]+ normal_lengths[n/2] )/2
        else:
            self.median =normal_lengths[len(normal_lengths)/2] 
        self.R = self.median + 3*sigma/(2**0.5)
        self.epsilon = self.median / 5.0

        if self.calc_hist:
            #Histogram calculation
            logger.info('Calculating histogram...')
            logger.debug('Normal fragments: [' + str(normal_lengths[0]) + ';' + str(normal_lengths[-1]) + ']')
            #bins = set(normal_lengths)
            historam_file = open(self.config['working_dir'] + self.config['length_histogram_file'], 'w')
            n, bins, patches = plt.hist(normal_lengths, len(set(normal_lengths)), normed=1, facecolor='green', alpha=0.75)
            for i in range(len(bins)-1):
                historam_file.write(str(round(bins[i])) + ' ' + str(round(n[i],9)) + '\n')
            historam_file.close()
            (self.mu, self.sigma) = norm.fit(normal_lengths)
            self.a = self.mu - (12**0.5) * 0.5 * self.sigma
            self.b = self.mu + (12**0.5) * 0.5 * self.sigma
            self.p = 1 / (self.b - self.a)

    # Print chr statistics
    def PrintStats(self):
        logger.info('Direction flag: ' + self.flag_direction)
        logger.info('Confidence interval: [' + str(self.smallest_normal) + ';' + str(self.biggest_normal) + ']')
        logger.info('Mean: ' + str(self.mean))
        logger.info('Mediana: ' + str(self.median))
        logger.info('R: ' + str(self.R))
        logger.info('Epsilon: ' + str(self.epsilon))
        if self.calc_hist:
            logger.info('Mu: ' + str(self.mu))
            logger.info('Sigma: ' + str(self.sigma))
            logger.info('a: ' + str(self.a))
            logger.info('b: ' + str(self.b))
            logger.info('p: ' + str(self.p))
        logger.info('------------')

    # Split normal and abnormal fragments
    def SplitNormAbnorm(self):
        logger.info('Splitting normal and abnormal fragments...')

        self.ff_abn = []
        self.fr_abn = []
        self.rf_abn = []
        self.rr_abn = []        
        self.normal_fragments = []

        self.rr_abn = [frag for frag in self.rr if frag.unique_flag]
        self.ff_abn = [frag for frag in self.ff if frag.unique_flag]
        if self.flag_direction == 'fr':
            self.rf_abn = [frag for frag in self.rf if frag.unique_flag]
            self.fr_abn = [frag for frag in self.fr if frag.unique_flag and frag.is_abnormal(self.smallest_normal, self.biggest_normal)]
            self.normal_fragments = [frag for frag in self.fr if not frag.is_abnormal(self.smallest_normal, self.biggest_normal)]
        else:
            self.fr_abn = [frag for frag in self.fr if frag.unique_flag and frag.length >= 1000]
            self.rf_abn = [frag for frag in self.rf if frag.unique_flag and frag.is_abnormal(self.smallest_normal, self.biggest_normal)]
            self.normal_fragments = [frag for frag in self.rf if not frag.is_abnormal(self.smallest_normal, self.biggest_normal)]

        logger.info('Done.')
        self.num_abnormal = len(self.fr_abn) + len(self.rf_abn) + len(self.rr_abn) + len(self.ff_abn)
        logger.info('Abnormal fragments: ' + str(self.num_abnormal))
        logger.info('Normal fragments: ' + str(len(self.normal_fragments)))
        logger.info('Lengths of abnormal arrays: ' + str(len(self.fr_abn)) + ' ' + 
                                                     str(len(self.ff_abn)) + ' ' + 
                                                     str(len(self.rr_abn)) + ' ' + 
                                                     str(len(self.rf_abn)))

    # Write normal fragments to a file
    def WriteNormal(self):
        logger.info('Writing normal fragments to file...')
        fr_norm = open(self.config['working_dir'] + self.config['normal_fragments_dir'] + 'normal_fragments_' + self.chr + '.txt', 'w')
        fr_norm.write( '| chr | begin | unique_flag | mapp_qul_flag | name | \n')
        for frag in self.normal_fragments:
            fr_norm.write(frag.first_read_chr + ' ' + \
                            str(frag.begin) + ' ' + \
                            str(frag.unique_flag) + ' ' + \
                            str(frag.mapp_qul_flag) + ' ' + 
                            frag.name + '\n')
        fr_norm.close()

    # Cluster fragments with different direction types
    def PerformClustering(self):
        utils.clust(self.fr_abn, self.chr, 'fr', self.biggest_normal, self.smallest_normal, self.config, self.flag_direction)
        utils.clust(self.rf_abn, self.chr, 'rf', self.biggest_normal, self.smallest_normal, self.config, self.flag_direction)
        utils.clust(self.ff_abn, self.chr, 'ff', self.biggest_normal, self.smallest_normal, self.config, self.flag_direction)
        utils.clust(self.rr_abn, self.chr, 'rr', self.biggest_normal, self.smallest_normal, self.config, self.flag_direction)
    
    def SerializeStatsToTmp(self):
        stats_to_serialize = dict()
        stats_to_serialize['num_all_abn'] = self.num_abnormal
        stats_to_serialize['R'] = self.R
        stats_to_serialize['smallest_normal'] = self.smallest_normal
        stats_to_serialize['biggest_normal'] = self.biggest_normal
        stats_to_serialize['median'] = self.median
        stats_to_serialize['flag_direction'] = self.flag_direction
        if self.calc_hist:
            stats_to_serialize['a'] = float(self.a)
            stats_to_serialize['b'] = float(self.b)
            stats_to_serialize['p'] = float(self.p)

        stats_temp_fname = '/tmp/stats_' + self.chr + '.tmp'
        stats_temp_file = open(stats_temp_fname, 'w')
        stats_temp_file.write(yaml.dump(stats_to_serialize))
        stats_temp_file.close()


    # Load .sam file into array of fragments (pysam-based)
    def __LoadSam(self, sam_file):
        num_unmapped_pairs = 0
        num_bad_quality = 0
        reads_processed = 0
        repeats = 0
        k=0
        logger.info('Reading ' + sam_file + ', creating a Read class from every line...')
        (name, extension) = os.path.splitext(sam_file)
        if extension == '.sam':
            sam_in = pysam.Samfile(sam_file, 'r')
        elif extension == '.gz': # .sam.gz
            sam_in = pysam.Samfile(sam_file, 'r')
        elif extension == '.bam':
            sam_in = pysam.Samfile(sam_file, 'rb')
        else:
            logger.info('Unsupported extension in load_sam_file')
            return
        
        mate_reads = dict()
        short_keys = set() 
        for read in sam_in.fetch():
            reads_processed += 1
            if reads_processed % 1000000 == 0:
                logger.info('Reads processed: ' + str(reads_processed))

            if read.is_unmapped or read.mate_is_unmapped:
                num_unmapped_pairs += 1
                continue

            # We use tuple (read 1 tid, read 1 pos, read 2 tid, read 2 pos) as a key to combine 2 reads
            # If a read with this key not observed yet, save it to a dictionaty
            # Otherwise combine 2 reads and remove key from dictionary
            # Short keys are keys without name - use them not to 
            if read.is_read1:
                key = (read.qname,read.tid, read.pos, read.rnext, read.pnext)
                short_key = (read.tid, read.pos, read.rnext, read.pnext)
            else:
                key = (read.qname,read.rnext, read.pnext, read.tid, read.pos)
                short_key = (read.rnext, read.pnext, read.tid, read.pos)
            
            
            if key not in mate_reads:
                mate_reads[key] = read
            else:
                # Protection against wrong duplicates processing
                if short_key in short_keys:
                    continue
                else:
                    short_keys.add(short_key)

                frag = fragment.Fragment()
                frag.from_reads(read, mate_reads[key], sam_in)
                
                if frag.unique_flag or frag.mapp_qul_flag:
        
                    self.__AppendFrag(frag)
                else:
                    continue

                del mate_reads[key]

        sam_in.close()
        logger.info('Finished loading ' + sam_file + ', sorting...')
        self.ff.sort(key = lambda frag: frag.middle)
        self.rr.sort(key = lambda frag: frag.middle)
        self.fr.sort(key = lambda frag: frag.middle)
        self.rf.sort(key = lambda frag: frag.middle)
        logger.info('Done.')
        # Some stats for loaded data
        frags_totally = len(self.rr) + len(self.fr) + len(self.rf) + len(self.ff) + len(self.translocations)
        logger.info('In mate_reads dictionary: ' + str(len(mate_reads)))
        logger.info('Unmapped pairs: ' + str(num_unmapped_pairs))
        logger.info('Repeats : ' + str(repeats))
        logger.info('Bad quality fragments : ' + str(num_bad_quality))
        logger.info('Translocations: ' + str(len(self.translocations)))
        logger.info('Pair ended (fr) fragments: ' + str(len(self.fr)))
        logger.info('Mate ended (rf) fragments: ' + str(len(self.rf)))
        logger.info('Double forward fragments: ' + str(len(self.ff)))
        logger.info('Double reverse fragments: ' + str(len(self.rr)))
        logger.info('Fragments totally: ' + str(frags_totally))
        logger.info('Memory consumed: ' + str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024) + \
            ' (it is Kb for OSX and Mb for Linux)')
        logger.info('------------')

    # Append frag to translocations list
    # Or to one of rr, fr, rf, ff lists
    def __AppendFrag(self, frag):
        if frag.first_read_chr != frag.second_read_chr:
            self.translocations.append(frag)
        elif frag.direction == 'rr':
            self.rr.append(frag)
        elif frag.direction == 'fr':
            self.fr.append(frag)
        elif frag.direction == 'rf':
            self.rf.append(frag)
        elif frag.direction == 'ff':
            self.ff.append(frag)
        else:
            logger.error('Can\'t determine list to append fragment')
