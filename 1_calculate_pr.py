#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2018 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run script for processing of YOUR PROJECT HERE




#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
#whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/bin/pipeline/'

#sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import argparse
import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
import subprocess
import pickle
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================







#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#RNA-Seq
data_file = '%sdata_tables/171207_SUM159_FKBP_DHX15_CRISPR_D5_FASTQ_TABLE.txt' % (projectFolder)

#gene gtf
gene_gtf_path = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'


#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

class Genome:
    __chrDict = dict()
    __featureDict = dict()

    #at its core, the genome has attributes of a build name, a fasta directory and an annotation file
    def __init__(self,name,genome_directory,annot_file):
        self._name = name
        self._directory = genome_directory
        self._annot = annot_file

    def name(self):
        return self._name
    
    def directory(self):
        return self._directory

    def annot(self):
        return self._annot

    def addFeature(self,feature,path):

        if feature in self.__featureDict:
            print('WARNING OVERRIDING %s PATH WITH %s' % (feature,path))
        self.__featureDict[feature] = path

    def returnFeature(self,feature):

        #tries to load the selected feature from the feature dictionary

        if feature not in self.__featureDict:
            print('ERROR: GENOME %s DOES NOT HAVE FEATURE %s' % (self.name(),feature))
            sys.exit()
        else:
            return self.__featureDict[feature]

    def hasFeature(self,feature):

        if feature in self.__featureDict:
            return True
        else:
            return False



#================================================================================
#=================================FUNCTIONS======================================
#================================================================================


def loadGenome(genome_build,config_file = '',verbose=False):

    '''
    loads annotation for a genome into a genome object
    '''

    #this nested dictionary has all of the useful information and likely will have to be
    #edited so it can be configured any time
    genome_build = string.upper(genome_build)

        
    genomeDict = {
        'HG38':{'annot_file':'%sannotation/hg38_refseq.ucsc' % (pipeline_dir),
                'genome_directory':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/',
                'gene_gtf':'/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
                }

        }

    #allow an optional config file to overwrite default paths
    if len(config_file) >0:
        config_table = utils.parseTable(config_file,'\t')
        for line in config_table[1:]:
            (build,field,feature_path) = line[0].split(':')
            genomeDict[string.upper(build)][string.lower(field)] = feature_path
    
    if genome_build not in genomeDict:
        print('ERROR: UNSUPPORTED GENOME BUILD %s. EXITING NOW' % (genome_build))
        sys.exit()
    else:
        if verbose:
            print('USING BUILD %s WITH FOLLOWING FIELDS:' % (genome_build))
            print(genomeDict[genome_build])

    #now attempt to load the genome
    genome = Genome(genome_build,genomeDict[genome_build]['genome_directory'],genomeDict[genome_build]['annot_file'])

    #adding additional optional features
    genome.addFeature('gene_gtf',genomeDict[genome_build]['gene_gtf'])

    return genome



#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():



    parser = argparse.ArgumentParser(usage='%(prog)s [additional options]  -g [GENOME] -d [DATA_FILE] -o [OUTPUTFOLDER]')

    # required flags
    parser.add_argument("-d", "--data_file", dest="data_file",type=str,
                        help="Enter path to a standard lin lab format data table", required=True)
    parser.add_argument("-g", "--genome", dest="genome", type=str,
                        help="specify a genome, HG38 is currently supported", required=True)
    # output flag
    parser.add_argument("-o", "--output", dest="output", type=str,
                        help="Enter the output folder.", required=True)

    # additional options
    parser.add_argument("-title", "--analysis_title", dest="analysis_title", default='',
                        help="specify a title for the analysis")
    parser.add_argument("-names", "--names", dest="names_list", default='',type=str,
                        help="Provide a comma separated list of dataset names to analyze. Default analyzes all datasets in the data table")
    parser.add_argument("-genes", "--gene_list", dest="gene_list_path", default='',type=str,
                        help="To analyze only a subest of genes. Provide a text file where the first column contains refseq NM or NR IDs ")
    parser.add_argument("-gtf", "--gene-gtf", dest="gene_gtf", default='',type=str,
                        help="Provide path of a gene gtf or a gene_gtf.pkl from a previous run. Default will generate a genome-wide gene dictionary pkl based off the genome choice")


    args = parser.parse_args()

    #print(args)

    print('\n\n')
    print('#======================================================================')
    print('#======================I. CHECKING INPUT PARAMETERS====================')
    print('#======================================================================')
    print('\n\n')

        
    #setting up the analysis_title for the project
    if args.analysis_title == '':
        analysis_name = args.data_file.split('/')[-1].split('.')[0]
    else:
        analysis_name = args.analysis_title
    print('Running processing ratio analysis for project %s\n' % (analysis_name))

    #getting genome
    genome_build = args.genome
    genome = loadGenome(genome_build,'',False)
    
    #getting gene gtf
    if args.gene_gtf == '':
        gene_gtf_path = genome.returnFeature('gene_gtf')
    else:
        gene_gtf_path = args.gene_gtf

    print('Using genome build %s and gene gtf %s\n' % (genome_build, gene_gtf_path))

    gene_list_path = args.gene_list_path
    if gene_list_path != '':

        gene_table= utils.parseTable(gene_list_path,'\t')
        ref_id_list = [line[0] for line in gene_table if line[0][0:2] == 'NM' or line[0][0:2] == 'NR']
        print('Using %s to get a gene list of %s genes for analylsis' % (gene_list_path,len(ref_id_list)))
    else:
        ref_id_list = []
        
    #setting up the output folder
    output_folder = utils.formatFolder(args.output,True)
    print('Writing output to %s\n' % (output_folder))
    os.chdir(projectFolder)

    #verifying the data table
    data_file = args.data_file
    pipeline_dfci.summary(data_file)
    data_dict=  pipeline_dfci.loadDataTable(data_file)

    #checking the names list for datasets to analyze
    if args.names_list == '':
        names_list = data_dict.keys()
    else:
        names_list = args.names_list.split(',')
    names_list.sort()
    print('Running processing ratio analysis on the following datasets:')
    print(names_list)

    print('\n\n')
    print('#======================================================================')
    print('#======================II. LOADING GENE GTF============================')
    print('#======================================================================')
    print('\n\n')

    
    #figure out if the gene gtf is a pickle or a text file
    if gene_gtf_path.split('.')[-1] == 'pkl':
        gene_dict_path = gene_gtf_path
        print('Loading gene gtf from an existing .pkl at %s' %(gene_dict_path))
        gene_dict = loadGeneDict(gene_dict_path)
    elif gene_gtf_path.split('.')[-1] == 'gtf':
        gene_dict_path = '%s%s_%s_gene_gtf.pkl' % (output_folder,genome_build,analysis_name)
        print('Making a new gene dict from %s and writing to %s' % (gene_gtf_path,gene_dict_path))
        gene_dict = makeGeneDict(gene_gtf_path,gene_dict_path,ref_id_list)
    else:
        print('Error: gene gtf file needs to end with either .pkl or .gtf')


    print('\n\n')
    print('#======================================================================')
    print('#================IV. PULLING READS FROM TRANSCRIPTS====================')
    print('#======================================================================')
    print('\n\n')


    calculatePR(genome,data_file,gene_dict,ref_id_list,names_list,analysis_name,output_folder)
    

#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~MAKING GENE DICT FROM A GTF~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def makeGeneDict(gene_gtf_path,output='',ref_list = []):

    '''
    makes a dictionary of genes keyed by NM ID from a gene gtf
    '''
    exon_dict =defaultdict(list)
    cds_dict = defaultdict(list)
    common_dict = defaultdict(list)
    sense_dict = defaultdict(list)
    chr_dict = defaultdict(list)
    #gene_gtf = utils.parseTable(gene_gtf_path,'\t')
    gene_gtf = open(gene_gtf_path,'r')
    ticker = 0
    
    if ref_list == []:
        filter_genes = False
    else:
        filter_genes = True

    print('parsing gene gtf file: %s' % (gene_gtf_path))
    print('Number of lines:')
    for line in gene_gtf:
        line = line.rstrip().split('\t')
        if ticker % 100000 == 0:
            print(ticker)
        if line[2] == 'exon' or line[2] == 'CDS':
            transcript_id = re.findall('transcript_id "N[^;]*',line[-1])[0].split('"')[1]
            if filter_genes and ref_list.count(transcript_id) > 0:

                gene_id = re.findall('gene_id "[^;]*',line[-1])[0].split('"')[1]
                sense_dict[transcript_id].append(line[6])
                chr_dict[transcript_id].append(line[0])
                common_dict[transcript_id].append(gene_id)
                if line[2] == 'exon':
                    exon_dict[transcript_id].append((line[3],line[4]))
                if line[2] == 'CDS':
                    cds_dict[transcript_id].append((line[3],line[4]))

        ticker+=1


    #now try to assemble a gene dict
    transcript_list = chr_dict.keys()
    print('making gene dict from %s genes' % (len(transcript_list)))
    print('number of genes:')
    gene_dict = {}
    ticker =0 
    for transcript_id in transcript_list:

        if ticker % 10000 == 0:
            print(ticker)
        ticker+=1
        chrom = utils.uniquify(chr_dict[transcript_id])[0]
        sense = utils.uniquify(sense_dict[transcript_id])[0]
        commonName = utils.uniquify(common_dict[transcript_id])[0]
        exStarts = [int(exon[0]) for exon in exon_dict[transcript_id]]
        exEnds = [int(exon[1]) for exon in exon_dict[transcript_id]]

        tx_coords = exStarts + exEnds
        txStart = min(tx_coords)
        txEnd = max(tx_coords)
        if len(cds_dict[transcript_id]) > 0: # for coding genes
            cdsStarts = [int(cds_exon[0]) for cds_exon in cds_dict[transcript_id]] 
            cdsEnds = [int(cds_exon[1]) for cds_exon in cds_dict[transcript_id]]
            cds_coords = cdsStarts +cdsEnds
            cdsStart = min(cds_coords)
            cdsEnd = max(cds_coords)
            gene = utils.Gene(transcript_id,chrom,sense,[txStart,txEnd],[cdsStart,cdsEnd],exStarts,exEnds,commonName)
        else: # for non coding genes
            gene = utils.Gene(transcript_id,chrom,sense,[txStart,txEnd],[0,0],exStarts,exEnds,commonName)
        gene_dict[transcript_id] = gene
    
    # print(gene_dict.keys())
    # gene = gene_dict['NR_024321']
    # print(gene.commonName())
    # print(gene.chr())
    # print([exon.coords() for exon in gene.cdExons()])
    # print([exon.coords() for exon in gene.txExons()])
    # print([intron.coords() for intron in gene.introns()])
    
    if len(output) >0: #i can pickle this

        output_file = open(output,'w')
        pickle.dump(gene_dict,output_file)
        output_file.close()
        print('wrote gene dict as a pickle to %s' % (output_file))
        return gene_dict
 

def makeSubDict(gene_dict_path,output_path,ref_list):
    '''
    takes a gene dict and makes a subsetted dict
    '''

    #first check if the output already exists
    if utils.checkOutput(output_path,0,0):
        dict_file = open(output_path,'rb')
        sub_dict = pickle.load(dict_file)
        print('FOUND EXISTING GENE DICTIONARY AT %s' % (output_path))
        return sub_dict

    gene_dict = loadGeneDict(gene_dict_path)

    sub_dict = {}
    for ref_ID in ref_list:
        sub_dict[ref_ID] = gene_dict[ref_ID]


    output_file = open(output_path,'w')
    pickle.dump(sub_dict,output_file)
    output_file.close()

    return sub_dict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~I CAN PICKLE THAT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def cached(cachefile):
    """
    A function that creates a decorator which will use "cachefile" for caching the results of the decorated function "fn".
    taken from https://datascience.blog.wzb.eu/2016/08/12/a-tip-for-the-impatient-simple-caching-with-python-pickle-and-decorators/
    """
    def decorator(fn):  # define a decorator for a function "fn"
        def wrapped(*args, **kwargs):   # define a wrapper that will finally call "fn" with all arguments            
            # if cache exists -> load it and return its content
            if os.path.exists(cachefile):
                    with open(cachefile, 'rb') as cachehandle:
                        print("using cached result from '%s'" % cachefile)
                        return pickle.load(cachehandle)

            # execute the function with all arguments passed
            res = fn(*args, **kwargs)

            # write to cache file
            with open(cachefile, 'wb') as cachehandle:
                print("saving result to cache '%s'" % cachefile)
                pickle.dump(res, cachehandle)

            return res

        return wrapped

    return decorator   # return this "customized" decorator that uses "cachefile"

#@cached('slow_function_cache.pickle')
def loadGeneDict(gene_dict_path):
    '''
    loads the pickle
    '''
    print('loading gene dict from %s' % (gene_dict_path))
    input_file = open(gene_dict_path,'r')
    gene_dict = pickle.load(input_file)
    return gene_dict



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~MAKING EXP DICTIONARY~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def makeExpTable(exp_file_list,output=''):

    '''
    makes an expression table w/ the average exp in each timepoint/condition
    '''
    
    #first set the table up
    output_table = [['GENE']]
    
    #do the first one manually to grab gene names
    exp_table = utils.parseTable(exp_file_list[0],'\t')

    output_table[0]+=exp_table[0]
    for line in exp_table[1:]:
        output_table.append(line)

    if len(exp_file_list) > 1:
        for exp_file in exp_file_list[1:]:
            exp_table = utils.parseTable(exp_file,'\t')
            output_table[0]+= exp_table[0]
            for i in range(1,len(exp_table)):
                output_table[i]+=exp_table[i][1:]
        
    if len(output) >0:
        utils.unParseTable(output_table,output,'\t')
        return output

    else:
        return output_table

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~GETTING NON OVERLAPPING GENES~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        

    
def getNonOverlappingGenes(exp_table_path,gene_dict,cutOff =10.0,output=''):
    '''
    gets all expressed genes that do not overlap a different detectable gene
    '''

    #if already exists returns output
    if len(output) > 0:
        if utils.checkOutput(output,0,0):
            print('identified previous output at %s' % (output))
            filtered_ref_table = utils.parseTable(output,'\t')
            return filtered_ref_table

    cutOff = float(cutOff)

    #first pull up the expression table and get a list of expressed genes
    #by common name
    exp_table = utils.parseTable(exp_table_path,'\t')
    exp_genes = [line[0] for line in exp_table[1:] if max([float(x) for x in line[1:]]) > cutOff]
    exp_genes = utils.uniquify(exp_genes)

    #also get detectable genes > 1 fpkm
    detectable_genes = [line[0] for line in exp_table[1:] if max([float(x) for x in line[1:]]) > 1]
    detectable_genes = utils.uniquify(detectable_genes)

    print('Identified %s expressed genes with expression > %s in at least 1 sample in %s' % (len(exp_genes),cutOff,exp_table_path))

    #now we want to get a refseq list of all expressed genes
    all_ref_id_list = gene_dict.keys()
    name_to_ref_dict = defaultdict(list)
    expressed_ref_list = [] #fpkm > cutOff
    detectable_ref_list = [] #fpkm > 1

    for ref_id in all_ref_id_list:
        common_name = gene_dict[ref_id].commonName()
        if exp_genes.count(common_name) >0:
            name_to_ref_dict[common_name].append(ref_id)
            expressed_ref_list.append(ref_id)
        if detectable_genes.count(common_name) >0:

            detectable_ref_list.append(ref_id)


    print('Identified expressed genes by refseq ID')
    expressed_ref_list = utils.uniquify(expressed_ref_list)
    print(len(expressed_ref_list))


    #now for each gene, we only want those that do not overlap another gene (period)
    #tx_loci = [gene_dict[refID].txLocus() for refID in all_ref_id_list]
    tx_loci = []
    for ref_id in detectable_ref_list:
        gene = gene_dict[ref_id]
        tx_locus = utils.Locus(gene.chr(),gene.txLocus().start(),gene.txLocus().end(),gene.sense(),gene.commonName())
        tx_loci.append(tx_locus)

    print('making transcript collection')
    tx_collection = utils.LocusCollection(tx_loci)

    filtered_ref_table = [['REF_ID','GENE']]
    for ref_id in expressed_ref_list:

        gene_locus = gene_dict[ref_id].txLocus()
        overlapping_genes = tx_collection.getOverlap(gene_locus,'both')
        overlapping_names = [locus.ID() for locus in overlapping_genes]

        if len(utils.uniquify(overlapping_names)) == 1:
            filtered_ref_table.append([ref_id,gene_dict[ref_id].commonName()])

    print('after filtering %s genes remaining' % (len(filtered_ref_table)))

    if len(output) >0:
        utils.unParseTable(filtered_ref_table,output,'\t')
    return filtered_ref_table


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~PULLING GENE READS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def decompose_sam_read(read,genomeDirectory,n_segment=10,debug=False):

    '''
    splits up a read
    '''
    import re
    import operator
    operator_dict = {'I': max, #this is a dirty hack to get it to do nothing since the span will always be less than
                     #the read coordinate unless you are magically at the very end of the chrom
                       'M': operator.add,
                       'D': operator.add,
                       'N': operator.add,
                       'S': max,
                       }

    sense_dict = { '99':'+',
                   '163':'+',
                   '83':'-',
                   '147':'-',
                   }

    
    read_chrom = read[2]
    cigar = read[5]
    read_seq = read[9]
    read_start = int(read[3])
    read_sense = sense_dict[read[1]]
    cigar_len_list = [x for x in re.split('[A-Z]',cigar) if x != '']
    cigar_operator_list = [x for x in re.split('[0-9]',cigar) if x != '']
    if debug:
        print(cigar_len_list)
        print(cigar_operator_list)
        print('CIGAR STRING IS: %s' % (cigar))
    
    read_loci = []
    segment_start = int(read_start)

    for i in range(len(cigar_len_list)):
        cig = cigar_operator_list[i]
        span = int(cigar_len_list[i])

        segment_stop = operator_dict[cig](segment_start,span) -1
        # if cig =='N':
        #     segment_stop +=1
        if debug:
            print(segment_start,segment_stop)

        if cig =='M' and span > n_segment:
            if debug:
                genome_seq = utils.fetchSeq(genomeDirectory,read_chrom,segment_start,segment_stop,False)
                print('this is what i think it is')
                print(genome_seq)
                print('this is the original sequence')
                print(read_seq)
                print('adding to the segment loci')
            segment_locus = utils.Locus(read_chrom,segment_start,segment_stop,read_sense)
            read_loci.append(segment_locus)
        segment_start = int(segment_stop) + 1

    #print(utils.fetchSeq(genomeDirectory,read_chrom,67198599,67198672))

    return read_loci



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAKE BED STRING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def makeBedLine(read1_loci,read2_loci,first_read_strand='.',read_id ='.'):

    '''
    red = first strand +
    blue = first strand -
    thick = read part
    thin = fragment
    dashed = gaps in read
    '''
    #get all of the read bits
    read1_lengths = [locus.len() for locus in read1_loci]
    read1_starts = [locus.start() for locus in read1_loci]
    read1_stops = [locus.end() for locus in read1_loci]
    read1_coords=  read1_starts + read1_stops
    read1_coords.sort()

    read2_lengths = [locus.len() for locus in read2_loci]
    read2_starts = [locus.start() for locus in read2_loci]
    read2_stops = [locus.end() for locus in read2_loci]
    read2_coords=  read2_starts + read2_stops
    read2_coords.sort()

    overall_coords = read1_coords + read2_coords
    
    #bed field 1
    chrom =read1_loci[0].chr() # if read is on multiple chroms, everything blows up
    
    #bed fields 2 and 3
    chromStart = min(overall_coords)
    chromEnd = max(overall_coords) + 1
    
    #for now leave name blank to save disk space
    #bed field 4 and 5
    name = read_id
    score = 0
    
    #for strand go by strand of first read
    #bed field 6
    strand = first_read_strand
    
    #for thicc parts, first read should always be thicc
    #bed fields 7,8,9
    if strand == '+' or strand == '.':
        thickStart = read1_coords[0]
        thickEnd = read1_coords[-1]
        itemRgb='255,0,0'
    else:
        thickStart = read2_coords[0]
        thickEnd = read2_coords[-1] + 1
        itemRgb='0,0,255'
    
    #bed field 10
    blockCount = len(read1_loci) + len(read2_loci)

    blockSizes_list = [int(x) for x in read1_lengths + read2_lengths]
    blockSizes = ','.join([str(x) for x in blockSizes_list])

    blockStarts_list = [int(x-chromStart) for x in read1_starts + read2_starts]
    blockStarts = ','.join([str(x) for x in blockStarts_list])

    blockEnds = [chromStart + blockStarts_list[i] + blockSizes_list[i] for i in range(len(blockStarts_list))]
    if blockEnds[-1] == chromEnd and blockStarts_list[0] == 0:

        bed_line = [chrom,chromStart,chromEnd,name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts]

        return bed_line
    else:
        return []
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~GET PROCESSING RATIO~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def get_processing_ratio(bam,gene,genomeDirectory,n_segment=10,return_unprocessed=False):

    '''
    calculates the ratio of unprocsesed to processed reads
    requires read segments > n_segment (default 10)
    '''
    
    intron_collection = utils.LocusCollection(gene.introns())
    exon_collection = utils.LocusCollection(gene.txExons())
    
    all_reads = bam.getRawReads(gene.txLocus(),'.',False,True,False)
    
    #order is important here for the flags because we will use this later
    #to always structure pairs that are +/- orientation on the actual genome
    #e.g. 163 is second in the pair, but a + strand read
    #for fr library types, + strand genes will be 163/83
    #- strand genes will be 99/147
    flag_filter_dict = {'+':['163','83'], #2nd is + strand #1st is - strand
                        '-':['99','147'], #1st is + strand #2nd is - strand
                        }

    #now filter to sense reads for this particular transcript
    filtered_reads = [read for read in all_reads if flag_filter_dict[gene.sense()].count(read[1]) >0]

    #now we need to decompose into pairs
    pair_dict = defaultdict(list)
    for read in filtered_reads:
        pair_dict[read[0]].append(read) #since reads in the bam are ordered along the DNA + strand, our desired pair orientation is preserved

    #now keep only full pairs
    read_id_list = pair_dict.keys()

    for read_id in read_id_list:

        if len(pair_dict[read_id]) != 2:
            pair_dict.pop(read_id) #note, popping the keys does not update the read_id_list

    read_id_list = pair_dict.keys() #updated to remove non pairs
    if len(read_id_list) < 10:
        print('Less than 10 total reads for %s %s' % (gene.name(),gene.commonName()))
        if return_unprocessed:
            return 0,0,0,0,[]
        else:
            return 0,0,0,0
    #create some temp buckets to store these
    processed_reads = []
    unprocessed_reads = []
    cryptic_reads = []
    ticker = 0

    # print(pair_dict['NS500589:294:H72CMBGX5:4:12408:11118:14417'])
    # sys.exit()

    #using the leftmost read to figure out pair orientation
    pair_orientation_dict = {'163':'-',
                             '99':'+',
                             }

    if return_unprocessed:
        unprocessed_bed = []

    for read_id in read_id_list:
        
        #stores for each segment whether it's exonic (E), intronic (I), or straddling a boundary (EI)
        #E and E consecutive segments would denote an EE boundary
        first_read_segments = []
        second_read_segments = []

        #decompose the reads into individual segments
        first_read_loci = decompose_sam_read(pair_dict[read_id][0],genomeDirectory,n_segment)
        second_read_loci = decompose_sam_read(pair_dict[read_id][1],genomeDirectory,n_segment)

        #if we want to spit back the bed, need to figure out direction of first strand
        if return_unprocessed == True:
            first_read_strand = pair_orientation_dict[pair_dict[read_id][0][1]]

        #first read
        #assume processed unless told otherwise

        for read_locus in first_read_loci:
            if intron_collection.getContainers(read_locus): #read segment entirely contained

                first_read_segments.append('I')
            elif intron_collection.getOverlap(read_locus):

                first_read_segments.append('EI')
            else:
                first_read_segments.append('E')

        #second read
        #assume processed unless told otherwise

        for read_locus in second_read_loci:
            if intron_collection.getContainers(read_locus): #read segment entirely contained

                second_read_segments.append('I')
            elif intron_collection.getOverlap(read_locus):

                second_read_segments.append('EI')
            else:
                second_read_segments.append('E')

        if first_read_segments.count('EI') > 0 or second_read_segments.count('EI') > 0:

            
            unprocessed_reads.append(read_id)
            if return_unprocessed:
                bed_line = makeBedLine(first_read_loci,second_read_loci,first_read_strand,'.')
                if bed_line != []:
                    unprocessed_bed.append(bed_line)
            # print(read_id)
            # print(first_read_segments)
            # print([locus.coords() for locus in first_read_loci])
            # print(second_read_segments)
            # print([locus.coords() for locus in second_read_loci])
            # print([locus.coords() for locus in gene.txExons()])
            # print('\n\n\n')
        elif utils.uniquify(first_read_segments + second_read_segments) == ['E']:
            processed_reads.append(read_id)

        else:
            cryptic_reads.append(read_id)
            if return_unprocessed:
                bed_line = makeBedLine(first_read_loci,second_read_loci,first_read_strand,'.')
                if bed_line != []:
                    unprocessed_bed.append(bed_line)


        ticker+=1

    total_reads = len(read_id_list)
    processed_ratio = round(len(processed_reads)/float(len(read_id_list)),3)
    unprocessed_ratio = round(len(unprocessed_reads)/float(len(read_id_list)),3)
    cryptic_ratio = round(len(cryptic_reads)/float(len(read_id_list)),3)
    if return_unprocessed:
        return total_reads,processed_ratio,unprocessed_ratio,cryptic_ratio,unprocessed_bed
    else:
        return total_reads,processed_ratio,unprocessed_ratio,cryptic_ratio
    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~MAKE READ PROCESSING TABLE~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def calculatePR(genome,data_file,gene_dict,ref_list = [],names_list = [],analysis_name = '',output_folder =''):

    '''
    spits out a table for processed, unprocessed, and cryptic for all genes/bams provided
    '''
    
    genome_build = genome.name()
    genomeDirectory = genome.directory()
    #set the output folder
    if len(output_folder) > 0:
        output_folder = utils.formatFolder(output_folder)
    else:
        output_folder = './'
        
    print('writing output to %s' % (output_folder))

    #load the data file
    data_dict = pipeline_dfci.loadDataTable(data_file)

    if len(names_list) == 0:
        names_list = data_dict.keys()

    #get all of the bams involved
    print('loading bams from %s' % (data_file))
    print(names_list)
    bam_list = [utils.Bam(data_dict[name]['bam']) for name in names_list]
    
    #set up the output tables
    total_table = [['REF_ID','GENE'] + names_list]
    processed_table = [['REF_ID','GENE'] + names_list]
    cryptic_table = [['REF_ID','GENE'] + names_list]
    unprocessed_table = [['REF_ID','GENE'] + names_list]

    #now we need to set up the bed files that we're writing to for each bam
    bed_path_list = ['%s%s_%s_unprocessed.bed' % (output_folder,genome_build,name) for name in names_list]
    for i in range(len(bed_path_list)):
        bed_path = bed_path_list[i]
        bed_file = open(bed_path,'w')
        #write the track header
        header = 'track name=%s description="%s unprocesed reads" itemRgb=On\n' % (names_list[i],names_list[i])
        bed_file.write(header)
        bed_file.close()
    
    #get your gene list
    if len(ref_list) == 0:
        ref_list = gene_dict.keys()

    print('iterating through each gene')
    #now iterate through genes as the first loop
    ticker = 0
    for ref_id in ref_list:
        #print(ticker)
        if ticker % 100 == 0:
            print(ticker)
        # if ticker == 10:
        #    break

        gene = gene_dict[ref_id]
        total_vector = [ref_id,gene.commonName()]
        processed_vector = [ref_id,gene.commonName()]
        unprocessed_vector =[ref_id,gene.commonName()]
        cryptic_vector = [ref_id,gene.commonName()]
        for i in range(len(bam_list)):
            bam = bam_list[i]
            
            total_reads,processed_ratio,unprocessed_ratio,cryptic_ratio,unprocessed_bed = get_processing_ratio(bam,gene,genomeDirectory,n_segment=10,return_unprocessed=True)

            if len(unprocessed_bed) > 0:
                bed_file = open(bed_path_list[i],'a')
                for line in unprocessed_bed:
                    bed_file.write('\t'.join([str(x) for x in line]) +'\n')
                bed_file.close()
            
            total_vector.append(total_reads)
            processed_vector.append(processed_ratio)
            unprocessed_vector.append(unprocessed_ratio)
            cryptic_vector.append(cryptic_ratio)
        
        #add this gene line for all of the tables
        total_table.append(total_vector)
        processed_table.append(processed_vector)
        unprocessed_table.append(unprocessed_vector)
        cryptic_table.append(cryptic_vector)
        ticker +=1


    #now write the output
    if len(analysis_name) ==0:
        #grab the name from the data_file
        analysis_name = data_file.split('/')[-1].split('.')[0]

        
    total_output = '%s%s_total.txt' % (output_folder,analysis_name)
    processed_output = '%s%s_processed.txt' % (output_folder,analysis_name)
    unprocessed_output = '%s%s_unprocessed.txt' % (output_folder,analysis_name)
    cryptic_output = '%s%s_cryptic.txt' % (output_folder,analysis_name)

    utils.unParseTable(total_table,total_output,'\t')
    utils.unParseTable(processed_table,processed_output,'\t')
    utils.unParseTable(unprocessed_table,unprocessed_output,'\t')
    utils.unParseTable(cryptic_table,cryptic_output,'\t')

        

    
    
#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()







