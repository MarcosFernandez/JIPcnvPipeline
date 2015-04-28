#!/usr/bin/env python
import os
import json
import argparse
import sys


class CreateConfigurationFile(object):
    """Class which manages Configuration file Manager"""
 
    def __init__(self):
        """Class constructor"""
        #PIPELINE PARAMETERS
        self.run_pcr_duplicates = True            #By default run PCR duplicates process
        self.input_data_se = False                #Input data is in single end mode
        self.json_file = None                     #Name of the json configuration file to be created.


        #BWA PARAMETERS
        self.bwa_reference = None                 #File were ara located BWA index files
        self.picard_path = "/apps/PICARD/1.95/"   #Picard Tools Path Location
        self.bwa_threads = "8"                    #BWA number of threads
        self.tmp_folder = "/scratch_tmp"          #Path to the temporary folder for the machine where the job is running
        self.bwa_time_job = "23:59:59"            #Job Time for the Lustre Job administrator
        
        #MARK DUPLICATES
        self.java_heap = "25g"                    #Memory to reserve for the Java Heap Memory Stack, when running MarkDuplicates
        self.mark_dup_time_job = "23:59:59"       #Job Time for the Lustre Job administrator
        self.mark_dup_threads = "5"               #Number of threads to run mark duplicates

        #BAM TO FASTQ
        self.chunks = "50"                        #Number of chunks to fragments the FASTQ input data to later perform a mapping job per chunk
        self.kmer_length = 36                     #Read length for the chopped reads
        self.windowing = 36                       #while chopping reads you can specify a bp window sliding,
                                                  #if you are not interested in any overlap leave this value as the kmer-length.
        self.first_position = 10                  #Position to start chopping your reads, if you are detecting bad patterns
                                                  #at the begining of the reads you can remove such bases. As a filtering criteria.
        self.bam_to_fastq_time_job = "15:59:59"   #BAM to FASTQ Job time for the lustre Job administrator
        self.bam_to_fastq_threads = "8"           #BAM to FASTQ threads

        #FRAGMENTATION AND CHOPPING
        self.frag_and_chop_time_job = "09:59:59"   #Fragmentation and chopping time job
        self.fragmentation_threads = "8"           #Fragmentation threads

        #CNV MAPPING
        self.gem_index = None                     #Path to the .gem genome index
        self.quality = "ignore"                   #FASTQ qualities encoding.
        self.mismatches = 2                       #Set the allowed mismatch ratio as 0 < mm < 1. 
        self.max_edit_distance = 4                #Maximum edit distance (ratio) allowed for an alignment 
        self.strata_after_best = 2                #The number of strata examined after the best one
        self.max_decoded_matches = 20             #Maximum decoded matches
        self.min_decoded_strata = 1               #Minimum decoded strata
        self.gem_threads = "8"                    #Number of threads to perform gem mapping
        self.gem_time_job = "05:59:59"            #Job time limit for the lustre Job administrator

        #CNV CALLING
        self.reference_conf = None                #Reference configuration file from mrcanavar PREP step
        self.cnv_call_time_job = "09:59:59"       #CNV Calling Job time limit for the lustre Job administrator
        self.cnv_threads = "1"                    #CNV Threads

        #PARAMETERS DICTIONARIES
        self.all_parameters = {}
        self.run_remove_duplicates_dict = {}
        self.bwa_mapping_dict = {}
        self.chunks_chop_dict = {}
        self.cnv_mapping_dict = {}
        self.cnv_calling_dict = {}

    def register_parameter(self, parser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_basic_mapping(parser)
        self.register_rm_duplicates(parser)
        self.register_fragmentation(parser)
        self.register_cnv_mapping(parser)
        self.register_cnv_calling(parser)

    def register_general(self, parser):
        """Register all general parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        general_group = parser.add_argument_group('Pipeline Parameters')
        general_group.add_argument('--no-duplicates', dest="run_pcr_duplicates", action="store_false", help='If specified, do not run remove duplicates steps.')
        general_group.add_argument('-se', dest="input_data_se", action="store_true", help='Input FASTQ data in single end mode.')
        general_group.add_argument('--json-file', dest="json_file", metavar="JSON", help='JSON configuration file to be created')


    def register_basic_mapping(self, parser):
        """Register all basic mapping parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        basic_mapping_group = parser.add_argument_group('Basic Mapping Parameters')
        basic_mapping_group.add_argument('-bwa-reference', dest="bwa_reference", metavar="fastaRef", help='Fasta reference file. In the same folder must be found BWA MEM index files')
        basic_mapping_group.add_argument('-picard-path', dest="picard_path", metavar="PICARD_PATH", help='Picard Tools Path Location. Default %s' % self.picard_path)
        basic_mapping_group.add_argument('-bwa-threads', dest="bwa_threads", metavar="BWA_THREADS", help='Number of threads for BWA. Default %s' % self.bwa_threads)
        basic_mapping_group.add_argument('-tmp-folder', dest="tmp_folder", metavar="TMPDIR",\
                                          help='Path to the temporary folder for the machine where the job is running. Default %s' % self.tmp_folder)
        basic_mapping_group.add_argument('-bwa-time-job', dest="bwa_time_job", metavar="TIME_JOB",\
                                          help='Job Time for the Lustre Job administrator. Default %s' % self.bwa_time_job)


    def register_rm_duplicates(self, parser):
        """Register all remove duplicates parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        rm_duplicates_group = parser.add_argument_group('Remove Duplicates Parameters')
        rm_duplicates_group.add_argument('-java-heap', dest="java_heap", metavar="XXg", help='Memory to reserve for the Java Heap Memory Stack, when running MarkDuplicates.\
                                          Default %s' % self.java_heap)
        rm_duplicates_group.add_argument('-mark-dup-time-job', dest="mark_dup_time_job", metavar="TIME_JOB", help='Job Time for the Lustre Job administrator.\
                                          Default %s' % self.mark_dup_time_job)
        rm_duplicates_group.add_argument('-mark-dup-threads', dest="mark_dup_threads", metavar="THREADS", help='Number of threads for PicardTools. Default %s' % self.mark_dup_threads)

    def register_fragmentation(self,parser):
        """Register all fragmentation parameters with the given
        argparse parser
        parser -- the argparse parser
        """
        fragmentation_group = parser.add_argument_group('Fragmentation and choppping')
        fragmentation_group.add_argument('-chunks', dest="chunks", metavar="CHUNKS", help='Number of chunks to fragment FASTQ input data to later perform a mapping job per chunk.\
                                          The more chunks, the more parallelization. Default %s' % self.chunks)
        fragmentation_group.add_argument('-kmer-length', dest="kmer_length", metavar="KMER_LEN", help='Read length for the chopped reads. Default %i' % self.kmer_length)
        fragmentation_group.add_argument('-windowing', dest="windowing", metavar="OVERLAP", help='While chopping reads you can specify a bp window sliding,\
                                          if you are not interested in any overlap leave this value as the kmer-length. Default %i' % self.windowing) 
        fragmentation_group.add_argument('-first-position', dest="first_position", metavar="POSITION", help='Position to start chopping your reads, \
                                          if you are detecting bad patterns at the begining of the reads you can remove such bases.\
                                          As a filtering criteria. Default %i' % self.first_position)
        fragmentation_group.add_argument('-bam_fastq_time_job', dest="bam_to_fastq_time_job", metavar="TIME", help='Time for a job where reads comes from a bam free of duplicates\
                                          then a tranformation to Fastq will be perform. Default %s' % self.bam_to_fastq_time_job)
        fragmentation_group.add_argument('-bam_fastq_threads', dest="bam_to_fastq_threads", metavar="THREADS", help='Threads for a job where reads comes from a bam free of duplicates\
                                          then a tranformation to Fastq will be perform. Default %s' % self.bam_to_fastq_threads)
        fragmentation_group.add_argument('-frag-chop-time-job', dest="frag_and_chop_time_job", metavar="TIME", help='Fragmentation time job, in case no PCR duplicates removal is perform.\
                                          Default %s' % self.frag_and_chop_time_job)
        fragmentation_group.add_argument('-fragmentation-threads', dest="fragmentation_threads", metavar="THREADS", help='Number of threads for Fragmentation.\
                                          Default %s' % self.fragmentation_threads)
       
    def register_cnv_mapping(self,parser):
        """Register all cnv parameters with the given
        argparse parser
        parser -- the argparse parser
        """  
        cnv_mapping_group = parser.add_argument_group('CNV Mapping')
        cnv_mapping_group.add_argument('-gem-index', dest="gem_index", metavar=".gem", help='Path to the .gem genome index. Default %s' % self.gem_index)
        cnv_mapping_group.add_argument('-mismatches',dest="mismatches", metavar="BASES", help='Maximum number of nucleotide substitutions allowed \
                                        while mapping each read. Default %i' % self.mismatches)
        cnv_mapping_group.add_argument('-max-edit-distance', dest="max_edit_distance", metavar="EDIT_DISTANCE", help='Maximum edit distance (ratio) allowed for an alignment.\
                                       Default %i' % self.max_edit_distance)
        cnv_mapping_group.add_argument('-strata-after-best', dest="strata_after_best", metavar="STRATAS", help='Number of strata examined after the best one.\
                                       Default %i' % self.strata_after_best)
        cnv_mapping_group.add_argument('-max-decoded-matches', dest="max_decoded_matches", help='Maximum decoded matches.\
                                       Default %i' % self.max_decoded_matches)
        cnv_mapping_group.add_argument('-min_decoded_strata', dest="min_decoded_strata", help='Minimum decoded strata.\
                                       Default %i' % self.min_decoded_strata)
        cnv_mapping_group.add_argument('-cnv-mapping-threads', dest="cnv_mapping_threads", metavar="THREADS", help='Number of threads for cnv mapping.\
                                          Default %s' % self.gem_threads)
        cnv_mapping_group.add_argument('-cnv-time-job', dest="cnv_time_job", metavar="TIME", help='Job time limit for cnv mapping in the lustre queue.\
                                          Default %s' % self.gem_time_job)       
      
    def register_cnv_calling(self,parser):
        """Register all cnv calling with the given
        argparse parser
        parser -- the argparse parser
        """  
        cnv_calling_group = parser.add_argument_group('CNV Calling')
        cnv_calling_group.add_argument('-reference-conf', dest="reference_conf", metavar="CONF_FILE", help='Reference configuration file from mrcanavar PREP step.')
        cnv_calling_group.add_argument('-cnv-call-time-job', dest="cnv_call_time_job", metavar="TIME", help='Job time limit for cnv calling in the lustre queue.\
                                        Default %s' % self.cnv_call_time_job)  
        cnv_calling_group.add_argument('-cnv-threads', dest="cnv_threads", metavar="THREADS", help='Number of threads for cnv calling.\
                                          Default %s' % self.cnv_threads)


    def check_parameters(self,args):
        """Check parameters consistency
        args -- set of parsed arguments
        return True if everithing is ok, otherwise false"""
        
        if not args.run_pcr_duplicates:
            self.run_pcr_duplicates = args.run_pcr_duplicates       

        if args.json_file == None:
            print "Sorry! JSON file was not specified!!"
            return False
        else:
            self.json_file = args.json_file

        if args.input_data_se:
            self.input_data_se = args.input_data_se

        if self.run_pcr_duplicates:
            #1. CHECK PARAMETER FOR BASIC MAPPING
            if args.bwa_reference == None:
                print "Sorry! No bwa reference was specified!!"
                return False
            else:
                if not os.path.exists(args.bwa_reference):
                    print args.bwa_reference + " not found"
                    return False
                else:
                    self.bwa_reference = os.path.abspath(args.bwa_reference)

            if args.picard_path:
                if not os.path.exists(args.picard_path):
                    print args.picard_path + " not found"
                    return False
                else:
                    self.picard_path = os.path.abspath(args.picard_path)
         
            if args.bwa_threads:
                self.bwa_threads = args.bwa_threads
            
            if args.tmp_folder:
                self.tmp_folder = args.tmp_folder
            
            if args.bwa_time_job:
                self.bwa_time_job = args.bwa_time_job        
  
            #2. REMOVE DUPLICATES
            if  args.java_heap:
                self.java_heap = args.java_heap
            
            if args.mark_dup_time_job:
                self.mark_dup_time_job = args.mark_dup_time_job
            
            if args.mark_dup_threads:
                self.mark_dup_threads = args.mark_dup_threads 

        #3. FRAGMENTATION
        if args.chunks:
            self.chunks = args.chunks
        if args.kmer_length:
            self.kmer_length = args.kmer_length
        if args.windowing:
            self.windowing = args.windowing
        if args.first_position:
            self.first_position = args.first_position
        if args.bam_to_fastq_time_job:
            self.bam_to_fastq_time_job = args.bam_to_fastq_time_job
        if args.bam_to_fastq_threads:
            self.bam_to_fastq_threads = args.bam_to_fastq_threads
        if args.frag_and_chop_time_job:
            self.frag_and_chop_time_job = args.frag_and_chop_time_job
        if args.fragmentation_threads:
            self.fragmentation_threads = args.fragmentation_threads

        #4. CNV MAPPING
        if args.gem_index == None:
            print "Sorry! No gem reference was specified!!"
            return False
        else:
            if not os.path.exists(args.gem_index):
                print args.gem_index + " not found"
                return False
            else:
                self.gem_index = os.path.abspath(args.gem_index)

        if args.mismatches:
            self.mismatches = args.mismatches
        if args.max_edit_distance:
            self.max_edit_distance = args.max_edit_distance
        if args.strata_after_best:
            self.strata_after_best = args.strata_after_best
        if args.max_decoded_matches:
            self.max_decoded_matches = args.max_decoded_matches
        if args.min_decoded_strata:
            self.min_decoded_strata = args.min_decoded_strata
        if args.cnv_mapping_threads:
            self.gem_threads = args.cnv_mapping_threads
        if args.cnv_time_job:
            self.gem_time_job = args.cnv_time_job

        #5. CNV CALLING
        if args.reference_conf == None:
            print "Sorry! No canavar configuration reference was specified!!"
            return False
        else:
            if not os.path.exists(args.reference_conf):
                print args.reference_conf + " not found"
                return False
            else:
                self.reference_conf = os.path.abspath(args.reference_conf)
        if args.cnv_call_time_job:
            self.cnv_call_time_job = args.cnv_call_time_job
        if args.cnv_threads:
            self.cnv_threads = args.cnv_threads

        return True


    def storePipelineParameters(self,args):
        """Updates general parameters to the map of parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.run_remove_duplicates_dict["run_pcr_duplicates"] = self.run_pcr_duplicates
        self.run_remove_duplicates_dict["input_data_se"] = self.input_data_se
        self.all_parameters ["Pipeline"] = self.run_remove_duplicates_dict

    def storeBasicMappingParameters(self,args):
        """Updates basic mapping parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.bwa_mapping_dict["bwa_reference"] = self.bwa_reference 
        self.bwa_mapping_dict["picard_path"] = self.picard_path 
        self.bwa_mapping_dict["bwa_threads"] = self.bwa_threads
        self.bwa_mapping_dict["tmp_folder"] = self.tmp_folder
        self.bwa_mapping_dict["bwa_time_job"] = self.bwa_time_job
        self.bwa_mapping_dict["java_heap"] = self.java_heap
        self.bwa_mapping_dict["mark_dup_time_job"] = self.mark_dup_time_job
        self.bwa_mapping_dict["mark_dup_threads"] = self.mark_dup_threads 
        self.all_parameters ["Bwa"] = self.bwa_mapping_dict

    def storeChoppingParameters(self,args):
        """Updates basic mapping parameters to be store in a JSON file
        args -- set of parsed arguments
        """
        self.chunks_chop_dict["chunks"] = self.chunks
        self.chunks_chop_dict["kmer_length"] = self.kmer_length
        self.chunks_chop_dict["windowing"] = self.windowing
        self.chunks_chop_dict["first_position"] = self.first_position
        self.chunks_chop_dict["bam_to_fastq_time_job"] = self.bam_to_fastq_time_job
        self.chunks_chop_dict["bam_to_fastq_threads"] = self.bam_to_fastq_threads
        self.chunks_chop_dict["frag_and_chop_time_job"] = self.frag_and_chop_time_job
        self.chunks_chop_dict["fragmentation_threads"] = self.fragmentation_threads
        self.all_parameters ["Fragment"] = self.chunks_chop_dict

    def storeCnvMappingParameters(self,args):
        """Updates cnv mapping parameters to be store in a JSON file
        args -- set of parsed arguments
        """

        self.cnv_mapping_dict["gem_index"] = self.gem_index 
         
        self.cnv_mapping_dict["mismatches"] = self.mismatches
        self.cnv_mapping_dict["max_edit_distance"] = self.max_edit_distance
        self.cnv_mapping_dict["strata_after_best"] = self.strata_after_best
        self.cnv_mapping_dict["max_decoded_matches"] = self.max_decoded_matches
        self.cnv_mapping_dict["min_decoded_strata"] = self.min_decoded_strata
        self.cnv_mapping_dict["gem_threads"] = self.gem_threads
        self.cnv_mapping_dict["gem_time_job"] = self.gem_time_job
        self.all_parameters ["CnvMapping"] = self.cnv_mapping_dict

    def storeCnvCallingParameters(self,args):
        """Updates cnv calling parameters to be store in a JSON file
        args -- set of parsed arguments
        """

        self.cnv_calling_dict["reference_conf"] = self.reference_conf  
        self.cnv_calling_dict["cnv_call_time_job"] = self.cnv_call_time_job
        self.cnv_calling_dict["cnv_threads"] = self.cnv_threads
        self.all_parameters ["CnvCalling"] = self.cnv_calling_dict


#1.Create object class Configuration File
configManager = CreateConfigurationFile()

#2.Create object for argument parsinng
parser = argparse.ArgumentParser(prog="createConfigurationFile",description="Create a configuration json file for the cnv pipeline.")
#2.1 Updates arguments and parsing
configManager.register_parameter(parser)
args = parser.parse_args()
#2.2 Check Parameters
if not configManager.check_parameters(args):
    sys.exit()

#3. store arguments to super map structure
configManager.storePipelineParameters(args)
configManager.storeBasicMappingParameters(args)
configManager.storeChoppingParameters(args)
configManager.storeCnvMappingParameters(args)
configManager.storeCnvCallingParameters(args)


#4. Store JSON file
with open(args.json_file, 'w') as of:
    json.dump(configManager.all_parameters, of, indent=2)



                              
                                        
