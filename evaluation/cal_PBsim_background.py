#! /usr/bin/python

import sys, os
import subprocess
import Annotation_formats

# To enable importing from samscripts submodulew
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
SCRIPT_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, ".."))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts'))

from fastqparser import read_fastq


class Report:
    def __init__(self):
        self.Total_expected_exons = 0
        self.Total_reads = 0
        self.Total_bases = 0
        self.Total_level2_reads = 0
        self.Total_level2_r_reads = 0
        self.Total_level3_SS_reads = 0
        self.Total_level3_AS_reads = 0

        self.Total_level2_bases = 0
        self.Total_level2_r_bases = 0
        self.Total_level3_SS_bases = 0
        self.Total_level3_AS_bases = 0

        self.Total_level2_expected_exons = 0
        self.Total_level2_r_expected_exons = 0
        self.Total_level3_SS_expected_exons = 0
        self.Total_level3_AS_expected_exons = 0


def processData(datafolder, annotationfile, ss_list, hq_read_file):

    #load annotation:
    annotations = Annotation_formats.Load_Annotation_From_File(annotationfile)
    annotation_dict = {}
    for annotation in annotations:
        if annotation.transcriptname in annotation_dict:
            pass
        else:
            # print(annotation.transcriptname)
            annotation_dict[annotation.transcriptname] = annotation

    #cal file count
    fFile = os.listdir(datafolder)
    file_count = int(len(fFile) / 2)

    SS_list = list()
    with open(ss_list, 'r') as f_ss:
        for line in f_ss:
            SS_list.append(line.strip())

    RR_list = list()
    with open(hq_read_file+".txt", 'r') as f_ss:
        for line in f_ss:
            RR_list.append(line.strip())
    # print(RR_list[0])

    RL_list = list()
    with open(hq_read_file+"_len.txt", 'r') as f_ss:
        for line in f_ss:
            RL_list.append(line.strip())
    # print(RL_list[0])

    report = Report()

    simFileSuffix = 'SimG2_S'

    for i in range(file_count):
        simFileName = simFileSuffix + '_%04d' % (i + 1)
        # print(simFileName)
        simRefFileName = simFileName + '.ref' #SimG1_S_0001.ref
        simMafFileName = simFileName + '.maf' #SimG1_S_0001.maf

        simFilePath = datafolder
        simRefFilePath = os.path.join(simFilePath, simRefFileName)
        simMafFilePath = os.path.join(simFilePath, simMafFileName)

        if not os.path.exists(simRefFilePath):
            raise Exception('Reference file for simulated read %s does not exist!' % simRefFilePath)
        if not os.path.exists(simMafFilePath):
            raise Exception('Sequence alignment (MAF) for simulated read %s does not exist!' % simMafFilePath)

        # Reading reference file
        [headers, seqs, quals] = read_fastq(simRefFilePath)
        simGeneName = headers[0]
        # if "transcript" in simGeneName:
        #     simGeneName = simGeneName.split(':')[1]
        annotation = annotation_dict[simGeneName]       # Getting the correct annotation
        
        maf_startpos = maf_length = 0
        i = 0
        len_i = 0
        l_c = 0
        sigA = False
        total_sim_bases = 0
        total_sim_exons = 0
        with open(simMafFilePath, 'rU') as maffile:
            for line in maffile:
                if line[0] == 's':
                    if line.split()[1] == 'ref': # sim ref
                        line_sim = maffile.readline()
                        line_sim_name = "SimG2_"+line_sim.split()[1]
                        # print(line_sim_name)
                        if line_sim_name in RR_list:

                            flag_wrong = 0
                            l_c += 1
                            elements = line.split()
                            maf_startpos = int(elements[2]) #0
                            maf_length = int(elements[3])   #3490
                            maf_reflen = int(int(elements[5])/3)   #3675

                            # Calculating expected partial alignmetns from MAF and annotations
                            #IMPORTANT:  if the reads were generated from an annotation on reverse strand, expected partial alignments must be reversed
                            if annotation.strand == Annotation_formats.GFF_STRANDRV:
                                maf_startpos = maf_reflen*3 - maf_length - maf_startpos
                                if maf_startpos > maf_reflen*2:
                                    maf_startpos = maf_startpos - maf_reflen*2
                                elif maf_startpos >maf_reflen:
                                    maf_startpos = maf_startpos - maf_reflen
                            # 1. Calculating the index of the first exon
                            # i - the index of exon currently being considered

                            i = 0
                            # print(annotation.items[i].getLength())
                            # print(maf_startpos)
                            while annotation.items[i].getLength() <= maf_startpos:
                                maf_startpos -= annotation.items[i].getLength()
                                i += 1
                                if len(annotation.items) == i:
                                    flag_wrong = 1
                                    break
                            if flag_wrong == 1:
                                continue
                            # Calculating expected partial alignments by filling up exons using maf_length
                            # maf_length = 
                            expected_partial_alignments = []
                            maf_length = int(maf_length/3)
                            while maf_length > 0:
                                # print(i)
                                start = annotation.items[i].start + maf_startpos
                                end = annotation.items[i].end
                                assert start <= end

                                # OLD: length = end-start+1
                                # KK: End is already indicating position after the last base, so adding one when callculating length is not correct
                                length = end - start
                                # print(length)
                                # print(maf_length)
                                if length <= maf_length:
                                    expected_partial_alignments.append((start, end))
                                    maf_length -= length
                                    i += 1
                                    if len(annotation.items) == i:
                                        maf_length = 0
                                else:
                                    expected_partial_alignments.append((start, start + maf_length))
                                    maf_length = 0
                                    i += 1

                                # Start position should only be considered for the first exon
                                maf_startpos = 0

                            report.Total_expected_exons += len(expected_partial_alignments)
                            num = len(expected_partial_alignments)
                            total_sim_exons += num
                            #level2
                            for ele in expected_partial_alignments[1:-1]:
                                if ele[1] - ele[0] < 30:
                                    report.Total_level2_reads += 1
                                    report.Total_level2_expected_exons += num
                                    sigA = True
                                    break
                            if sigA == False:
                                report.Total_level2_r_reads += 1
                                report.Total_level2_r_expected_exons += num
                    # else: #sim read
                            sim_bases = int(RL_list[len_i])
                            len_i += 1
                            report.Total_bases += sim_bases
                            total_sim_bases += sim_bases
                            #level2
                            if sigA == True:
                                report.Total_level2_bases += sim_bases
                                sigA = False
                            else:
                                report.Total_level2_r_bases += sim_bases
            #level3
            #print simGeneName
            if simGeneName in SS_list:
                report.Total_level3_SS_reads += l_c
                report.Total_level3_SS_bases += total_sim_bases
                report.Total_level3_SS_expected_exons += total_sim_exons
            else:
                report.Total_level3_AS_reads += l_c
                report.Total_level3_AS_bases += total_sim_bases
                report.Total_level3_AS_expected_exons += total_sim_exons

    report.Total_reads = report.Total_level3_SS_reads + report.Total_level3_AS_reads
    # print(report.Total_reads, report.Total_bases, report.Total_expected_exons)
    # print(report.Total_level2_reads, report.Total_level2_bases, report.Total_level2_expected_exons)
    # print(report.Total_level2_r_reads, report.Total_level2_r_bases, report.Total_level2_r_expected_exons)
    # print(report.Total_level3_AS_reads, report.Total_level3_AS_bases, report.Total_level3_AS_expected_exons)
    # print(report.Total_level3_SS_reads, report.Total_level3_SS_bases, report.Total_level3_SS_expected_exons)
    return report


def process(data_folder, group_list, annotationfile, ss_list, hq_read_file):
    group_path = data_folder + group_list #../simulation/ONT2D/group1
    report = processData(group_path, annotationfile, ss_list, hq_read_file)
    return report


if __name__ == '__main__':
    data_folder = sys.argv[1]
    group_list = sys.argv[2]
    annotationfile = sys.argv[3]
    ss_list = sys.argv[4]

    g_list = []
    with open(group_list, 'r') as g:
        lines = g.readlines()
        for line in lines:
            g_list.append(data_folder + group_folder)

    report_total = Report()
    for group_path in g_list:
        report = processData(group_path, annotationfile, ss_list, hq_read_file)
        #merge
        merge_static(report_total, report)
