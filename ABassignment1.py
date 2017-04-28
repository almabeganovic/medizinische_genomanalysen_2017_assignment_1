#! /usr/bin/env python3

import mysql.connector
import sys
import pysam
import pybedtools
import subprocess

__author__ = 'Alma Beganovic'


class Assignment1:
    def __init__(self):
        ## my gene of interest
        self.gene = "BDNF"


        self.geneinfo = ['BDNF', 'NM_170731', 'chr11', 27676440, 27743605, '-', 2, b'27676440,27742958,', b'27680132,27743605']

    def fetch_gene_coordinates(self, genome_reference, file_name):
        print ("Connecting to UCSC to fetch data")

        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password',
                                      db=genome_reference)

        ## Get cursor
        cursor = cnx.cursor()

        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]

        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)

        ## Execute query
        cursor.execute(query)


        with open(file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")

        ## Close cursor & connection
        cursor.close()
        cnx.close()

        print ("Done fetching data")

    def get_sam_header(self):
        print("-sam_header:")
        cmd = ["samtools_0.1.18  view -H {}".format('HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam')]
        subprocess.call(cmd, shell=True)
        print ()



    # http://pysam.readthedocs.io/en/latest/api.html
    # For the function, .bam file must indexed (samtools index HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam )
    def get_properly_paired_reads_of_gene(self):
        samFile = pysam.AlignmentFile('HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam', "rb")
        line_count = 0
        print('The proper paired reads are:')
        for read in samFile.fetch("11", self.geneinfo[3], self.geneinfo[4]):
            if read.is_proper_pair:
                line_count += 1
                print(read)
        samFile.close()
        print('Count paired_reads_of_gene', line_count)


    def get_gene_reads_with_indels(self):
        samFile = pysam.AlignmentFile('HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam', "rb")
        line_count = 0
        print('Gene_reads_with_indels:')
        for read in samFile.fetch("11", self.geneinfo[3], self.geneinfo[4]):
            columns = str(read).split("\t")
            if "I" in str(columns[5]) or "D" in str(columns[5]):
                print("-gene_reads_with_indels:")  #for print indels please comment
                print(read)                        #for print indels please comment
                line_count += 1
        samFile.close()
        print('Count Reads_with_indels', line_count)

    def calculate_total_average_coverage(self):
        a = pybedtools.BedTool('HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam')
        cov = a.genome_coverage(bg=True)
        line_count = 0
        total_coverage = 0
        for line in cov:
            line_count += 1
            total_coverage += int(line[3])

        total_average_coverage = total_coverage / line_count
        print("-total_average_coverage:")
        print(total_average_coverage)

    def calculate_gene_average_coverage(self):
        a = pybedtools.BedTool('HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam')
        cover = a.genome_coverage(bg=True)
        line_count = 0
        total_coverage = 0
        for line in cover:
            if int(line[1]) >= int(self.geneinfo[3]) and int(line[2]) <= int(self.geneinfo[4]):
                line_count += 1
                total_coverage += int(line[3])

        gene_average_coverage = total_coverage / line_count
        print("-gene_average_coverage:")
        print(gene_average_coverage)

    def get_number_mapped_reads(self):
        samFile = pysam.AlignmentFile('HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam', "rb")
        n = 0
        for read in samFile:
            if not read.is_unmapped:
                n += 1
        if n == 0:
            print ("-number_nonmapped_reads:")
        else:
            print("-number_mapped_reads:")
            print(n)
        samFile.close()

    ## Gene Symbol and the first position will be printed out
    def get_gene_symbol(self):
        print("-gene_symbol:")
        print(self.geneinfo[0])


    ## Region of gene will be printed out
    def get_region_of_gene(self):
        print("-region_of_gene:")
        print("{} starts with {} and ends with {}".format(self.geneinfo[2], self.geneinfo[3], self.geneinfo[4]))
        print("-region_of_gene on {} starts with {} and ends with {}".format(self.geneinfo[2], self.geneinfo[3], self.geneinfo[4]))

    ## Number of exon and the first position will be printed out
    def get_number_of_exons(self):
        print("-number_of_exons:")
        print (self.geneinfo[6])


    def print_summary(self):
        print("All results:")
        self.get_sam_header()
        self.get_properly_paired_reads_of_gene()      # 3056
        self.get_gene_reads_with_indels()             # 30 run takes long
        self.calculate_total_average_coverage()       # 5.608294198569065 run takes long
        self.calculate_gene_average_coverage()        # 5.359311201846263 run takes long
        self.get_number_mapped_reads()                # 6396581
        self.get_gene_symbol()                        # BDNF
        self.get_region_of_gene()                     # chr11 starts with 27676440 and ends with 27743605
        self.get_number_of_exons()                    # 2


if __name__ == '__main__':
    print ("Assignment 1")
    assignment1 = Assignment1()
    assignment1.print_summary()
