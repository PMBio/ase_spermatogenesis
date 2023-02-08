import pysam
import sys
import gzip
import os


# codes used by pysam for aligned read CIGAR strings
BAM_CMATCH     = 0 # M
BAM_CINS       = 1 # I
BAM_CDEL       = 2 # D
BAM_CREF_SKIP  = 3 # N
BAM_CSOFT_CLIP = 4 # S
BAM_CHARD_CLIP = 5 # H
BAM_CPAD       = 6 # P
BAM_CEQUAL     = 7 # =
BAM_CDIFF      = 8 # X

BAM_CIGAR_DICT = {0 : 'M',
                  1 : 'I',
                  2 : 'D',
                  3 : 'N',
                  4 : 'S',
                  5 : 'H',
                  6 : 'P',
                  7 : '=',
                  8 : 'X'}

# this script processes a cellranger possorted bam file with snp tags into two independent matrices, one for the reference and alternative alleles
# the output should be a expression matrix in long format, i.e. a barcode and feature file and a list with barcode, cell, nReads entries
# there's a couple of open questions how cellranger does the processing, update while we figure that out
# for now, we'll filter out reads that 
# have duplicate UMIs (shouldnt be many, count that) # do that later, no so easy to do efficiently - our duplicate rate should be very low rn
# aren't in a pre-defined barcode derived from processing the full dataset (first) - done
# mapq and flag filters?
# are unclear wrt snp status (one snp supports maternal/one paternal allele - also shouldnt be many, count)

# this is a modification of the original script that counts across the whole sample, coverage per SNP

# input files
input_dir=sys.argv[1]
input_filename=input_dir + "/outs/allele_specific_counts/possorted_genome_bam.keep.bam.bam"
input_bam=pysam.Samfile(input_filename, "r")

print("Reading real cell barcodes...")

real_barcodes=gzip.open(input_dir + "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", "rb")
real_barcodes_arr=[str(line.rstrip().decode("utf-8")) for line in real_barcodes.readlines()]

# output files
output_dir=input_dir+"/outs/allele_specific_counts/ase_feature_matrix/"
if not os.path.exists(output_dir):
        os.makedirs(output_dir)

#barcodes_ref=open(output_dir + "./barcodes.reference.tsv", "w+")
#barcodes_alt=open(output_dir + "./barcodes.alternative.tsv", "w+")

#genes_ref=open(output_dir + "./genes.reference.tsv", "w+")
#genes_alt=open(output_dir + "./genes.alternative.tsv", "w+")

#matrix_ref=open(output_dir + "./matrix.reference.tsv", "w+")
#matrix_alt=open(output_dir + "./matrix.alternative.tsv", "w+")

output=open(output_dir + "./per_snp_coverage.tsv", "w+")

# hash of hashes to save umis against barcodes
# structure: {cell_bc1 = {gene1 = [umi1, ....]}

# structure of results: 
# {gene = {snp1: [nref, nalt]}

gene_snp_output = {}

#cell_barcode_ref = {bc : {} for bc in real_barcodes_arr}
#cell_barcode_alt = {bc : {} for bc in real_barcodes_arr}

# count statistics
read_count=0
reads_no_barcode=0
reads_ref=0
reads_alt=0
reads_conflict=0
duplicate_umis=0
no_umi=0

# current chromosome - process the umi hash after each chromosome to reduce size

print("Starting to process reads...")

def process_read(read, gene_snp_output):

    global duplicate_umis
    global no_umi

    gene_ensembl = read.get_tag("GX")

    if ";" in gene_ensembl:
        return(gene_snp_output)

    cell_barcode = read.get_tag("CB")

    try:
        umi = read.get_tag("UB")
    except KeyError:
        no_umi += 1
        return(gene_snp_output)

    try:
        gene_hash = gene_snp_output[gene_ensembl]
    except KeyError:
        gene_hash = {}
        #gene_snp_output[gene_ensembl] = []

    position_offset = 0
    discontinue_read = False

    # Get position offset incase part of read doesnt align: 
    for cigar in read.cigar:
        op = BAM_CIGAR_DICT[cigar[0]] # CIGAR 'operation'
        op_len  = cigar[1] # length of operation
        if op != "M" and op != "S":
            discontinue_read = True
        if op != "M":
            position_offset += op_len

    if discontinue_read:
        return(gene_snp_output)

    allele_tag = read.get_tag("ST").split("__")

    position_tags = []

    for snp in allele_tag:
        split_allele_tag = snp.split(":")
        position = int(split_allele_tag[0]) + read.reference_start - position_offset + 1

        if split_allele_tag[3] == "reference_base":
            which_base = 0
        else: 
            which_base = 1

        if position in gene_hash.keys():
            gene_hash[position][which_base] += 1
        else:
            if which_base == 0:
                new_counts = [1, 0]
            else:
                new_counts = [0, 1]
            gene_hash[position] = new_counts

    gene_snp_output[gene_ensembl] = gene_hash

    return(gene_snp_output)


def process_hash(barcode_hash):
    # some genes might already be processed in this way
    for cell in barcode_hash:
        for gene in barcode_hash[gene]:
            if not type(gene) == "numeric":
                barcode_hash[cell][gene] = len(barcode_hash[cell][gene])

    return(barcode_hash)


def process_allele_tag(allele_tag, mapping_position):

    for snp in allele_tag:
        split_allelle_tag = allele_tag.split(":")
        position = split_allele_tag[0]
        if split_allele_tag[3] == "reference_base":
            which_base = 0
        else: 
            which_base = 1


for read in input_bam:

    read_count+=1

    # this collapses the umi-list to read counts for the previous chromosome - should avoid creating too massive objects
    #if read.get_reference_name() != current_chromosome:
#       cell_barcode_ref = process_hash(cell_barcode_ref)
#       cell_barcode_alt = process_hash(cell_barcode_alt)

    # check if read barcode is in list of valid barcodes
    # this might be inefficient..


    try:
        if not read.get_tag("CB") in real_barcodes_arr:
            reads_no_barcode += 1
            continue
    except KeyError:
        reads_no_barcode += 1
        continue # this should give a key error if cell barcode is not set, in that case it also just continue

    try:
        gene = read.get_tag("GX")
    except KeyError:
        continue


    gene_snp_output = process_read(read, gene_snp_output)

    #n_reference_alleles = len([s for s in allele_tag if "reference" in s])
    #n_alternative_alleles = len([s for s in allele_tag if "alternative" in s])

    #if (n_reference_alleles > 0) & (n_alternative_alleles == 0):
#       # write to reference files
#       reads_ref += 1
#       cell_barcode_ref = process_read(read, cell_barcode_ref)

#   elif (n_reference_alleles == 0) & (n_alternative_alleles > 0):
#       # write to alternative files
#       reads_alt += 1
#       cell_barcode_alt = process_read(read, cell_barcode_alt)
#   else:
#       reads_conflict += 1

    #if read_count > 50000:
#       break

# now process hashes into read count files

#unique_genes_ref={}
#unique_genes_alt={}

for gene in gene_snp_output:
    for snp in gene_snp_output[gene]:
        gene = gene
        snp = snp
        count_ref = gene_snp_output[gene][snp][0]
        count_alt = gene_snp_output[gene][snp][1]
        output.write(gene + "\t" + str(snp) + "\t" + str(count_ref) + "\t" + str(count_alt) + "\n")
    
#for gene in unique_genes_ref.keys():
#   genes_ref.write(gene + "\n")

#for cell in real_barcodes_arr:
#   barcodes_ref.write(str(cell) + "\n")

# same for alternative alleles
#for cell in cell_barcode_alt:
#   for gene in cell_barcode_alt[cell]:
#       unique_genes_alt[gene] = 0
#       gene_count=len(cell_barcode_alt[cell][gene])
#       matrix_alt.write(cell + "\t" + gene + "\t" + str(gene_count) + "\n")

#for gene in unique_genes_alt.keys():
#   genes_alt.write(gene + "\n")


#print("Total reads %s" % read_count)
#print("No barcode reads %s" % reads_no_barcode)
#print("Reads reference %s" % reads_ref)
#print("Reads alternative %s" % reads_alt)
#print("Reads conflict %s" % reads_conflict)
#print("UMI-duplicates %s" % duplicate_umis)
#print("No UMI %s" % no_umi)

#%%MatrixMarket matrix coordinate integer general
#%metadata_json: {"format_version": 2, "software_version": "3.1.0"}
#54532 4378 8127238


input_bam.close()
real_barcodes.close()
output.close()
