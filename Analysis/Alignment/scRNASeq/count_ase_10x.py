import pysam
import sys
import gzip
import os

# this script processes a cellranger possorted bam file with snp tags into two independent matrices, one for the reference and alternative alleles
# the output should be a expression matrix in long format, i.e. a barcode and feature file and a list with barcode, cell, nReads entries

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

barcodes_ref=open(output_dir + "./barcodes.reference.tsv", "w+")
barcodes_alt=open(output_dir + "./barcodes.alternative.tsv", "w+")

genes_ref=open(output_dir + "./genes.reference.tsv", "w+")
genes_alt=open(output_dir + "./genes.alternative.tsv", "w+")

matrix_ref=open(output_dir + "./matrix.reference.tsv", "w+")
matrix_alt=open(output_dir + "./matrix.alternative.tsv", "w+")

# hash of hashes to save umis against barcodes
# structure: {cell_bc1 = {gene1 = [umi1, ....]}

cell_barcode_ref = {bc : {} for bc in real_barcodes_arr}
cell_barcode_alt = {bc : {} for bc in real_barcodes_arr}

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

def process_read(read, barcode_hash):

	global duplicate_umis
	global no_umi

	gene_ensembl = read.get_tag("GX")

	if ";" in gene_ensembl:
		return(barcode_hash)

	cell_barcode = read.get_tag("CB")

	try:
		umi = read.get_tag("UB")
	except KeyError:
		no_umi += 1
		return(barcode_hash)

	curr_bc = barcode_hash[cell_barcode]

	try:
		gene_hash = curr_bc[gene_ensembl]
		if umi in gene_hash:
			duplicate_umis += 1
			return(barcode_hash)
		gene_hash.append(umi)
		curr_bc[gene_ensembl] = gene_hash
	except KeyError:
		curr_bc[gene_ensembl] = [umi]

	barcode_hash[cell_barcode] = curr_bc

	return(barcode_hash)


def process_hash(barcode_hash):
	# some genes might already be processed in this way
	for cell in barcode_hash:
		for gene in barcode_hash[gene]:
			if not type(gene) == "numeric":
				barcode_hash[cell][gene] = len(barcode_hash[cell][gene])

	return(barcode_hash)


for read in input_bam:

	read_count+=1

	# this collapses the umi-list to read counts for the previous chromosome - should avoid creating too massive objects
	#if read.get_reference_name() != current_chromosome:
#		cell_barcode_ref = process_hash(cell_barcode_ref)
#		cell_barcode_alt = process_hash(cell_barcode_alt)

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

	allele_tag = read.get_tag("ST").split("__")

	n_reference_alleles = len([s for s in allele_tag if "reference" in s])
	n_alternative_alleles = len([s for s in allele_tag if "alternative" in s])

	if (n_reference_alleles > 0) & (n_alternative_alleles == 0):
		# write to reference files
		reads_ref += 1
		cell_barcode_ref = process_read(read, cell_barcode_ref)

	elif (n_reference_alleles == 0) & (n_alternative_alleles > 0):
		# write to alternative files
		reads_alt += 1
		cell_barcode_alt = process_read(read, cell_barcode_alt)
	else:
		reads_conflict += 1

	#if read_count > 50000:
#		break



# now process hashes into read count files

unique_genes_ref={}
unique_genes_alt={}

for cell in cell_barcode_ref:
	for gene in cell_barcode_ref[cell]:
		unique_genes_ref[gene] = 0
		gene_count=len(cell_barcode_ref[cell][gene])
		matrix_ref.write(cell + "\t" + gene + "\t" + str(gene_count) + "\n")

for gene in unique_genes_ref.keys():
	genes_ref.write(gene + "\n")

for cell in real_barcodes_arr:
	barcodes_ref.write(str(cell) + "\n")

# same for alternative alleles
for cell in cell_barcode_alt:
	for gene in cell_barcode_alt[cell]:
		unique_genes_alt[gene] = 0
		gene_count=len(cell_barcode_alt[cell][gene])
		matrix_alt.write(cell + "\t" + gene + "\t" + str(gene_count) + "\n")

for gene in unique_genes_alt.keys():
	genes_alt.write(gene + "\n")

for cell in real_barcodes_arr:
	barcodes_alt.write(str(cell) + "\n")


print("Total reads %s" % read_count)
print("No barcode reads %s" % reads_no_barcode)
print("Reads reference %s" % reads_ref)
print("Reads alternative %s" % reads_alt)
print("Reads conflict %s" % reads_conflict)
print("UMI-duplicates %s" % duplicate_umis)
#print("No UMI %s" % no_umi)

#%%MatrixMarket matrix coordinate integer general
#%metadata_json: {"format_version": 2, "software_version": "3.1.0"}
#54532 4378 8127238


input_bam.close()
real_barcodes.close()

barcodes_ref.close()
barcodes_alt.close()
genes_ref.close()
genes_alt.close()
matrix_ref.close()
matrix_alt.close()
