import glob
import math
import numpy as np
import os
import subprocess


def bash(cmd, run_desc):
	return_str = ''
	''' Call the bash command passed as the argument

	Parameters:
	-----------
	cmd: String
		the command to be passed to the shell
	run_desc: String
		the description to be printed while command is running
	'''

	return_str += (run_desc+'\n')
	subprocess.call(cmd, shell=True, executable='/bin/bash')
	return (return_str)


def process(name, index):
	''' Call all bash scripts together from this wrapper function

	Parameters:
	-----------
	name: String
		name of the sra file
	index: String
		name for the index we create (not important)
	'''

	return_str = ''
	dir_ = '/home/vcm/'

	return_str += bash(
		'prefetch ' + name, 'Prefetching read metadata')
	return_str += bash(
		'sudo '+dir_+'sratoolkit.2.9.2-centos_linux64/bin/fastq-dump \
		-O . --split-spot '+name, 'Running fastq-dump')

	return_str += bash(
		'bwa index -p '+index+' -a bwtsw '+dir_+'wg.fa', 'Index reference')

	return_str += bash(
		'cp '+dir_+'site/'+name+'.fastq '+dir_, '')
	for ext in ['amb', 'ann', 'bwt', 'pac', 'sa']:
		return_str += bash(
			'cp '+dir_+'site/'+index+'.'+ext+' '+dir_, '')

	return_str += bash(
		'bwa aln -t 2 '+index+' '+dir_+name+'.fastq > '+dir_+name+'.fastq.bwa',
		'Aligning using bwa')
	return_str += bash(
		'bwa samse '+index+' '+dir_+name+'.fastq.bwa '+dir_+name+'.fastq > \
		'+dir_+name+'.sam', 'Obtaining results in SAM format')

	return_str += bash(
		'samtools view -S -b '+dir_+name+'.sam > '+dir_+name+'.bam',
		'Converting to BAM')
	return_str += bash(
		'samtools sort '+dir_+name+'.bam '+dir_+name+'-sorted',
		'Sorting BAM file')
	return_str += bash(
		'samtools faidx /home/vcm/wg.fa', 'Index using samtools')
	return_str += bash(
		'samtools mpileup -g -f '+dir_+'wg.fa '+dir_+name+'-sorted.bam > \
		'+dir_+name+'-raw.bcf', 'Use mpileup to generate vcf format')

	return_str += bash(
		'bcftools view -bvcg '+dir_+name+'-raw.bcf > \
		'+dir_+name+'-var.bcf', 'Call SNPs')
	return_str += bash(
		'bcftools view '+dir_+name+'-var.bcf | vcfutils.pl varFilter - > \
		'+dir_+name+'-var-final.vcf', 'Filter SNPs')

	print('all scripts done')

	return (return_str)


def check_high_qual(char, qual_threshold):
	''' Check if a quality score crosses the specified quality qual_threshold

	Parameters:
	-----------
	char: character
		quality score to be evaluated
	qual_threshold: integer
		the definition of a high quality read
	'''

	qual_str = '!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI'
	high_qual_str = qual_str[qual_threshold:]

	if char in high_qual_str:
		return(True)

	return(False)


def count_high_qual(name, len_threshold, qual_threshold):
	''' Count the number of high quality pairs in the dataset

	Parameters:
	-----------
	name: String
		name of the sra file
	len_threshold: integer
		the minimum length of a high quality subsequence in a sequence
	qual_threshold: integer
		the definition of a high quality read
	'''

	print('entered high quality')
	line_count = 0
	high_count = 0

	# open file
	fname = '/home/vcm/' + name + '.fastq'
	with open(fname, 'r') as f:
		for line in f:

			# ignore every first three of four lines
			if line_count < 3:
				line_count += 1
				continue

			# work with the fourth line
			if line_count == 3:
				line_count = 0
				high_qual_len = 0

				# iterate through fourth line
				for char in line:

					if check_high_qual(char, qual_threshold):
						high_qual_len += 1
					else:
						high_qual_len = 0

					if high_qual_len >= len_threshold:
						high_count += 1
						break

	return_str = 'High quality reads: ' + str(high_count) + '\n'
	cov = round((high_count*1.0)/130000, 2)
	return_str += 'Assuming 13 million bases, coverage %: ' + str(cov) + '\n'
	return(return_str)


def fill_bins(name):
	''' Initialize empty bin counts for each bins

	Parameters:
	-----------
	name: String
		the name of the sra
	'''

	bin_x_counts = dict()
	max_chr_pos = dict()
	fname = '/home/vcm/' + name + '-var-final.vcf'
	with open(fname, 'r') as f:
		for line in f:
			if line[0] == '#':
				continue

			line_split = line.split('\t')
			if line_split[0] not in max_chr_pos:
				max_chr_pos[line_split[0]] = int(line_split[1])
				continue
			else:
				max_chr_pos[line_split[0]] = max(
					int(line_split[1]),	max_chr_pos[line_split[0]])

	for chrom in max_chr_pos:
		for i in range((max_chr_pos[chrom]/10000)+1):
			bin_x_counts[chrom+':'+str(i*10000)+'-'+str((i+1)*10000)] = 0

	return(bin_x_counts)


def poisson(k, mean):
	''' Calculate poisson for X=x

	Parameters:
	-----------
	k: int
		the integer for which the poisson is to be calculated
	mean: int
		the mean/std dev for the distribution
	'''

	return(((math.exp(-mean))*(mean**k))/(1.0*math.factorial(k)))


def sum_poisson(num):
	''' Call the poisson for all number upto the provided num

	Parameters:
	-----------
	num: int
		the number till which the poissons are to be calculated
	'''
	sum_poisson = 0
	for i in range(num+1):
		sum_poisson += poisson(i, num)
	return sum_poisson


def calc_snp_types(name, print_num=False, print_type=False,
		print_window=False, print_stats=False, print_exp_obs=False):
	''' Calculate the number of SNPs and number of each type

	Parameters:
	-----------
	name: String
		the name of the sra
	print_num: boolean
		defines whether number of SNPs is to be printed
	print_type: boolean
		defines whether counts of types of indels are to be printed
	print_window: boolean
		defines whether counts of windows of indels are to be printed
	print_stats: boolean
		defines whether mean, std dev of windowed counts are to be printed
	print_exp_obs: boolean
		defines whether expected vs observed number of windows is to be printed
	'''

	return_str = ''

	net_count = 0
	individual_counts = {
		'A -> C': 0, 'A -> T': 0, 'A -> G': 0, 'T -> C': 0,
		'T -> A': 0, 'T -> G': 0, 'C -> G': 0, 'C -> A': 0,
		'C -> T': 0, 'G -> C': 0, 'G -> A': 0, 'G -> T': 0}

	bin_snp_counts = fill_bins(name)

	fname = '/home/vcm/' + name + '-var-final.vcf'
	with open(fname, 'r') as f:
		for line in f:
			if line[0] == '#':
				continue

			line_split = line.split('\t')
			if line_split[7][:2] == 'DP':
				net_count += 1
				z = int(line_split[1])/10000
				b = line_split[0]+':'+str(z*10000)+'-'+str((z+1)*10000)
				bin_snp_counts[b] += 1
				if ',' in line_split[4]:
					for each in line_split[4].split(','):
						individual_counts[line_split[3]+' -> '+each] += 1
				else:
					individual_counts[line_split[3]+' -> '+line_split[4]] += 1

	if print_num:
		return_str += ('\nTotal number of SNPs: ' + str(net_count) + '\n')

	if print_type:
		for snp in individual_counts:
			return_str += (
				'Number of ' + snp + ' SNPs is: ' +
				str(individual_counts[snp]) + '\n')

	if print_window:
		return_str += ('SNPs in each 10 kb window are as follows:\n')

		for key in sorted(bin_snp_counts.iterkeys()):
			return_str += ("%s: %s" % (key, bin_snp_counts[key]) + '\n')

	if print_stats:
		mean = np.mean(bin_snp_counts.values())
		std = np.std(bin_snp_counts.values())
		return_str += ('\nMean of SNPs in each bin is: %.2f \n' % mean)
		return_str += ('Std dev of SNPs in each bin is: %.2f \n' % std)

	if print_exp_obs:
		mean = np.mean(bin_snp_counts.values())
		std = np.std(bin_snp_counts.values())
		expected_ge_4_devs = 1-sum_poisson(int(4*math.sqrt(mean)))
		return_str += '\nExpected number of windows with SNPs 4 sigma greater \
			than mean is: ' + str(int(expected_ge_4_devs*len(bin_snp_counts))) + '\n'
		obtained_ge_4_devs = [x for x in bin_snp_counts.values() if x > mean]
		return_str += 'Obtained number of windows with SNPs 4 sigma greater than \
			mean is: ' + str(len(obtained_ge_4_devs)) + '\n'

	return (return_str)


def calc_indel_types(name, print_num=False, print_type=False,
		print_window=False, print_stats=False, print_exp_obs=False):
	''' Calculate the number of single base indels and number of each type

	Parameters:
	-----------
	name: String
		the name of the sra
	print_num: boolean
		defines whether number of SNPs is to be printed
	print_type: boolean
		defines whether counts of types of indels are to be printed
	print_window: boolean
		defines whether counts of windows of indels are to be printed
	print_stats: boolean
		defines whether mean, std dev of windowed counts are to be printed
	print_exp_obs: boolean
		defines whether expected vs observed number of windows is to be printed
	'''

	return_str = ''

	def count_changer(str1, str2):
		''' Update frequencies of each character

		Parameters:
		-----------
		str1: String
			Final sequence
		str2: String
			Original sequence
		'''
		if count_occurences(str1, 'A') - count_occurences(str2, 'A') == 1:
			individual_indel_counts['A+'] += 1
		if count_occurences(str1, 'A') - count_occurences(str2, 'A') == -1:
			individual_indel_counts['A-'] += 1
		if count_occurences(str1, 'T') - count_occurences(str2, 'T') == 1:
			individual_indel_counts['T+'] += 1
		if count_occurences(str1, 'T') - count_occurences(str2, 'T') == -1:
			individual_indel_counts['T-'] += 1
		if count_occurences(str1, 'G') - count_occurences(str2, 'G') == 1:
			individual_indel_counts['G+'] += 1
		if count_occurences(str1, 'G') - count_occurences(str2, 'G') == -1:
			individual_indel_counts['G-'] += 1
		if count_occurences(str1, 'C') - count_occurences(str2, 'C') == 1:
			individual_indel_counts['C+'] += 1
		if count_occurences(str1, 'C') - count_occurences(str2, 'C') == -1:
			individual_indel_counts['C-'] += 1

	def count_occurences(str1, char):
		''' Count occurence of a character in the string

		Parameters:
		-----------
		str1: String
			Sequence in which frequency is to be calculated
		char: String
			character whose frequency is to be calculated
		'''
		count = 0
		for ch in str1:
			if ch == char:
				count += 1
		return(count)

	net_indel_count = 0
	individual_indel_counts = {
		'A+': 0, 'A-': 0, 'T+': 0, 'T-': 0, 'G+': 0, 'G-': 0, 'C+': 0, 'C-': 0}

	bin_indel_counts = fill_bins(name)

	fname = '/home/vcm/' + name + '-var-final.vcf'
	with open(fname, 'r') as f:
		for line in f:
			if line[0] == '#':
				continue

			line_split = line.split('\t')
			if line_split[7][:2] == 'IN':
				if ',' in line_split[4]:
					for each in line_split[4].split(','):
						if abs(len(each)-len(line_split[3])) != 1:
							continue
						net_indel_count += 1
						z = int(line_split[1])/10000
						b = line_split[0]+':'+str(z*10000)+'-'+str((z+1)*10000)
						bin_indel_counts[b] += 1
						count_changer(each, line_split[3])

				else:
					if abs(len(line_split[4])-len(line_split[3])) != 1:
						continue
					net_indel_count += 1
					z = int(line_split[1])/10000
					b = line_split[0]+':'+str(z*10000)+'-'+str((z+1)*10000)
					bin_indel_counts[b] += 1
					count_changer(line_split[4], line_split[3])

	if print_num:
		return_str += (
			'\nTotal number of single base indels:' +
			str(net_indel_count) + '\n')

	if print_type:
		for indel in individual_indel_counts:
			return_str += (
				'Number of ' + indel + ' indels is: ' +
				str(individual_indel_counts[indel]) + '\n')

	if print_window:
		return_str += ('Indels in each 10 kb window are as follows:\n')

		for key in sorted(bin_indel_counts.iterkeys()):
			return_str += ("%s: %s" % (key, bin_indel_counts[key]) + '\n')

	if print_stats:
		mean = np.mean(bin_indel_counts.values())
		std = np.std(bin_indel_counts.values())
		return_str += ('\nMean of single base indels in each bin is: %.3f \n' % mean)
		return_str += ('Std dev of single base indels in each bin is: %.3f \n' % std)

	if print_exp_obs:
		mean = np.mean(bin_indel_counts.values())
		std = np.std(bin_indel_counts.values())
		expected_ge_4_devs = 1-sum_poisson(int(4*math.sqrt(mean)))
		return_str += '\nExpected number of windows with Indels 4 sigma greater \
			than mean is: ' + str(int(expected_ge_4_devs*len(bin_indel_counts))) + '\n'
		obtained_ge_4_devs = [x for x in bin_indel_counts.values() if x > mean]
		return_str += 'Obtained number of windows with Indels 4 sigma greater than \
			mean is: ' + str(len(obtained_ge_4_devs)) + '\n'

	return return_str


def delete_files(name, index):
	''' Delete files for fresh start

	name: String
		name of the sra file
	index: String
		name for the index we create (not important)
	'''

	sra = 'ncbi/public/sra/'
	for prefix in [name, index, sra+name, 'site/'+name, 'site/'+index]:
		for f in glob.glob('/home/vcm/' + prefix + '*'):
			os.remove(f)

	print('All already downloaded files deleted!')
	return ('All already downloaded files deleted\n')
