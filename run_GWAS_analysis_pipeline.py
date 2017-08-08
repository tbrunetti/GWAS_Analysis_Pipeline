from chunkypipes.components import *
import collections
import multiprocessing
import math
import datetime
import sys
import os
import time

class Pipeline(BasePipeline):
	def dependencies(self):
		return ['pandas', 'numpy', 'matplotlib', 'fpdf', 'Pillow', 'pypdf2', 'statistics', 'xlrd']

	def description(self):
		return 'Pipeline made for analyzing GWAS data after QC cleanup'

	def configure(self):
		return {
			'plink':{
				'path':'Full path to PLINK executable'
			},
			'king':{
				'path':'Full path to KING executable'
			},
			'thousand_genomes':{
				'path': 'Full path to PLINK BED file format LD pruned phase 3 1000 genomes file'
			},
			'R_libraries':{
				'path': 'Full path to directory where R libraries are stored'
			}
		}

	def add_pipeline_args(self, parser):
		# should other input options be available??
		parser.add_argument('-inputPLINK', required=True, type=str, help='Full path to PLINK file ending in .BED or .PED')
		parser.add_argument('-phenoFile', required=True, type=str, help='Full path to phenotype file, see argument readme for more details on format')
		parser.add_argument('--sampleRemoval', default=None, help='Full path for samples to remove before analysis (i.e. those with sex discrepenences, poor QC, etc...) see readme for more details on format')
		parser.add_argument('--outDir', default=os.getcwd(), type=str, help='[default=current working directory] Full path of existing directory to output results')
		parser.add_argument('--projectName', default=str(datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")), type=str, help='[default=date time stamp] Name of project')
		parser.add_argument('--startStep', default='hwe', type=str, help='The part of the pipeline you would like to start with; if user does not specify this parameter, begins analyis at Hardy-Weinberg equilibrium')
		parser.add_argument('--endStep', default=None, type=str, help='Point of the pipeline where you would like to stop analysis, if none specified, stops after PCA is done or if GENanalysis is specified as start, finishes after GENanalysis is done')
		parser.add_argument('--hweThresh', default=1e-6, help='Filters out SNPs that are smaller than this threshold due to liklihood of genotyping error')
		parser.add_argument('--LDmethod', default='indep', type=str, help='[default=indep, options:indep, indep-pairwise, indep-pairphase] Method to calculate linkage disequilibrium')
		parser.add_argument('--VIF', default=2, type=int, help='[default=2] variant inflation factor for indep method LD pruning')
		parser.add_argument('--rsq', default=0.50, type=float, help='[default=0.50] r squared threshold for indep-pairwise or indep-pairphase LD pruning method')
		parser.add_argument('--windowSize', default=50, type=int, help='[default=50] the window size in kb for LD analysis')
		parser.add_argument('--stepSize', default=5, type=int, help='[default=5] variant count to shift window after each interation')
		parser.add_argument('--maf', default=0.05, type=float, help='[default=0.05], filter remaining LD pruned variants by MAF')
		parser.add_argument('--hetThresh', default=0.10, type=float, help='[default=0.10], filter out samples where inbreeding coefficient is greater than threshold (heterozygosity filtering)')
		parser.add_argument('--hetThreshMin', default=-0.10, type=float, help='')
		parser.add_argument('--reanalyze', action='store_true', help='by adding this flag, it means you are going to pass a dataset through the pipeline that has already been partially/fully analyzed by this pipeline. WARNING! May over write exisiting data!!')
		parser.add_argument('--usePCs', default='pc1,pc2,pc3', type=str, help='[default: pc1,pc2,pc3], the user can pass any pc from 1-50 in a comma separated list to regress out in the GENanalysis step; must be used with --reanalyze flag; format should be the following: pc1,pc2,pc3,pc5')

	@staticmethod
	def check_steps(order, start, stop):
		if start == 'hwe' and stop == None:
			return order
		elif start != 'GENanalysis': # reanalyze flag should be used here
			index_of_start =order.index(start)
			if stop == None:
				return order[index_of_start:]
			else:
				index_of_end = order.index(stop)
				return order[index_of_start:index_of_end+1]
		elif start == 'GENanalysis': # reanalyze flag should be used here
			return ['GENanalysis']	


	# checks the type of PLINK file input
	@staticmethod
	def check_plink_format(plinkFile, plink):
		if plinkFile[-4:].lower() == '.ped':
			print "Input .ped, converting to binary"
			convert = subprocess.call([str(plink), '--file', str(plinkFile[:-4]), '--make-bed', '--out', str(plinkFile[:-4])])
		
		elif plinkFile[-4:].lower() == '.bed':
			print "Input seems to follow bed format"

		else:
			sys.exit("Error!! Input not recognized, please input .ped or .bed PLINK file" )


	# create files to input into plink for samples to keep
	@staticmethod
	def ethnic_plinks_lists(phenotype, plinkFileName, famFile, removeSamples, outDir):
		import pandas as pd

		phenotype_table = pd.read_excel(phenotype, sheetname="Sheet2", header=0, converters={'FID':str, 'IID':str}) 
		phenotype_table['IID'] = phenotype_table['IID'].str.strip() # remove whitespace for proper merging
		phenotype_table['FID'] = phenotype_table['FID'].str.strip() # remove whitespace for proper merging
		
		# update fam file with phenotype information
		original_fam = pd.read_table(famFile, delim_whitespace=True, names = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'], converters={'FID':str, 'IID':str, 'Gender':str}) 
		original_fam['IID'] = original_fam['IID'].str.strip() # remove whitespace for proper merging
		original_fam['FID'] = original_fam['FID'].str.strip() # remove whitespace for proper merging

		merged_dataframe = original_fam.merge(phenotype_table, how='left', on=['FID', 'IID']) # only merge information where FID and IID are available in the original fam file

		race_subsets = list(set(list(merged_dataframe['Race']))) # gets all racial groups listed in phenotype file
		race_subsets_cleaned = [i for i in race_subsets if str(i)!='nan'] # removes samples that do not have a race listed 
		make_plinks = {} # stores file name of samples to keep for each group
		

		for ethnic_group in race_subsets_cleaned:
			os.mkdir(outDir + '/' + '_'.join(ethnic_group.split()))
			keep_file = open(outDir + '/' + '_'.join(ethnic_group.split()) + '/' + plinkFileName +'_'+str('_'.join(ethnic_group.split())) + '_keepIDs.txt', 'w')
			
			if removeSamples != None:
				removeSamples_dataframe = pd.read_table(removeSamples, delim_whitespace=True, names=['FID', 'IID'])
				# get a list of sample IDs to remove
				sampleIDs_remove = list(removeSamples_dataframe['IID'])
				subset_by_race_only = merged_dataframe.loc[(merged_dataframe['Race'] == ethnic_group) & (pd.isnull(merged_dataframe['FID'].str.strip()) == False)]
				# need to check for last part of conditional in case FID is missing, remove the row...causes problems with PLINK
				subset_by_race = subset_by_race_only[(subset_by_race_only['IID'].isin(sampleIDs_remove) == False)]
				
			else:
				# need to check for last part of conditional in case FID is missing, remove the row...causes problems with PLINK
				subset_by_race = merged_dataframe.loc[(merged_dataframe['Race'] == ethnic_group) & (pd.isnull(merged_dataframe['FID'].str.strip()) == False)]

			subset_by_race[['FID', 'IID']].to_csv(keep_file.name, sep='\t', index=False, header=False) # format it FID <tab> IID <new line>
			merged_dataframe[['FID', 'IID', 'PAT', 'MAT', 'Gender', 'Phenotype']].to_csv(famFile, sep=' ', index=False, header=False)
			make_plinks['_'.join(ethnic_group.split())]=keep_file.name
			keep_file.flush()
			keep_file.close()

		# also create a directory to sort through merged ethnic results
		os.mkdir(outDir + '/merged_group_files')
		
		return make_plinks


	def run_pipeline(self, pipeline_args, pipeline_config):
		import pandas as pd
		sys.path.append(".")
		import subprocess
		import summary_stats
		import statistics as stats
		from fpdf import FPDF
		import PyPDF2

		# keeps track of all files that can be deleted after the pipeline finishes
		stage_for_deletion = []

		reduced_plink_name = pipeline_args['inputPLINK'].split('/')[-1][:-4] #only get plink file name not full absolute path and removes suffix
		
		# specifying output location and conflicting project names of files generated	
		try:
			if pipeline_args['reanalyze'] == True:
				outdir = pipeline_args['outDir']+'/'+pipeline_args['projectName']
				print "Reanalyzing data from an existing project"
			else:
				os.stat(pipeline_args['outDir']+'/'+pipeline_args['projectName'])
				sys.exit("project already exists!!")
		except:
			print "Making new directory called "+str(pipeline_args['projectName']) + ' located in ' + str(pipeline_args['outDir'])
			outdir = pipeline_args['outDir']+'/'+pipeline_args['projectName'] # new output directory
			os.mkdir(outdir)
			#self.settings.logger.set(
			#	destination=pipeline_args['outDir'] + '/' + pipeline_args['projectName'] +'/stdout.log',
			#	destination_stderr=pipeline_args['outDir'] + '/' + pipeline_args['projectName'] + '/stderr.log')



		step_order = self.check_steps(
			order = ['hwe', 'LD', 'maf', 'merge', 'het', 'ibd', '1000_genomes', 'KING', 'PCA'], 
			start = pipeline_args['startStep'],
			stop = pipeline_args['endStep']
			)
		

		#step_order = ['hwe', 'LD', 'maf', 'merge', 'het', 'ibd', '1000_genomes', 'KING', 'PCA'] # order of pipeline if full suite is used
		#step_order = ['hwe', 'LD', 'maf', 'merge', 'het', 'ibd', '1000_genomes', 'KING', 'PCA']
		#step_order = ['GENanalysis']
		
		# initialize PLINK and KING software
		general_plink = Software('plink', pipeline_config['plink']['path'])
		general_king = Software('king', pipeline_config['king']['path'])
		# make sure plink input format is correct
		self.check_plink_format(
			plinkFile = pipeline_args['inputPLINK'],
			plink = pipeline_config['plink']['path']
			)

		
		if pipeline_args['reanalyze'] == False:
			keep_files = self.ethnic_plinks_lists(
				phenotype = pipeline_args['phenoFile'],
				plinkFileName = reduced_plink_name,
				famFile = pipeline_args['inputPLINK'][:-4] + '.fam',
				removeSamples = pipeline_args['sampleRemoval'],
				outDir = outdir
				)


			# will make separate plink files for each ethnic group and use autosomes only
			for key, value in keep_files.iteritems():
				general_plink.run(
					Parameter('--bfile', pipeline_args['inputPLINK'][:-4]),
					Parameter('--keep', value),
					Parameter('--autosome'),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/' + str(key) + '/' + reduced_plink_name + '_' + str(key))
					)


		while len(step_order) != 0:
			
			# hardy-weinberg equilibrium filtering
			if step_order[0] == 'hwe':
				print "running HWE step"
				print multiprocessing.cpu_count()
				hwe_passing = {}
				for directories in os.listdir(outdir):
					if ((os.path.isdir(os.path.join(outdir, directories))) and (directories != 'merged_group_files')):

						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories),
							Parameter('--hardy'),
							Parameter('--hwe', pipeline_args['hweThresh']),
							Parameter('midp'),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_hweFiltered')
							)
						
						total_snps_analyzed_hwe = subprocess.check_output(['wc', '-l',  outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.bim'])
						total_snps_passing_hwe = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_hweFiltered.bim'])
						hwe_passing[directories] = [total_snps_analyzed_hwe.split()[0]] + [total_snps_passing_hwe.split()[0]] # store total analyzed and passing for hwe step
						hwe_passing[directories] = hwe_passing[directories] + [outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_hweFiltered.hwe']
				
				hwe_stats = summary_stats.hwe(dictHWE=hwe_passing, thresh=pipeline_args['hweThresh'], outDir = outdir)
				hwe_stats.output(outdir + '/' + pipeline_args['projectName'] + '_hweStats.pdf', 'F') # output results to PDF
				step_order.pop(0)

			
			# LD pruning
			elif step_order[0] == 'LD':
				print "running LD pruning step"
				ld_passing = {}

				if pipeline_args['LDmethod'] == 'indep': # gets lists of variants to keep and to prune using VIP
					for directories in os.listdir(outdir):
						if ((os.path.isdir(os.path.join(outdir, directories))) and (directories != 'merged_group_files')):
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
								Parameter('--'+pipeline_args['LDmethod'], str(pipeline_args['windowSize'])+'kb', str(pipeline_args['stepSize']), str(pipeline_args['VIF'])),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories)
								)
							
							total_snps_analyzed_ld = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered.bim'])
							total_snps_passing_ld = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.prune.in'])
							ld_passing[directories] = [total_snps_analyzed_ld.split()[0]] + [total_snps_passing_ld.split()[0]] # store total analyzed and passing for LD pruning step
							
							# creates new PLINK files with excluded variants removed
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
								Parameter('--exclude', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.prune.out'),
								Parameter('--make-bed'),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned')
								)
							
				
				else:
					for directories in os.listdir(outdir): # get lists of variants to keep and to prune using rsq
						if ((os.path.isdir(os.path.join(outdir, directories))) and (directories != 'merged_group_files')):
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
								Parameter(pipeline_args['LDmethod'], str(pipeline_args['windowSize'])+'kb', str(pipeline_args['stepSize']), str(pipeline_args['rsq'])),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories)
								)

							total_snps_analyzed_ld = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered.bim'])
							total_snps_passing_ld = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.prune.in'])
							ld_passing[directories] = [total_snps_analyzed_ld] + [total_snps_passing_ld] # store total analyzed and passing for LD pruning step
							
							# creates new PLINK files with excluded variants removed
							general_plink.run(
								Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
								Parameter('--exclude', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.prune.out'),
								Parameter('--make-bed'),
								Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned')
								)
					
				ld_stats = summary_stats.pruning(dictLD=ld_passing)
				ld_stats.output(outdir + '/' + pipeline_args['projectName'] + '_ldStats.pdf', 'F') # output results to PDF
				step_order.pop(0)

			# filters pruned variants by MAF
			elif step_order[0] == 'maf':
				print "running maf step"
				list_of_merge_maf_greater_thresh = open(outdir + '/' + reduced_plink_name + '_greater_mafs.txt', 'w')
				list_of_merge_maf_less_thresh = open(outdir + '/' + reduced_plink_name + '_small_mafs.txt', 'w')
				
				maf_passing = {}
				for directories in os.listdir(outdir):
					if ((os.path.isdir(os.path.join(outdir, directories))) and (directories != 'merged_group_files')):
						
						# filters out variants below set maf threshold
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned'),
							Parameter('--maf', str(pipeline_args['maf'])),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf')
							)

						total_snps_analyzed_maf = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned.bim'])
						total_snps_greater_maf = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf.bim'])
						maf_passing[directories] = [total_snps_analyzed_maf.split()[0]] + [total_snps_greater_maf.split()[0]]

						# stores name of file for future merging (must be space delimited ordered bed, bim, fam files)
						list_of_merge_maf_greater_thresh.write(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf.bed ' + 
							outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf.bim ' +
							outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf.fam' + '\n')

						# filters out variants set maf threshold
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned'),
							Parameter('--max-maf', str(pipeline_args['maf'])),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf')
							)
						
						total_snps_less_maf = subprocess.check_output(['wc', '-l', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf.bim'])
						maf_passing[directories] = maf_passing[directories] + [total_snps_less_maf.split()[0]]
	

						# stores name of file for future merging
						list_of_merge_maf_less_thresh.write(outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf.bed ' +
							outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf.bim ' +
							outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf.fam' +'\n')
					
					
				list_of_merge_maf_greater_thresh.flush() # push out buffer
				list_of_merge_maf_less_thresh.flush() # push out buffer
				
				maf_stats = summary_stats.minor_allele_freq(dictMAF=maf_passing, thresh=pipeline_args['maf'])
				maf_stats.output(outdir + '/' + pipeline_args['projectName'] + '_mafStats.pdf', 'F') # output results to PDF

				step_order.pop(0)
				
		

			# merge all sets of data together
			elif step_order[0] == 'merge':
				# merge all ethnic groups together with mafs greater than threshold
				print "merging all files"
				first_file_large = ''
				with open(list_of_merge_maf_greater_thresh.name, 'r') as base_file:
					remainder = base_file.read().splitlines(True)
					first_file_large = remainder[0]
				with open(list_of_merge_maf_greater_thresh.name, 'w') as merge_list:
					merge_list.writelines(remainder[1:])
					merge_list.flush()

				prefixFirstLarge = first_file_large.split()[0][:-4] 
				
				general_plink.run(
					Parameter('--bfile', prefixFirstLarge),
					Parameter('--merge-list', list_of_merge_maf_greater_thresh.name),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged')
					)

				# recalculate freq stats for merged set
				general_plink.run(
					Parameter('--bfile', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged'),
					Parameter('--freq'),
					Parameter('--out', outdir + '/merged_group_files/' + reduced_plink_name + '_freq_recalculated_greater_all_merged')
					)


				# merge all ethnic groups together with mafs smaller than threshold
				first_file_small = ''
				with open(list_of_merge_maf_less_thresh.name, 'r') as base_file:
					remainder = base_file.read().splitlines(True)
					first_file_small = remainder[0]
				with open(list_of_merge_maf_less_thresh.name, 'w') as merge_list:
					merge_list.writelines(remainder[1:])
					merge_list.flush()

				prefixFirstSmall = first_file_small.split()[0][:-4]
				
				general_plink.run(
					Parameter('--bfile', prefixFirstSmall),
					Parameter('--merge-list', list_of_merge_maf_less_thresh.name),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_less_thresh_all_ethnic_groups_merged')
					)


				# recalculate freq stats for merged set
				general_plink.run(
					Parameter('--bfile', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_less_thresh_all_ethnic_groups_merged'),
					Parameter('--freq'),
					Parameter('--out', outdir + '/merged_group_files/' + reduced_plink_name + '_freq_recalculated_less_all_merged')
					)

				step_order.pop(0)


			elif step_order[0] == 'het':
				print "checking heterozygosity"
				general_plink.run(
					Parameter('--bfile', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged'),
					Parameter('--het'),
					Parameter('--out', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged')
					)


				het_dataframe = pd.read_table(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged.het', delim_whitespace=True)
				samples_failing_het, het_pdf = summary_stats.heterozygosity(het_dataframe = het_dataframe, thresh = pipeline_args['hetThresh'], minThresh=pipeline_args['hetThreshMin'],outDir = outdir)
				
				general_plink.run(
					Parameter('--bfile', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged'),
					Parameter('--remove', samples_failing_het),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged')
					)
				
				het_pdf.output(outdir + '/' + pipeline_args['projectName'] + '_hetStats.pdf', 'F') # output results to PDF

				step_order.pop(0)
				


			# determine the pairwise relationship of all merged samples
			# DUPLICATES MUST BE REMOVED BEFORE USING GENESIS PIPELINE!
			elif step_order[0] == 'ibd':
				# only gets run on the bed file with MAF > specified tresh
				print "running PLINK ibd step"
				general_plink.run(
					Parameter('--bfile', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged'),
					Parameter('--genome'),
					Parameter('--read-freq', outdir + '/merged_group_files/' + reduced_plink_name + '_freq_recalculated_greater_all_merged.frq'),
					Parameter('--out', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged')
					)
				
				ibd_results = pd.read_table(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged.genome', delim_whitespace=True)
				relatedness_stats, remove_samples = summary_stats.relatedness(ibd_dataframe = ibd_results, outDir=outdir)

				# remove samples that are duplicates
				general_plink.run(
					Parameter('--bfile', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged'),
					Parameter('--remove', remove_samples),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed')
					)
				
				relatedness_stats.output(outdir + '/' + pipeline_args['projectName'] + '_relatedness.pdf', 'F') # output results to PDF
				step_order.pop(0)
	

			
			elif step_order[0] == '1000_genomes':
				
				keep_these_snps_in_1000 = open(outdir + '/all_passing_snps_from_project.txt', 'w')
				bim_file = pd.read_table(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed.bim', names=['chr', 'SNP_name', 'pos1', 'pos2', 'allele1', 'allele2'])
				bim_file_1000 = pd.read_table(pipeline_config['thousand_genomes']['path'][:-4]+'.bim',  names=['chr', 'SNP_name', 'pos1', 'pos2', 'allele1', 'allele2'])
				extract_for_PCA = list(set(list(bim_file_1000['SNP_name'])) & set(list(bim_file['SNP_name'])))
				keep_these_snps_in_1000.write('\n'.join(extract_for_PCA)) # input file for PLINK extraction of snps in 1000 genomes
				keep_these_snps_in_1000.close() # flushes and closes file

				no_suffix = pipeline_config['thousand_genomes']['path'][:-4]
				

				general_plink.run(
					Parameter('--bfile', no_suffix),
					Parameter('--extract', keep_these_snps_in_1000.name),
					Parameter('--make-bed'),
					Parameter('--out', no_suffix + '_extracted_from_passing_project_SNPs')
					)

				general_plink.run(
					Parameter('--bfile', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed'),
					Parameter('--extract', keep_these_snps_in_1000.name),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed' +'_SNPs_for_use_with_1000_merge')
					)


				general_plink.run(
					Parameter('--bfile', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed' +'_SNPs_for_use_with_1000_merge'),
					Parameter('--bmerge', no_suffix + '_extracted_from_passing_project_SNPs.bed', no_suffix + '_extracted_from_passing_project_SNPs.bim', no_suffix + '_extracted_from_passing_project_SNPs.fam'),
					Parameter('--allow-no-sex'),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen')
					)

				step_order.pop(0)


			elif step_order[0] == 'KING':
				
				print "running KING step"
				phenoFile_thous_Genesis = open(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen_phenoGENESIS.txt', 'w')
				phenoFile_Genesis = open(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_phenoGENESIS.txt', 'w')

				# run KING and output file as -b prefix name ending in .kin, .kin0 with 1000 genomes
				general_king.run(
					Parameter('-b', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen.bed'),
					Parameter('--prefix', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen')
					)
				
				# run KING and output file as -b prefix name ending in .kin, .kin0 WITHOUT 1000 genomes
				general_king.run(
					Parameter('-b', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed.bed'),
					Parameter('--prefix', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed')
					)

				# generate phenotype table for input into GENESIS setup analysis pipeline with 1000 genomes
				pheno_1000_Genesis = pd.read_table(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen.fam', delim_whitespace=True, names = ['FID,', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'], dtype=str)
				pheno_1000_Genesis.loc[(pheno_1000_Genesis['AFF']!='1') & (pheno_1000_Genesis['AFF']!='2'), 'AFF']='NA' # required to put NA for association analysis downstream or miscalculates and fails
				pheno_1000_Genesis[['IID', 'AFF']].to_csv(phenoFile_thous_Genesis.name, sep='\t', index=False, header=False) # format it FID <tab> IID <new line>
				
				# generate phenotype table for input into GENESIS setup analysis pipeline WITHOUT 1000 genomes
				pheno_Genesis = pd.read_table(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed.fam', delim_whitespace=True, names = ['FID,', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'], dtype=str)
				pheno_Genesis.loc[(pheno_Genesis['AFF']!='1') & (pheno_Genesis['AFF']!='2'), 'AFF']='NA' # required to put NA for association analysis downstream or miscalculates and fails
				pheno_Genesis[['IID', 'AFF']].to_csv(phenoFile_Genesis.name, sep='\t', index=False, header=False) # format it FID <tab> IID <new line>

				step_order.pop(0)


			elif step_order[0] == 'PCA':
				
				print "running PCA step"
				#Popen should launch jobs in parallel
				processes = []
				processes.append(subprocess.Popen(['Rscript', 'GENESIS_setup_ANALYSIS_PIPELINE.R', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen', phenoFile_thous_Genesis.name, pipeline_config['R_libraries']['path']]))
				processes.append(subprocess.Popen(['Rscript', 'GENESIS_setup_ANALYSIS_PIPELINE.R', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed', phenoFile_Genesis.name, pipeline_config['R_libraries']['path']]))

				for job in processes:
					job.wait() # wait for all parallel jobs to finish before proceeding to next step

				step_order.pop(0)


			# this is the step at which analysis will be restarted so as to add PCs from
			# previous step
			elif step_order[0] == 'GENanalysis':
				final_results_merged = open(outdir +'/final_results_merged.txt', 'a+')
				
				# calculates the number of chunked files that will be produced and what the names of those files will be

				num_snps = sum(1 for line in open(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed.bim'))
				num_files = int(math.ceil((float(num_snps)/float(100000))))

				# stores name of all chunked files
				total_files = []
				for split_range in range(0, num_files):
					if split_range < 10: # required because linux split uses 00, 01, 02, 03, format therefore when using python iterator a 0 must be appended to front
						total_files.append(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removedsplit0'+str(split_range)+'.results.txt')
					else:
						total_files.append(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removedsplit'+str(split_range)+'.results.txt')


				# run a shell script which will submit slurm script
				try:
					subprocess.call(['./export_var_slurm_streamlined.sh', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed',pipeline_config['R_libraries']['path'], pipeline_args['usePCs']])
				# concatenate all results together with only one line of header
				except:
					sys.exit("Problem submitting slurm script, system exiting...")


				## FILE CHECK HERE -- DO NOT PROCEED UNTIL ALL FILES CREATED

				while True:
					if all(os.path.isfile(chunk_file) for chunk_file in total_files):
						break;
					else:
						print "waiting for split files..."
						time.sleep(60) # wait 1 minute before checking if files are done

				subprocess.call(['head', '-n', '1', outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removedsplit00.results.txt'], stdout=final_results_merged)
				final_results_merged.flush()
				for filename in total_files:
					subprocess.call(['tail', '-n', '+2', '-q', filename], stdout=final_results_merged)
					final_results_merged.flush()
	
				case_controls = pd.read_table(outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed.fam', delim_whitespace=True, names = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'])
				total_case_controls = collections.Counter(list(case_controls['AFF']))
				try:
					controls = str(total_case_controls[1])
					cases = str(total_case_controls[2])
				except KeyError:
					controls = 'NA'
					cases = 'NA'
				# creates Manhattan and qqplots of data
				subprocess.call(['Rscript', 'genesis_clean_qqman_ANALYSIS_PIPELINE.R', final_results_merged.name, outdir + '/merged_group_files/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed.bim', pipeline_config['R_libraries']['path'], pipeline_args['projectName'], outdir + '/merged_group_files/', cases, controls])
				step_order.pop(0)
		
		print "writing results to PDF"
		paramsThresh = summary_stats.parameters_and_thresholds(params=pipeline_args)


		
		# output PDFs -- need to make this compatible with --reanalyze
		# put these under each of the steps
		paramsThresh.output(outdir + '/' + pipeline_args['projectName'] + '_parameters_and_thresholds.pdf', 'F') # output results to PDF
