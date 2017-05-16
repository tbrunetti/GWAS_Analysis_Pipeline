from chunkypipes.components import *
import datetime
import sys
import os

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
		parser.add_argument('--projectName', default=str(datetime.datetime.now()), type=str, help='[default=date time stamp] Name of project')
		parser.add_argument('--startStep', default='hwe', type=str, help='The part of the pipeline you would like to start with')
		parser.add_argument('--endStep', default=None, type=str, help='Point of the pipeline where you would like to stop analysis, if none specified, stops after start step is completed')
		parser.add_argument('--hweThresh', default=1e-6, help='Filters out SNPs that are smaller than this threshold due to liklihood of genotyping error')
		parser.add_argument('--LDmethod', default='indep', type=str, help='[default=indep, options:indep, indep-pairwise, indep-pairphase] Method to calculate linkage disequilibrium')
		parser.add_argument('--VIF', default=2, type=int, help='[default=2] variant inflation factor for indep method LD pruning')
		parser.add_argument('--rsq', default=0.50, type=float, help='[default=0.50] r squared threshold for indep-pairwise or indep-pairphase LD pruning method')
		parser.add_argument('--windowSize', default=50, type=int, help='[default=50] the window size in kb for LD analysis')
		parser.add_argument('--stepSize', default=5, type=int, help='[default=5] variant count to shift window after each interation')
		parser.add_argument('--maf', default=0.05, type=float, help='[default=0.05], filter remaining LD pruned variants by MAF')
		parser.add_argument('--hetThresh', default=0.10, type=float, help='[default=0.10], filter out samples where inbreeding coefficient is greater than threshold (heterozygosity filtering)')
		parser.add_argument('--reanalyze', action='store_true', help='by adding this flag, it means you are going to pass a dataset through the pipeline that has already been partially/fully analyzed by this pipeline. WARNING! May over write exisiting data!!')
		parser.add_argument('--user', default=None, type=str, help='user name for HPC charge')
		parser.add_argument('--projectCharge', default=None, type=str, help='project name on HPC to charge')

	@staticmethod
	def check_steps(order, start, stop):
		pass

	# checks the type of PLINK file input
	@staticmethod
	def check_plink_format(plinkFile, plink):
		if plinkFile[-4:].lower() == '.ped':
			print "Input .ped, converting to binary"
			convert = subprocess.call(['./'+str(plink), '--file', str(plinkFile[:-4]), '--make-bed', '--out', str(plinkFile[:-4])])
		
		elif plinkFile[-4:].lower() == '.bed':
			print "Input seems to follow bed format"

		else:
			sys.exit("Error!! Input not recognized, please input .ped or .bed PLINK file" )


	# creates files to input into plink for samples to keep
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
		
		return make_plinks


	def run_pipeline(self, pipeline_args, pipeline_config):
		import pandas as pd
		sys.path.append(".")
		import subprocess
		import summary_stats
		import statistics as stats
		from fpdf import FPDF
		import PyPDF2

		
		reduced_plink_name = pipeline_args['inputPLINK'].split('/')[-1][:-4] #only get plink file name not full absolute path and removes suffix
		
		# specifying output location and conflicting project names of files generated	
		try:
			if pipeline_args['reanalyze'] == True:
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




		#step_order = ['hwe', 'LD', 'maf', 'merge', 'het', ibd', '1000_genomes', 'KING', 'PCA'] # order of pipeline if full suite is used
		#step_order = ['hwe', 'LD', 'maf', 'merge', 'het', 'ibd', '1000_genomes', 'KING', 'PCA']
		step_order = ['GENanalysis']
		# initialize PLINK and KING software
		general_plink = Software('plink', pipeline_config['plink']['path'])
		general_king = Software('king', pipeline_config['king']['path'])
		# make sure plink input format is correct
		self.check_plink_format(
			plinkFile = pipeline_args['inputPLINK'],
			plink = pipeline_config['plink']['path']
			)

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
				hwe_passing = {}
				for directories in os.listdir(outdir):
					if os.path.isdir(os.path.join(outdir, directories)):

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
			
				step_order.pop(0)

			
			# LD pruning
			elif step_order[0] == 'LD':
				print "running LD pruning step"
				ld_passing = {}

				if pipeline_args['LDmethod'] == 'indep': # gets lists of variants to keep and to prune using VIP
					for directories in os.listdir(outdir):
						if os.path.isdir(os.path.join(outdir, directories)):
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
						if os.path.isdir(os.path.join(outdir, directories)):
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
					
				step_order.pop(0)

			# filters pruned variants by MAF
			elif step_order[0] == 'maf':
				print "running maf step"
				list_of_merge_maf_greater_thresh = open(outdir + '/' + reduced_plink_name + '_greater_mafs.txt', 'w')
				list_of_merge_maf_less_thresh = open(outdir + '/' + reduced_plink_name + '_small_mafs.txt', 'w')
				
				maf_passing = {}
				for directories in os.listdir(outdir):
					if os.path.isdir(os.path.join(outdir, directories)):
						
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
					Parameter('--out', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged')
					)

				# recalculate freq stats for merged set
				general_plink.run(
					Parameter('--bfile', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged'),
					Parameter('--freq'),
					Parameter('--out', outdir + '/' + reduced_plink_name + '_freq_recalculated_greater_all_merged')
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
					Parameter('--out', outdir + '/' + reduced_plink_name + '_maf_less_thresh_all_ethnic_groups_merged')
					)


				# recalculate freq stats for merged set
				general_plink.run(
					Parameter('--bfile', outdir + '/' + reduced_plink_name + '_maf_less_thresh_all_ethnic_groups_merged'),
					Parameter('--freq'),
					Parameter('--out', outdir + '/' + reduced_plink_name + '_freq_recalculated_less_all_merged')
					)

				step_order.pop(0)


			elif step_order[0] == 'het':
				print "checking heterozygosity"
				general_plink.run(
					Parameter('--bfile', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged'),
					Parameter('--het'),
					Parameter('--out', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged')
					)


				het_dataframe = pd.read_table(outdir + '/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged.het', delim_whitespace=True)
				samples_failing_het, het_pdf = summary_stats.heterozygosity(het_dataframe = het_dataframe, thresh = pipeline_args['hetThresh'], outDir = outdir)
				
				general_plink.run(
					Parameter('--bfile', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_all_ethnic_groups_merged'),
					Parameter('--remove', samples_failing_het),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged')
					)
				
				step_order.pop(0)
				


			# determine the pairwise relationship of all merged samples
			# DUPLICATES MUST BE REMOVED BEFORE USING GENESIS PIPELINE!
			elif step_order[0] == 'ibd':
				# only gets run on the bed file with MAF > specified tresh
				print "running PLINK ibd step"
				general_plink.run(
					Parameter('--bfile', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged'),
					Parameter('--genome'),
					Parameter('--read-freq', outdir + '/' + reduced_plink_name + '_freq_recalculated_greater_all_merged.frq'),
					Parameter('--out', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged')
					)
				
				ibd_results = pd.read_table(outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged.genome', delim_whitespace=True)
				relatedness_stats, remove_samples = summary_stats.relatedness(ibd_dataframe = ibd_results, outDir=outdir)

				# remove samples that are duplicates
				general_plink.run(
					Parameter('--bfile', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged'),
					Parameter('--remove', remove_samples),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed')
					)
				
				step_order.pop(0)
	

			
			elif step_order[0] == '1000_genomes':
				
				no_suffix = pipeline_config['thousand_genomes'][:-4]
				general_plink.run(
					Parameter('--bfile', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed'),
					Parameter('--bmerge', no_suffix + '.bed', no_suffix + '.bim', no_suffix + '.fam'),
					Parameter('--make-bed'),
					Parameter('--out', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen')
					)


				step_order.pop(0)


			elif step_order[0] == 'KING':
				
				print "running KING step"
				phenoFile_Genesis = open(outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen_phenoGENESIS.txt', 'w')
				# run KING and output file as -b prefix name ending in .kin, .kin0
				general_king.run(
					Parameter('--bfile', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen.bed'),
					Parameter('--prefix', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen')
					)

				# generate phenotype table for input into GENESIS analysis  pipeline
				pheno_Genesis = pd.read_table(outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen.fam', delim_whitespace=True, names = ['FID,', 'IID', 'PAT', 'MAT', 'SEX', 'AFF'])
				pheno_Genesis[['IID', 'AFF']].to_csv(phenoFile_Genesis.name, sep='\t', index=False, header=False) # format it FID <tab> IID <new line>

				step_order.pop(0)


			elif step_order[0] == 'PCA':
				#TO DO merge LD pruned 1000 genomes
				print "running PCA step"
				subprocess.call(['Rscript', 'GENESIS_setup_ANALYSIS_PIPELINE.R', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen', phenoFile_Genesis.name])
				step_order.pop(0)


			# this is the step at which analysis will be restarted so as to add PCs from
			# previous step
			elif step_order[0] == 'GENanalysis':
				final_results_merge = open(outdir +'/final_results_merged.txt', 'a+')
				# run a shell script which will submit slurm script
				subprocess.call(['./export_var_slurm_streamlined.sh', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen'])
				# concatenate all results together with only one line of header
				subprocess.call(['head', '-n', '1', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGensplit00.results.txt'], stdout=final_results_merged.name)
				subprocess.call(['tail', '-n', '+2', '-q', outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGensplit0*.results.txt'], stdout=final_results_merged.name)
				
				# creates Manhattan and qqplots of data
				subprocess.call(['Rscript', 'genesis_clean_qqman_ANALYSIS_PIPELINE.R', final_results_merged.name, outdir + '/' + reduced_plink_name + '_maf_greater_thresh_hetFiltered_all_ethnic_groups_merged_dups_removed_thousGen.bim'])

		
		print "writing results to PDF"
		paramsThresh = summary_stats.parameters_and_thresholds(params=pipeline_args)


		# output PDFsS
		relatedness_stats.output(outdir + '/' + pipeline_args['projectName'] + '_relatedness.pdf', 'F')
		paramsThresh.output(outdir + '/' + pipeline_args['projectName'] + '_parameters_and_thresholds.pdf', 'F')
		hwe_stats.output(outdir + '/' + pipeline_args['projectName'] + '_hweStats.pdf', 'F')
		ld_stats.output(outdir + '/' + pipeline_args['projectName'] + '_ldStats.pdf', 'F')
		maf_stats.output(outdir + '/' + pipeline_args['projectName'] + '_mafStats.pdf', 'F')
		het_pdf.output(outdir + '/' + pipeline_args['projectName'] + '_hetStats.pdf', 'F')
		