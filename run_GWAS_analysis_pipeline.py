from chunkypipes.components import *
import datetime
import sys
import os

class Pipeline(BasePipeline):
	def dependencies(self):
		return ['pandas']

	def description(self):
		return 'Pipeline made for analyzing GWAS data after QC cleanup'

	def configure(self):
		return {
			'plink':{
				'path':'Full path to PLINK executable'
			}
		}

	def add_pipeline_args(self, parser):
		# should other input options be available??
		parser.add_argument('-inputPLINK', required=True, type=str, help='Full path to PLINK file ending in .BED or .PED')
		parser.add_argument('-phenoFile', required=True, type=str, help='Full path to phenotype file, see argument readme for more details on format')
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
	def ethnic_plinks_lists(phenotype, plinkFileName, outDir):
		import pandas as pd

		phenotype_table = pd.read_table(phenotype)
		race_subsets = list(set(list(phenotype_table['Race']))) # gets all racial groups listed in phenotype file
		race_subsets_cleaned = [i for i in race_subsets if str(i)!='nan'] # removes samples that do not have a race listed 
		make_plinks = {} # stores file name of samples to keep for each group
		
		for ethnic_group in race_subsets_cleaned:
			os.mkdir(outDir + '/' + '_'.join(ethnic_group.split()))
			keep_file = open(outDir + '/' + '_'.join(ethnic_group.split()) + '/' + plinkFileName +'_'+str('_'.join(ethnic_group.split())) + '_keepIDs.txt', 'w')
			subset_by_race = phenotype_table.loc[phenotype_table['Race'] == ethnic_group]
			subset_by_race[['FID', 'IID']].to_csv(keep_file.name, sep='\t', index=False, header=False) # format it FID <tab> IID <new line>
			make_plinks['_'.join(ethnic_group.split())]=keep_file.name
			keep_file.flush()
			keep_file.close()
		
		return make_plinks


	def run_pipeline(self, pipeline_args, pipeline_config):
		import pandas as pd

		reduced_plink_name = pipeline_args['inputPLINK'].split('/')[-1][:-4] #only get plink file name not full absolute path and removes suffix
		
		# specifying output location and conflicting project names of files generated	
		try:
			os.stat(pipeline_args['outDir']+'/'+pipeline_args['projectName'])
			sys.exit("project already exists!!")
		except:
			print "Making new directory called "+str(pipeline_args['projectName']) + ' located in ' + str(pipeline_args['outDir'])
			outdir = pipeline_args['outDir']+'/'+pipeline_args['projectName'] # new output directory
			os.mkdir(outdir)

	
		#step_order = ['hwe', 'LD', 'maf', het', 'ibs' 'GENESIS'] # order of pipeline if full suite is used
		step_order = ['hwe', 'LD', 'maf']
		# initialize PLINK software
		general_plink=Software('plink', pipeline_config['plink']['path'])
		
		# make sure plink input format is correct
		self.check_plink_format(
			plinkFile = pipeline_args['inputPLINK'],
			plink = pipeline_config['plink']['path']
			)

		keep_files = self.ethnic_plinks_lists(
			phenotype = pipeline_args['phenoFile'],
			plinkFileName = reduced_plink_name,
			outDir = outdir
			)


		# will make separate plink files for each ethnic group
		for key, value in keep_files.iteritems():
			general_plink.run(
				Parameter('--bfile', pipeline_args['inputPLINK'][:-4]),
				Parameter('--keep', value),
				Parameter('--make-bed'),
				Parameter('--out', outdir + '/' + str(key) + '/' + reduced_plink_name + '_' + str(key))
				)


		while len(step_order) != 0:
			
			# hardy-weinberg equilibrium filtering
			if step_order[0] == 'hwe':
				for directories in os.listdir(outdir):
					general_plink.run(
						Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories),
						Parameter('--hardy'),
						Parameter('--hwe', pipeline_args['hweThresh']),
						Parameter('midp'),
						Parameter('--make-bed'),
						Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_hweFiltered')
						)
				step_order.pop(0)

			
			# LD pruning
			elif step_order[0] == 'LD':
				if pipeline_args['LDmethod'] == 'indep': # gets lists of variants to keep and to prune using VIP
					for directories in os.listdir(outdir):
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
							Parameter('--'+pipeline_args['LDmethod'], str(pipeline_args['windowSize'])+'kb', str(pipeline_args['stepSize']), str(pipeline_args['VIF'])),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories)
							)
						
						# creates new PLINK files with excluded variants removed
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
							Parameter('--exclude', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.prune.out'),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned')
							)
						
				else:
					for directories in os.listdir(outdir): # get lists of variants to keep and to prune using rsq
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
							Parameter(pipeline_args['LDmethod'], str(pipeline_args['windowSize'])+'kb', str(pipeline_args['stepSize']), str(pipeline_args['rsq'])),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories)
							)
						
						# creates new PLINK files with excluded variants removed
						general_plink.run(
							Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name + '_' + directories + '_hweFiltered'),
							Parameter('--exclude', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '.prune.out'),
							Parameter('--make-bed'),
							Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned')
							)
				
				
				step_order.pop(0)

			# filters pruned variants by MAF
			elif step_order[0] == 'maf':
				for directories in os.listdir(outdir):
					# filters out variants below set maf threshold
					general_plink.run(
						Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned'),
						Parameter('--maf', str(pipeline_args['maf'])),
						Parameter('--make-bed'),
						Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_greater_than_'+str(pipeline_args['maf'])+'_maf')
						)

					# filters out variants set maf threshold
					general_plink.run(
						Parameter('--bfile', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_LDpruned'),
						Parameter('--max-maf', str(pipeline_args['maf'])),
						Parameter('--make-bed'),
						Parameter('--out', outdir + '/' + directories + '/' + reduced_plink_name+ '_' + directories + '_less_than_or_equal_'+str(pipeline_args['maf'])+'_maf')
						)

				step_order.pop(0)



			print step_order

'''

				elif step_order[0] == 'het':
					step_order.pop(0)

				elif step_order[0] == 'ibs':
					step_order.pop(0)

				elif step_order[0] == 'GENESIS':
					step_order.pop(0)
'''