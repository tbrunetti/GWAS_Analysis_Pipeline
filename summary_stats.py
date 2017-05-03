import pandas as pd
import statistics as stats
import matplotlib.pyplot as plt
import numpy as np

from fpdf import FPDF

def parameters_and_thresholds(params):
	
	pdf = FPDF()
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Parameters and Thresholds", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', '', 16)
	for key in params:
		if key not in ['inputPLINK', 'phenoFile', 'outDir', 'projectName', 'config']:
			pdf.multi_cell(0, 8, str(key)+':     '+str(params[key]), 0, 1, 'L')

	return pdf



def hwe(dictHWE, thresh, outDir):
	pdf = FPDF() # create new PDF
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Hardy-Weinberg Equilibrium", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_x(20)
	pdf.set_font('Arial', '', 12)
	pdf.multi_cell(0, 8, 'Hardy-Weinberg equilibrium is only used to remove SNPs with extreme p-values are that are likely \
		to occur due to sequencing, genotyping, or study-design errors.  This calculation is sensitive to different ethinic groups \
		and races.  Therefore, it is independently calculated for each ethnic group.  The current p-value threshold that was used to determine \
		whether a SNP was removed was  ' + str(thresh) + '.  This calculation will only consider founders nonfounders are ignored.' , 0, 1, 'J')
	pdf.set_font('Arial', 'B', 16)
	pdf.set_fill_color(200)
	
	# iterate through all ethnic groups for HWE stats
	for key, value in dictHWE.iteritems():
		pdf.multi_cell(0, 8, str(key), 1, 'L', True)
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs analyzed:  ' +  str(value[0]), 1, 1, 'L')
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs Passing:  ' +  str(value[1]) + ' (' + str("%.2f" % round((float(value[1])/float(value[0]))*100, 2)) + '%)', 1, 1, 'L')
		pdf.multi_cell(0, 8, '\n\n', 0, 1, 'J')

		# NOTE hweFile is before filtering by HWE threshold and prior to removal of SNPs failing threshold
		# these plot the observed vs expected from pre-filter HWE and the associated p-values
		# red fail threhold
		hweFile_dataframe = pd.read_table(value[2], delim_whitespace=True)	
		for phenotypes in list(set(list(hweFile_dataframe['TEST']))):
			pheno_subset = hweFile_dataframe.loc[hweFile_dataframe['TEST'] == phenotypes]
			colors = np.where(pheno_subset.P < thresh, 'r', 'k')
			plt.scatter(pheno_subset['E(HET)'], pheno_subset['O(HET)'], c=colors, s=8)
			plt.xlabel('expected(het)', fontsize=14)
			plt.ylabel('observed(het)', fontsize=14)
			plt.suptitle(phenotypes + ':  observed(het) vs expected(het) of HWE', fontsize=14)
			plt.tight_layout(pad=2, w_pad=2, h_pad=2)
			plt.savefig(outDir+'/'+'hwe_'+str(key)+'_'+str(phenotypes)+'.png', bbox_inches='tight')
			plt.close()
	
	return pdf

def pruning(dictLD):
	pdf = FPDF() # create new PDF
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "LD Pruning", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', 'B', 16)
	pdf.set_fill_color(200)
	
	# iterate through all ethnic groups for LD pruning stats
	for key, value in dictLD.iteritems():
		pdf.multi_cell(0, 8, str(key), 1, 'L', True)
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs analyzed:  ' +  str(value[0]), 1, 1, 'L')
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs Passing:  ' +  str(value[1]) + ' (' + str("%.2f" % round((float(value[1])/float(value[0]))*100, 2)) + '%)', 1, 1, 'L')
		pdf.multi_cell(0, 8, '\n\n', 0, 1, 'J')


	return pdf

def minor_allele_freq(dictMAF, thresh):
	pdf = FPDF() # create new PDF
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Minor Allele Frequency", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)
	pdf.set_font('Arial', 'B', 16)
	pdf.set_fill_color(200)

	# iterate through all ethnic groups for LD pruning stats
	for key, value in dictMAF.iteritems():
		pdf.multi_cell(0, 8, str(key), 1, 'L', True)
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs analyzed:  ' +  str(value[0]), 1, 1, 'L')
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs MAF >= ' + str(thresh) + ': ' +  str(value[1]) + ' (' + str("%.2f" % round((float(value[1])/float(value[0]))*100, 2)) + '%)', 1, 1, 'L')
		pdf.set_x(30)
		pdf.multi_cell(0, 8, 'Total Number of SNPs MAF <= ' + str(thresh) + ': ' +  str(value[2]) + ' (' + str("%.2f" % round((float(value[2])/float(value[0]))*100, 2)) + '%)', 1, 1, 'L')
		pdf.multi_cell(0, 8, '\n\n', 0, 1, 'J')


	return pdf


def relatedness(ibd_dataframe, outDir):
	
	dups_text = open(outDir + '/' + 'duplicate_pairs.txt', 'w')
	pdf = FPDF() # create new PDF
	pdf.add_page()
	pdf.set_margins(20, 10, 20)
	pdf.set_font('Arial', 'B', 24)
	pdf.set_x(20)
	pdf.multi_cell(0, 30, "Relatedness", 0, 1, 'L')
	pdf.line(20, 32, 190, 32)


	pdf.set_font('Arial', 'B', 16)
	pdf.set_fill_color(200)
	pdf.multi_cell(0, 10, 'Total Number of Sample Pairs Analyzed:  ' +  str(len(ibd_dataframe.index)), 1, 'L', True)

	
	duplicates = ibd_dataframe.loc[ibd_dataframe['Z2'] > 0.97]
	parent_offspring = ibd_dataframe.loc[(ibd_dataframe['Z1'] > 0.97) & (ibd_dataframe['Z0'] < 0.05) & (ibd_dataframe['Z2'] < 0.05)]
	full_sibs = ibd_dataframe.loc[(ibd_dataframe['Z0'] < 0.40) & (ibd_dataframe['Z2'] > 0.16)]
	half_sibs = ibd_dataframe.loc[(ibd_dataframe['Z0'] < 0.60) & (ibd_dataframe['Z1'] < 0.58) & (ibd_dataframe['Z2'] < 0.05)]
	cousins = ibd_dataframe.loc[(ibd_dataframe['Z0'] > 0.60) & (ibd_dataframe['Z1'] < 0.40) & (ibd_dataframe['Z2'] < 0.02)]
	unrelated = ibd_dataframe.loc[ibd_dataframe['Z0'] > 0.78]
	
	pdf.set_font('Arial', '', 16)
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of duplicate pairs:  '+str(len(duplicates.index)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of parent-offspring pairs:  '+str(len(parent_offspring.index)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of full siblings pairs:  '+str(len(full_sibs.index)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of half sibling pairs:  '+str(len(half_sibs.index)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of cousin pairs:  '+str(len(cousins.index)), 1, 1, 'L')
	pdf.set_x(30)
	pdf.multi_cell(0, 10, '# of unrelated pairs:  '+str(len(unrelated.index)), 1, 1, 'L')


	duplicates.to_csv(dups_text.name, sep='\t', index=False)

	return pdf

