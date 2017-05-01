import pandas as pandas
import statistics as stats
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
		if key not in ['inputPLINK', 'phenoFile', 'outDir', 'projectName']:
			pdf.multi_cell(0, 8, str(key)+':     '+str(params[key]), 0, 1, 'L')

	return pdf



def relatedness(ibd_dataframe):

	#relatedness_info = open('related_samples.txt', 'w')
	
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


	return pdf

