import argparse
import pandas as pd
from math import nan
import xlsxwriter

def main(report, file):
    # format missing values 
    report['Cadd_score'].replace('None', nan, inplace=True)
    report['Sift_score'].replace('None', nan, inplace=True)
    report['Polyphen_score'].replace('None', nan, inplace=True)
    report['Gnomad_oe_lof_score'].replace(nan, '', inplace=True)

    # convert column types to float for conditional formatting
    report = report.astype({'Cadd_score': float,'Sift_score': float, 'Polyphen_score': float})

    report['Gnomad_oe_lof_score'].replace('', '.', inplace=True)
    report['Sift_score'].replace(nan, '.', inplace=True)
    report['Polyphen_score'].replace(nan, '.', inplace=True)

    # add blank column for notes
    report['Notes'] = ['']*len(report)

    # place notes and genes columns at beginning of report
    report_cols = [col for col in report.columns if col != 'Notes' and col != 'Gene']
    report = report[['Notes', 'Gene'] + report_cols]

    # made a dictionary of column names to indices
    index = 0 
    col_dict = {}
    for col in report.columns:
        col_dict[col] = index
        index += 1

    # grab positions of columns of interest
    cadd = col_dict['Cadd_score']
    gnomad_af_popmax = col_dict['Gnomad_af_popmax']
    gnomad_hom = col_dict['Gnomad_hom']
    try:    
        c4r = col_dict['C4R_WES_counts']
    except KeyError:
        c4r = col_dict['Frequency_in_C4R']
    lof_score = col_dict['Gnomad_oe_lof_score']
    pli_score = col_dict['Exac_pli_score']
    sift = col_dict['Sift_score']
    polyphen = col_dict['Polyphen_score']
    quality = col_dict['Quality']

    # create excel writer
    writer = pd.ExcelWriter(f'{file}.xlsx', engine='xlsxwriter')

    # convert report dataframe to an XlsxWriter object
    report.to_excel(writer, sheet_name='Variants', index=False)


    workbook  = writer.book
    worksheet = writer.sheets['Variants']

    # add a separate blank tab for HPO terms
    hpo = workbook.add_worksheet('HPO')

    # light red fill with dark red text
    format1 = workbook.add_format({'bg_color': '#FFC7CE',
                                'font_color': '#9C0006'})

    # green fill with dark green text
    format2 = workbook.add_format({'bg_color': '#C6EFCE',
                                'font_color': '#006100'})

    # yellow fill with dark yellow text
    format3 = workbook.add_format({'bg_color': '#fad97f',
                                'font_color': '#d1a52c'})

    # add a default format for blanks
    format4 = workbook.add_format()

    # apply a conditional format to the required cell range.
    # this will prevent blanks in columns with conditional formatting for values less than a threshold from being highlighted
    worksheet.conditional_format('A1:XFD1048576',  {'type': 'blanks',
                                            'stop_if_true': True,
                                            'format': format4})
    # CADD
    worksheet.conditional_format(0, cadd, 1048575, cadd,{'type': 'cell', 'criteria': '>=', 'value': 15, 'format': format2})
    # C4R counts
    worksheet.conditional_format(0, c4r, 1048575, c4r,{'type': 'cell', 'criteria': '>', 'value': 10, 'format': format1})
    # gnomAD AF popmax
    worksheet.conditional_format(0, gnomad_af_popmax, 1048575, gnomad_af_popmax,{'type': 'cell', 'criteria': '>', 'value': 0.005, 'format': format1})
    # gnomAD hom
    worksheet.conditional_format(0, gnomad_hom, 1048575, gnomad_hom,{'type': 'cell', 'criteria': '>', 'value': 2, 'format': format1})
    #gnomAD oe lof
    worksheet.conditional_format(0, lof_score, 1048575, lof_score,{'type': 'cell', 'criteria': '<', 'value': 0.35, 'format': format3})
    #plI
    worksheet.conditional_format(0, pli_score, 1048575, pli_score,{'type': 'cell', 'criteria': '>', 'value': 0.9, 'format': format3})
    #SIFT
    worksheet.conditional_format(0, sift, 1048575, sift,{'type': 'cell', 'criteria': '<', 'value': 0.05, 'format': format3})
    #polyphen
    worksheet.conditional_format(0, polyphen, 1048575, polyphen,{'type': 'cell', 'criteria': '>', 'value': 0.9, 'format': format3})
    # quality
    worksheet.conditional_format(0, quality, 1048575, quality, {'type': 'cell', 'criteria': '<', 'value': 500, 'format': format1})


    # Close the Pandas Excel writer and output the Excel file.
    writer.save()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert csv to xlsx and add conditional formatting')
    parser.add_argument('-report', type=str, help='input report csv')
    args = parser.parse_args()

    file = args.report.replace('.csv', '')
    report = pd.read_csv(args.report, encoding="ISO-8859-1")
    main(report, file)