import pandas as pd 
import argparse 
from group_variants import variants

def main(report, file):
    summaries = []
    # add summary for each variant describing pathogenicity predictions and gnomad frequency
    for index, variant in report.iterrows():
        summary = variants.summary_field(variant['Cadd_score'], variant['Sift_score'], variant['Polyphen_score'], variant['Vest3_score'],
                            variant['Revel_score'], variant['Gnomad_ac'], variant['Gnomad_hom'], variant['Refseq_change'],
                            variant['Variation'], variant['Gene'])
        summaries.append(summary)
    report['Summary'] = summaries
    cols = list(report.columns)
    cols = [cols[-1]] + cols[:-1]
    report = report[cols]


    report.to_csv('%s_with_summaries.csv' % file, index=False, encoding="ISO-8859-1")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generates summary column from raw variant report')
    parser.add_argument('-report', type=str, help='input report csv')
    args = parser.parse_args()

    file = args.report.strip('.csv')
    report = pd.read_csv(args.report, encoding="ISO-8859-1")
    main(report, file)

    