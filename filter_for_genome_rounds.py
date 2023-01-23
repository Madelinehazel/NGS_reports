import pandas as pd
import argparse
from group_variants import variants


def main(report, file):
    """
    From a C4R variant csv, creates an excel document with the following tabs:
    - One tab containing ALL variants
    - One tab containing variants filtered by gnomAD homozygous count, C4R counts, and quality
    - One tab with above filters containing variants only in OMIM genes (also includes ClinVar pathogenic variants >1% AF)
    """

    # add summary for each variant describing pathogenicity predictions and gnomad frequency
    summaries = []
    for index, variant in report.iterrows():
        try:
            summary = variants.summary_field(
                variant["Cadd_score"],
                variant["Sift_score"],
                variant["Polyphen_score"],
                variant["Vest3_score"],
                variant["Revel_score"],
                variant["Gnomad_ac"],
                variant["Gnomad_hom"],
                variant["Refseq_change"],
                variant["Variation"],
                variant["Gene"],
                variant["C4R_WES_counts"],
                variant["Quality"],
                variant["Exac_pli_score"],
            )
        except KeyError:
            summary = variants.summary_field(
                variant["Cadd_score"],
                variant["Sift_score"],
                variant["Polyphen_score"],
                variant["Vest3_score"],
                variant["Revel_score"],
                variant["Gnomad_ac"],
                variant["Gnomad_hom"],
                variant["Refseq_change"],
                variant["Variation"],
                variant["Gene"],
                variant["Frequency_in_C4R"],
                variant["Quality"],
                variant["Exac_pli_score"],
            )
        summaries.append(summary)
    report["Summary"] = summaries
    cols = list(report.columns)
    cols = [cols[-1]] + cols[:-1]
    report = report[cols]

    # first get clinvar path > 1% so can add to OMIM tab
    clinvar_greater_than_1 = report[report["Gnomad_af_popmax"] > 0.01]

    # apply gnomAD, C4R counts, quality
    try:
        report_filter = report[report["Frequency_in_C4R"] < 10]
    except KeyError:
        report_filter = report[report["C4R_WES_counts"] < 10]

    report_filter = report_filter[
        (report_filter["Gnomad_hom"] == 0) & (report_filter["Quality"] >= 300)
    ]

    # get variants in OMIM genes
    omim = report_filter[
        (report_filter["omim_phenotype"] != ".")
        & (
            (report_filter["omim_phenotype"] == report_filter["omim_phenotype"])
            | (report_filter["omim_phenotype"].notnull())
        )
    ]

    omim_clinvar = pd.concat([omim, clinvar_greater_than_1], ignore_index=True)

    # write report
    with pd.ExcelWriter("%s_for_exome_rounds.xlsx" % file) as writer:
        workbook = writer.book
        report.to_excel(writer, sheet_name="all", index=False)
        report_filter.to_excel(writer, sheet_name="rare_high_qual", index=False)
        omim_clinvar.to_excel(writer, sheet_name="rare_high_qual_omim", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter variants for streamlined SickKids WES analysis"
    )
    parser.add_argument("-report", type=str, help="input report tsv")
    args = parser.parse_args()

    file = args.report.replace(".tsv", "")
    report = pd.read_csv(args.report, encoding="ISO-8859-1", sep="\t")
    main(report, file)
