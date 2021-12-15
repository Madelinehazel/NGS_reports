import pandas as pd
import argparse
from group_variants import variants


def main(
    report, proband_id, report_type, family_id, file, maternal_id, paternal_id, panel
):
    summaries = []

    proband_id = variants.parse_id(args.family_id, proband_id)
    print(proband_id)

    if args.maternal_id and args.paternal_id:
        maternal_id = variants.parse_id(args.family_id, args.maternal_id)
        paternal_id = variants.parse_id(args.family_id, args.paternal_id)

    # add summary for each variant describing pathogenicity predictions and gnomad frequency
    for index, variant in report.iterrows():
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
        )
        summaries.append(summary)
    report["Summary"] = summaries
    cols = list(report.columns)
    cols = [cols[-1]] + cols[:-1]
    report = report[cols]

    # filter out variants with quality < 300
    report = report[report["Quality"] >= 300]

    # if it's a singleton, apply gnomAD, C4R counts, and impact filters to all variants
    if report_type == "singleton":
        report = report[report["Frequency_in_C4R"] < 100]
        report = report[report["Gnomad_hom"] == 0]
        report = report[
            report.apply(
                lambda x: variants.filter_impact(x.Variation, x.Clinvar), axis=1
            )
        ]

    # get variants by inheritance pattern
    AR = variants.autosomal_recessive(
        report, proband_id, maternal_id, paternal_id, report_type
    )
    comp_het = variants.compound_het(
        report, proband_id, maternal_id, paternal_id, report_type
    )
    hemi = variants.hemizygous(
        report, proband_id, maternal_id, paternal_id, report_type
    )
    omim = variants.dominant_OMIM(report, proband_id)
    if report_type == "trio":
        de_novo = variants.denovo(report, proband_id, maternal_id, paternal_id)
    else:
        dominant_nonOMIM = variants.dominant_nonOMIM(report, proband_id)
        dominant_nonOMIM = dominant_nonOMIM[
            dominant_nonOMIM["Exac_pli_score"].apply(variants.filter_pli)
        ]
        dominant_nonOMIM = dominant_nonOMIM[
            dominant_nonOMIM["Variation"].apply(variants.filter_lof)
        ]
    # if variants are annotated with a panel, add a tab containing all variants falling in panel
    if panel:
        panel_variants = variants.panel(report, proband_id)
        panel_variants = panel_variants[
            panel_variants["Cadd_score"].apply(variants.filter_cadd)
        ]

    # write report
    with pd.ExcelWriter("%s_formatted.xlsx" % file) as writer:
        workbook = writer.book
        if report_type == "trio":
            de_novo.to_excel(writer, sheet_name="De_novo", index=False)
        else:
            summary = pd.read_csv(
                "/hpf/largeprojects/ccmbio/mcouse/tools/NGS_reports/summary_pages/singleton.csv"
            )
            summary.to_excel(writer, sheet_name="Summary", index=False)
            dominant_nonOMIM.to_excel(
                writer, sheet_name="Dominant_nonOMIM", index=False
            )
        AR.to_excel(writer, sheet_name="Autosomal_recessive", index=False)
        comp_het.to_excel(writer, sheet_name="Compound_heterozygous", index=False)
        hemi.to_excel(writer, sheet_name="Hemizygous", index=False)
        omim.to_excel(writer, sheet_name="Dominant_OMIM", index=False)
        if panel:
            panel_variants.to_excel(writer, sheet_name="Panels", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates excel report from wes csv with variants grouped by inheritance pattern"
    )
    parser.add_argument("-report", type=str, help="input report csv")
    parser.add_argument("-report_type", type=str, help="singleton or trio")
    parser.add_argument("-family_id", type=str, help="family id")
    parser.add_argument("-proband_id", type=str, help="proband sample id")
    parser.add_argument(
        "-maternal_id", type=str, help="maternal sample id", default=None
    )
    parser.add_argument(
        "-paternal_id", type=str, help="paternal sample id", default=None
    )
    parser.add_argument(
        "-panel",
        type=bool,
        help="True if annotated with gene panel (e.g. immunopanel), default False",
        default=False,
    )
    args = parser.parse_args()

    file = args.report.strip(".csv")
    report = pd.read_csv(args.report, encoding="latin1")

    main(
        report,
        args.proband_id,
        args.report_type,
        args.family_id,
        file,
        args.maternal_id,
        args.paternal_id,
        args.panel,
    )
