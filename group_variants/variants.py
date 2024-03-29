import pandas as pd
import argparse


def get_zygosity(sample_id):
    zygosity = "Zygosity.%s" % sample_id
    return zygosity


def get_burden(sample_id):
    burden = "Burden.%s" % sample_id
    return burden


def autosomal_recessive(variants, proband, mother, father, report_type):
    if report_type == "singleton":
        variants_filtered = variants[
            (variants[get_zygosity(proband)] == "Hom")
            & (variants["Position"].str.contains("X") == False)
        ]
        variants_filtered = variants_filtered.sort_values(["omim_phenotype", "Gene"])
    elif report_type == "trio":
        variants_filtered = variants[
            (variants[get_zygosity(proband)] == "Hom")
            & (variants[get_zygosity(mother)] == "Het")
            & (variants[get_zygosity(father)] == "Het")
        ]
    return variants_filtered


def hemizygous(variants, proband, mother, father, report_type):
    if report_type == "singleton":
        variants_filtered = variants[
            (variants[get_zygosity(proband)] == "Hom")
            & (variants["Position"].str.contains("X") == True)
        ]
    elif report_type == "trio":
        variants_filtered = variants[(variants[get_zygosity(proband)] == "Hom")]
        # hemizygous variants inherited from mom
        variants_filtered_mat = variants_filtered[
            (variants_filtered[get_zygosity(mother)] == "Het")
            & (variants_filtered[get_zygosity(father)] == "-")
        ]
        # hemizygous variants inherited from dad
        variants_filtered_pat = variants_filtered[
            (variants_filtered[get_zygosity(mother)] == "-")
            & (variants_filtered[get_zygosity(father)] == "Het")
        ]
        variants_filtered = pd.concat(
            [variants_filtered_mat, variants_filtered_pat], ignore_index=True
        )
    variants_filtered = variants_filtered.sort_values(["omim_phenotype", "Gene"])
    return variants_filtered


def dominant_nonOMIM(variants, proband):
    variants_filtered = variants[(variants[get_zygosity(proband)] == "Het")]
    variants_filtered = variants_filtered[
        (variants_filtered["omim_phenotype"] != variants_filtered["omim_phenotype"])
        | (variants_filtered["omim_phenotype"] == ".")
    ]
    variants_filtered = variants_filtered.sort_values(["Gnomad_ac", "Gene"])
    return variants_filtered


def denovo(variants, proband, mother, father):
    variants_filtered = variants[
        (variants[get_zygosity(proband)] == "Het")
        | (variants[get_zygosity(proband)] == "Hom")
    ]
    variants_filtered = variants_filtered[
        (variants_filtered[get_zygosity(mother)] == "-")
        & (variants_filtered[get_zygosity(father)] == "-")
    ]
    variants_filtered = variants_filtered.sort_values(["omim_phenotype", "Gene"])
    return variants_filtered


def compound_het(variants, proband, mother, father, report_type):
    # this function captures quite a bit of junky misaligned variants where
    # the burden is high in all samples.
    variants_burden = variants[(variants[get_burden(proband)] >= 2)]
    variants_burden_het = variants_burden[
        (variants_burden[get_zygosity(proband)] == "Het")
    ]

    if report_type == "trio":
        variants_filtered_mat = variants_burden[
            (variants_burden[get_zygosity(mother)] == "Het")
            & (variants_burden[get_zygosity(father)] == "-")
        ]
        variants_filtered_pat = variants_burden[
            (variants_burden[get_zygosity(father)] == "Het")
            & (variants_burden[get_zygosity(mother)] == "-")
        ]

        # now check to see if gene has het variants in both mother and father
        denovo_genes = []
        index_list = []
        for index, variant in variants_burden_het.iterrows():
            gene = variant["Gene"]
            if (
                gene in variants_filtered_mat["Gene"].tolist()
                and gene in variants_filtered_pat["Gene"].tolist()
            ):
                index_list.append(index)
            # include de novo variants_burden
            elif (
                variant[get_zygosity(proband)] == "Het"
                and variant[get_zygosity(mother)] == "-"
                and variant[get_zygosity(father)] == "-"
            ):
                index_list.append(index)
                denovo_genes.append(variant["Gene"])
            elif gene in denovo_genes:
                index_list.append(index)
            else:
                pass
        variants_burden_het = variants_burden.loc[index_list, :]

    variants_burden_het = variants_burden_het.sort_values(["omim_phenotype", "Gene"])
    return variants_burden_het


def dominant_OMIM(variants, proband):
    variants_filtered = variants[
        (variants["omim_phenotype"] != ".")
        & (
            (variants["omim_phenotype"] == variants["omim_phenotype"])
            | (variants["omim_phenotype"].notnull())
        )
        & (variants[get_zygosity(proband)] == "Het")
    ]

    variants_filtered = variants_filtered.sort_values(["omim_inheritance", "Gnomad_ac"])
    return variants_filtered


def panel(variants, proband):
    variants_filtered = variants[
        variants[get_zygosity(proband)].isin(["Het", "Hom"])
        & (variants["Panels"] == variants["Panels"])
    ]

    variants_filtered = variants_filtered.sort_values(["omim_inheritance", "Gnomad_ac"])
    return variants_filtered


def summary_field(
    cadd,
    sift,
    polyphen,
    vest3,
    revel,
    gnomad_ac,
    gnomad_hom,
    variant_transcript,
    variant_type,
    gene,
    c4r,
    qual,
    pli
):
    num_tools = 0
    pathogenic_count = 0
    # cadd
    if cadd == "None" or cadd == ".":
        pass
    elif float(cadd) >= 15:
        num_tools += 1

        pathogenic_count += 1
    else:
        num_tools += 1
    # sift
    if sift == "None" or sift == ".":
        pass
    elif float(sift) < 0.05:
        num_tools += 1
        pathogenic_count += 1
    else:
        num_tools += 1
    # polyphen
    if polyphen == "None" or polyphen == ".":
        pass
    elif float(polyphen) > 0.446:
        num_tools += 1
        pathogenic_count += 1
    else:
        num_tools += 1
    # vest3
    if vest3 == "None" or vest3 == ".":
        pass
    elif float(vest3) > 0.5:
        num_tools += 1
        pathogenic_count += 1
    else:
        num_tools += 1
    # revel
    if revel == "None" or revel == ".":
        pass
    elif float(revel) > 0.5:
        num_tools += 1
        pathogenic_count += 1
    else:
        num_tools += 1

    try:
        qual = round(float(qual))
    except:
        pass

    summary = "CADD = {}; {}/{} tools predict an impact. {} alleles and {} homozygote(s) in gnomAD. Seen {} time(s) in C4R. Quality: {}. {} is a {} variant in {}. Gene plI: {}".format(
        cadd,
        pathogenic_count,
        num_tools,
        gnomad_ac,
        gnomad_hom,
        c4r,
        qual,
        variant_transcript,
        variant_type,
        gene,
        pli
    )
    return summary


def filter_cadd(cadd):
    if cadd == "None":
        return True
    elif float(cadd) >= 15:
        return True
    else:
        return False


def filter_impact(impact, clinvar):
    impacts = [
        "frameshift_variant",
        "inframe_deletion",
        "inframe_insertion",
        "missense_variant",
        "protein_altering_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "start_lost",
        "stop_gained",
        "stop_lost",
    ]
    if clinvar.lower().find("pathogenic") != -1:
        return True
    elif impact in impacts:
        return True
    else:
        return False


def parse_id(family_id, sample_id):
    family_id = family_id.replace("-", "_")
    sample_id = sample_id.replace("-", "_")
    id = family_id + "_" + sample_id
    return id


def filter_pli(pli):
    if pli == ".":
        return False
    elif float(pli) >= 0.95:
        return True
    else:
        return False


def filter_lof(impact):
    impacts = [
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "start_lost",
        "stop_gained",
    ]
    if impact in impacts:
        return True
    else:
        return False
