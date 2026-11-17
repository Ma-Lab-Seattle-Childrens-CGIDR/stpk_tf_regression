"""
Script to evaluate the impact of phosphosite phosphorylation state
on TF activity
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
import sys
from collections import defaultdict

# External Imports
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

# Local Imports
from common_functions import (
    fdr_with_nan,
    find_connected_components,
    rpkm_to_tpm,
)
from sklearn.preprocessing import StandardScaler

# Path Setup
try:
    BASE_PATH = pathlib.Path(__file__).parent.parent
except NameError:
    BASE_PATH = pathlib.Path(".").absolute()
DATA_PATH = BASE_PATH / "data"
RESULTS_PATH = BASE_PATH / "results"
LOG_PATH = BASE_PATH / "logs"


if __name__ == "__main__":
    # Create the directories if needed
    RESULTS_PATH.mkdir(parents=True, exist_ok=True)
    LOG_PATH.mkdir(parents=True, exist_ok=True)

    # Logging setup
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename=LOG_PATH / "01_kinase_tf_phosphosite_regression.log",
        filemode="w",
        level=logging.INFO,
    )
    # SECTION: Check if data downloaded
    phosphorylation_data_path = DATA_PATH / "41564_2022_1313_MOESM3_ESM.xlsx"
    if not phosphorylation_data_path.exists():
        phosphorylation_data_path = (
            DATA_PATH / "NIHMS1988244-supplement-Supplementary_Table_2.xlsx"
        )
    if not phosphorylation_data_path.exists():
        print(
            "Unable to access phosphorylation data, please download 'Supplementary Table 2' "
            "from 'The Mycobacterium tuberculosis protein O-phosphorylation landscape' (PMID: 36690861, "
            "doi: https://doi.org/10.1038/s41564-022-01313-7), and place it in the 'data' directory",
            file=sys.stderr,
        )
        sys.exit(-1)

    stpk_diff_expr_data_path = DATA_PATH / "41564_2022_1313_MOESM4_ESM.xlsx"
    if not stpk_diff_expr_data_path.exists():
        stpk_diff_expr_data_path = (
            DATA_PATH / "NIHMS1988244-supplement-Supplementary_Table_3.xlsx"
        )
    if not stpk_diff_expr_data_path.exists():
        print(
            "Unable to access differential gene expression data for the STPK altered activity strains, "
            "please download 'Supplementary Table 3' from 'The Mycobacterium tuberculosis protein "
            "O-phosphorylation landscape' (PMID: 36690861, doi: https://doi.org/10.1038/s41564-022-01313-7), "
            "and place it in the 'data' directory",
            file=sys.stderr,
        )
        sys.exit(-1)

    stpk_mutant_gene_expr_data_path = (
        DATA_PATH / "GSE195959_Expression.MTB.GeneData.csv.gz"
    )
    if not stpk_mutant_gene_expr_data_path.exists():
        stpk_mutant_gene_expr_data_path = (
            DATA_PATH / "GSE195959_Expression.MTB.GeneData.csv"
        )
    if not stpk_mutant_gene_expr_data_path.exists():
        print(
            "Unable to access gene expression data for the STPK mutant strains, "
            "please download 'GSE195959_Expression.MTB.GeneData.csv.gz' from "
            "GEO GSE195959 at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195959",
            file=sys.stderr,
        )
        sys.exit(-1)

    # SECTION: TF Targets
    # Read in the TFOE data to identify genes targeted by
    # the various transcirption factors
    logger.info("Reading in data for the TF targets")
    # fhaA (Rv0020c) and fhaB (Rv0019c) are not well
    # supported as transcriptional regulators,
    # so they will be removed from the TFs
    columns_to_remove = ["Rv0019c", "Rv0020c"]

    fc_df = pd.read_excel(
        DATA_PATH / "tfoe_targets.xlsx",
        sheet_name="SupplementaryTableS2",
        skiprows=list(range(8)) + [9],
        usecols="A,E:HB",
        index_col=0,
    ).drop(columns_to_remove, axis=1)
    pval_df = pd.read_excel(
        DATA_PATH / "tfoe_targets.xlsx",
        sheet_name="SupplementaryTableS2",
        skiprows=list(range(8)) + [9],
        usecols="A,HC:OZ",
        index_col=0,
    )
    pval_df.columns = pval_df.columns.str.replace(".1", "")
    pval_df = pval_df.drop(columns_to_remove, axis=1)
    tf_target_df: pd.DataFrame = (fc_df.abs() >= 1) & (pval_df <= 0.05)

    # SECTION: Kinase Expression Data
    logger.info("Reading in the kinase perturbation expression data")
    kinase_expr_rpkm = (
        pd.read_csv(stpk_mutant_gene_expr_data_path)
        .drop("GrowthPhase", axis=1)
        .set_index("SampleID")
    )
    logger.info("Converting expression data to tpm")
    kinase_expr_tpm = rpkm_to_tpm(kinase_expr_rpkm)
    kinase_expr_log2_tpm = pd.DataFrame(
        np.log2(kinase_expr_tpm + 1),
        index=kinase_expr_tpm.index,
        columns=kinase_expr_tpm.columns,
    )
    scaler = StandardScaler()
    scaler.set_output(transform="pandas")
    kinase_expr_log2_tpm = scaler.fit_transform(kinase_expr_log2_tpm)

    # WT10 and WT11 had some experimental issues, so removing them
    kinase_expr_log2_tpm = kinase_expr_log2_tpm.drop(["WT10", "WT11"], axis=0)
    kinase_expr_log2_tpm.index = kinase_expr_log2_tpm.index.str.replace(
        "KO", "LOF"
    )
    kinase_expr_log2_tpm.index = kinase_expr_log2_tpm.index.str.replace(
        "KD", "LOF"
    )
    kinase_expr_log2_tpm.index = kinase_expr_log2_tpm.index.str.replace(
        "IND", "OE"
    )

    # SECTION: Identify common genes for filtering out genes/tfs not shared by
    # both datasets
    logger.info("Identifying genes common to both datasets")
    frando_stpk_genes = set(kinase_expr_log2_tpm.columns)
    # NOTE: TFOE data has TFs as columns, genes as rows
    rustad_tfoe_genes = set(fc_df.index)
    common_genes = frando_stpk_genes.intersection(rustad_tfoe_genes)

    # Create a list of transcription factors
    logger.info("Getting a list of transcription factors")
    tf_list = [tf for tf in fc_df.columns if tf in common_genes]
    tf_set = set(tf_list)  # For faster membership testing

    # SECTION: Kinase Information
    logger.info("Creating kinase information lists")
    kinase_id_list = [
        "PknB",
        "PknD",
        "PknE",
        "PknF",
        "PknG",
        "PknH",
        "PknI",
        "PknJ",
        "PknK",
        "PknL",
    ]
    kinase_id_to_locus = {
        "PknB": "Rv0014c",
        "PknD": "Rv0931c",
        "PknE": "Rv1743",
        "PknF": "Rv1746",
        "PknG": "Rv0410c",
        "PknH": "Rv1266c",
        "PknI": "Rv2914c",
        "PknJ": "Rv2088",
        "PknK": "Rv3080c",
        "PknL": "Rv2176",
    }
    kinase_locus_list = [locus for _, locus in kinase_id_to_locus.items()]

    # SECTION: Target dictionaries

    # Create a dict of TF to phosphosites
    # also a list of phosphosites, and a list of TFs targeted
    # by at least 1 kinase
    phosphorylated_tf_set = set()
    tf_phosphosite_set = set()
    tf_to_phosphosite_set_dict: dict[str, set] = defaultdict(set)
    phosphosite_to_kinase_dict: dict[str, set[str]] = defaultdict(set)
    for kinase in kinase_id_list:
        kinase_target_df = pd.read_excel(
            phosphorylation_data_path, sheet_name=kinase
        )
        for mutant, mutant_df in kinase_target_df.groupby("Mutant"):
            direct_target_df = mutant_df[mutant_df["Indirect?"] != "Indirect"]
            indirect_target_df = mutant_df[mutant_df["Indirect?"] == "Indirect"]
            # Go through each TFs direct targets
            for phosphosite in direct_target_df[
                direct_target_df["Rv Number"].isin(tf_set)
            ]["Phosphosite"]:
                tf, *rest = phosphosite.split("_")
                site_location = "_".join(rest)
                direction = "hyper" if mutant == "OE" else "hypo"
                phosphorylated_tf_set.add(tf)
                psite = f"{tf}_{direction}_{site_location}"
                tf_phosphosite_set.add(psite)
                tf_to_phosphosite_set_dict[tf].add(psite)
                phosphosite_to_kinase_dict[
                    f"{tf}_{direction}_{site_location}"
                ].add(kinase)

            for phosphosite in indirect_target_df[
                indirect_target_df["Rv Number"].isin(tf_set)
            ]["Phosphosite"]:
                tf, *rest = phosphosite.split("_")
                site_location = "_".join(rest)
                direction = "hypo" if mutant == "OE" else "hyper"
                phosphorylated_tf_set.add(tf)
                psite = f"{tf}_{direction}_{site_location}"
                tf_phosphosite_set.add(psite)
                tf_to_phosphosite_set_dict[tf].add(psite)
                phosphosite_to_kinase_dict[
                    f"{tf}_{direction}_{site_location}"
                ].add(kinase)
    # Create a dataframe with the phosphosite and the kinases that target it as columns
    psite_to_kinase_df_psite_list = []
    psite_to_kinase_df_kinase_list = []
    for psite, kinase_set in phosphosite_to_kinase_dict.items():
        psite_to_kinase_df_psite_list.append(psite)
        psite_to_kinase_df_kinase_list.append(", ".join(sorted(kinase_set)))
    psite_to_kinase_df = pd.DataFrame(
        {
            "phosphosite": psite_to_kinase_df_psite_list,
            "kinases_targeting": psite_to_kinase_df_kinase_list,
        }
    )

    # Create a dict of gene to TF which target them
    gene_to_tf_dict = defaultdict(list)
    for tf in phosphorylated_tf_set:
        for gene in tf_target_df[tf_target_df[tf]].index:
            if gene not in common_genes:
                continue
            if gene in tf_set or gene in kinase_id_list:
                continue  # Don't want to fit TFs or Kinases to themselves
            gene_to_tf_dict[gene].append(tf)
    targeted_genes_list = list(gene_to_tf_dict.keys())

    # SECTION: Create a dataframe describing the phophorylation state of the
    # different TF phosphosites in each condition
    tf_phosphorylation_df = pd.DataFrame(
        0.0,
        index=kinase_expr_log2_tpm.index,
        columns=pd.Index(sorted(tf_phosphosite_set)),
    )
    for kinase in kinase_id_list:
        kinase_target_df = pd.read_excel(
            phosphorylation_data_path, sheet_name=kinase
        )
        for mutant, mutant_df in kinase_target_df.groupby("Mutant"):
            direct_target_df = mutant_df[mutant_df["Indirect?"] != "Indirect"]
            indirect_target_df = mutant_df[mutant_df["Indirect?"] == "Indirect"]
            # Go through each TFs direct targets
            phosphosite_set = set()
            for phosphosite in direct_target_df[
                direct_target_df["Rv Number"].isin(tf_set)
            ]["Phosphosite"]:
                tf, *rest = phosphosite.split("_")
                site_location = "_".join(rest)
                direction = "hyper" if mutant == "OE" else "hypo"
                phosphosite_set.add(f"{tf}_{direction}_{site_location}")
            # And the indirect targets
            for phosphosite in indirect_target_df[
                indirect_target_df["Rv Number"].isin(tf_set)
            ]["Phosphosite"]:
                tf, *rest = phosphosite.split("_")
                site_location = "_".join(rest)
                direction = "hypo" if mutant == "OE" else "hyper"
                phosphosite_set.add(f"{tf}_{direction}_{site_location}")
            tf_phosphorylation_df.loc[
                tf_phosphorylation_df.index.str.startswith(
                    f"{kinase}_{mutant}"
                ),
                list(phosphosite_set),
            ] = 1.0

    # Some of these phosphorylations modulations always occur together,
    # and those which perfectly line up will be filtered down to a single
    # representative to increase the power of the linear regression

    # Create a function which will find connected components of
    # the phosphorylation sites for a given transcription factor

    def tf_phosphosite_grouping(tf_name: str, tf_phosphorylation: pd.DataFrame):
        """
        Find the phosphosite groups for the provided tf

        Parameters
        ----------
        tf_name : str
            Name of the TF to find the phosphosite groupings for
        tf_phosphorylation : pd.DataFrame
            Dataframe descibing the phosphorylation state of the TF,
            columns should represent different phosphosites


        Returns
        -------
        phosphosite_groups : dict of str to set of str
            Phosphosite groupings, keyed by one member of the group
        """
        tf_psite_data = tf_phosphorylation.loc[
            :, tf_phosphorylation.columns.str.startswith(tf_name)
        ]
        psite_correlations = tf_psite_data.corr()
        node_list = tf_psite_data.columns.to_list()
        edge_list: list[tuple[str, str]] = [
            (psite_correlations.index[i], psite_correlations.index[j])
            for i, j in np.argwhere(psite_correlations > 0.95)
        ]
        components = {
            min(c): c - {min(c)}
            for c in find_connected_components(
                node_list=node_list, edge_list=edge_list
            )
        }
        return components

    # Create the dict holding all the TF phosphosite groupings
    tf_psites_to_remove: set[str] = set()
    tf_psite_group_dict: dict[str, set[str]] = {}
    for tf in phosphorylated_tf_set:
        tf_psite_group_dict.update(
            tf_phosphosite_grouping(tf, tf_phosphorylation_df)
        )
    for sites in tf_psite_group_dict.values():
        tf_psites_to_remove.update(sites)

    # Remove the extraneous p-sites from the dataframe and the psite dict
    tf_phosphorylation_df = tf_phosphorylation_df.drop(
        list(tf_psites_to_remove), axis=1
    )
    tf_to_phosphosite_set_dict = {
        k: v - tf_psites_to_remove
        for k, v in tf_to_phosphosite_set_dict.items()
    }

    # SECTION: Linear Model Fitting
    # Use the gene_to_tf_dict and the tf_to_phosphosite_set_dict to create
    # and fit linear models for each gene of the form:
    # gene ~ tf + tf*psite
    def fit_gene_tf_phosphosite_model(
        gene: str,
        expr_data: pd.DataFrame,
        phosphorylation_state: pd.DataFrame,
        gene_to_tf_dict: dict[str, list[str]],
        tf_to_phosphosite_set_dict: dict[str, set[str]],
    ) -> pd.DataFrame:
        """
        Fit a linear model with a single gene as the target, and features for each TF targeting that gene,
        including interaction terms between the TF expression and the phosphorylation state

        Parameters
        ----------
        gene : str
            Rv Number of the gene to fit the model for
        expr_data : pd.DataFrame of floats
            Expression level of the genes across different conditions,
            indexed by condition with columns for all targeted genes
        phosphorylation_state : pd.DataFrame of floats
            The phosphorylation state of the phosphorylation sites for
            the transcription factors across the same conditions as the expr_data,
            with a row for each condition, and columns for each phosphosite
        gene_to_tf_dict : dict of str to list of str
            Dict to identify which TFs target each gene, keys are
            the Rv numbers for the targeted genes, and values
            are lists of Rv numbers of the TFs that target the
            particular gene
        tf_to_phosphosite_set_dict : dict of str to set of str
            Dict describing the phosphorylations of each TF,
            should be indexed by the Rv number of the gene, with
            values that are the location of the phosphosite on that TF

        Returns
        -------
        pd.DataFrame
            Dataframe describing the results of the regression,
            it includes columns for:
                * gene: The gene targeted by the TF
                * TF: The TF involved in the interaction
                * phosphosite: The phosphorylation site
                * tf_coef: The coefficient for the TF alone
                * interaction_coef: The interaction coefficient
                * interaction_pvalue: The pvalue of the interaction term
                * interaction_pvalue_adj: The p-value of the interaction
                  term, adjusted to account for false discovery rate using
                  the Benjamani-Hochberg method
                * rsquared: The R^2 of the linear regression model
        """
        # Start by combining the expr_data and the phosphorylation_state
        data_df = pd.concat([expr_data, phosphorylation_state], axis=1)
        # Get a list of TFs that must be included
        tfs_targeting_gene = gene_to_tf_dict[gene]
        # Get a list of phosphorylation interaction tuples
        # of the form (tf, phosphosite)
        interaction_term_list: list[tuple[str, str]] = []
        for tf in tfs_targeting_gene:
            for psite in tf_to_phosphosite_set_dict[tf]:
                interaction_term_list.append((tf, psite))
        # Build the formula representing the linear regression model
        tf_terms = "+".join(tfs_targeting_gene)
        interaction_terms = "+".join(
            (f"{interaction_tuple[0]}:{interaction_tuple[1]}")
            for interaction_tuple in interaction_term_list
        )
        ols_formula = f"{gene} ~ {tf_terms} + {interaction_terms}"
        # Perform the OLS regression
        ols_res = smf.ols(formula=ols_formula, data=data_df).fit()
        # Extract the coefficients and the p-values
        ols_coefs = ols_res.params
        ols_pvals = ols_res.pvalues
        # Construct the results dataframe
        res_df = pd.DataFrame(
            {
                "gene": pd.Series(dtype="string"),
                "TF": pd.Series(dtype="string"),
                "phosphosite": pd.Series(dtype="string"),
                "tf_coef": pd.Series(dtype="float"),
                "interaction_coef": pd.Series(dtype="float"),
                "interaction_pvalue": pd.Series(dtype="float"),
                "interaction_pvalue_adj": pd.Series(dtype="float"),
                "rsquared": pd.Series(dtype="float"),
            },
            index=pd.Index(f"{it[0]}:{it[1]}" for it in interaction_term_list),
        )
        res_df["gene"] = gene
        res_df["rsquared"] = ols_res.rsquared
        for tf, phosphosite in interaction_term_list:
            int_term = f"{tf}:{phosphosite}"
            res_df.loc[int_term, "TF"] = tf
            res_df.loc[int_term, "phosphosite"] = phosphosite
            res_df.loc[int_term, "tf_coef"] = ols_coefs[tf]
            res_df.loc[int_term, "interaction_coef"] = ols_coefs[int_term]
            res_df.loc[int_term, "interaction_pvalue"] = ols_pvals[int_term]
        # Correct for false discovery rate
        corrected_pvals = fdr_with_nan(res_df["interaction_pvalue"])
        res_df.loc[corrected_pvals.index, "interaction_pvalue_adj"] = (
            corrected_pvals
        )
        return res_df.reset_index(drop=True)

    # Perform the model fitting for all the targeted genes
    logger.info("Starting to fit linear models")
    ols_res_list = []
    for gene in targeted_genes_list:
        ols_res_list.append(
            fit_gene_tf_phosphosite_model(
                gene=gene,
                expr_data=kinase_expr_log2_tpm,
                phosphorylation_state=tf_phosphorylation_df,
                gene_to_tf_dict=gene_to_tf_dict,
                tf_to_phosphosite_set_dict=tf_to_phosphosite_set_dict,
            )
        )
    ols_results = pd.concat(ols_res_list, axis=0, ignore_index=True)
    logger.info("Finished fitting models, adding additional information")
    # Add a column for TF adjusted FDR p-values
    tf_adj_pval_df_list = []
    for tf, df in ols_results.groupby("TF"):
        df["interaction_pvalue_tf_adj"] = fdr_with_nan(df["interaction_pvalue"])
        tf_adj_pval_df_list.append(df)
    ols_results = pd.concat(tf_adj_pval_df_list, axis=0, ignore_index=True)

    # Add a column for phosphosite adjusted FDR p-values
    site_adj_pval_df_list = []
    for site, df in ols_results.groupby("phosphosite"):
        df["interaction_pvalue_site_adj"] = fdr_with_nan(
            df["interaction_pvalue"]
        )
        site_adj_pval_df_list.append(df)
    ols_results = pd.concat(site_adj_pval_df_list, axis=0, ignore_index=True)

    # Add a column for globally adjusted FDR p-values
    ols_results["interaction_pvalue_global_adj"] = fdr_with_nan(
        ols_results["interaction_pvalue"]
    )

    # Add information about which kinases target which phosphosites
    ols_results = pd.merge(
        ols_results,
        psite_to_kinase_df,
        how="left",
        left_on="phosphosite",
        right_on="phosphosite",
    )
    # Add information about the vulnerability index and essentiality of each gene
    vi_df = pd.read_excel(
        DATA_PATH / "bosch_vulnerability_index.xlsx", sheet_name="(1) Mtb H37Rv"
    )
    vi_df["locus_tag"] = vi_df["locus_tag"].str.replace("RVBD", "Rv")
    vi_df = vi_df.rename(
        {
            "locus_tag": "gene",
        },
        axis=1,
    )
    vi_df = vi_df[["gene", "tnseq_ess", "crispr_ess", "Vulnerability Index"]]
    ols_results = pd.merge(
        ols_results, vi_df, how="left", left_on="gene", right_on="gene"
    )
    # Add columns for the individual kinases to indicate which
    # phosphosites they target
    for kinase in kinase_id_list:
        ols_results[kinase] = ols_results["kinases_targeting"].str.contains(
            kinase
        )

    # Add column for if the phosphorylation seems to potentiate the TF effect
    ols_results["phosphorylation_potentiates"] = (
        np.sign(ols_results["tf_coef"])
        * np.sign(ols_results["interaction_coef"])
        > 0.0
    )
    # Add column for the impact on the gene of the phosphorylation event
    ols_results["gene_up"] = ols_results["interaction_coef"] > 0.0

    # Add information for the differential expression by the kinase conditions
    for kinase in kinase_id_list:
        diff_gene_df = pd.read_excel(
            stpk_diff_expr_data_path, sheet_name=kinase
        )
        # First get the OE genes diff
        kinase_gene_diff_oe = diff_gene_df.loc[diff_gene_df["Mutant"] == "OE"]
        kinase_gene_diff_oe = (
            kinase_gene_diff_oe[["DEG", "Fold-change (log2)"]]
            .rename({"DEG": "gene", "Fold-change (log2)": "log2fc"}, axis=1)
            .set_index("gene")
            .add_prefix(f"{kinase}_OE_")
        )
        ols_results = pd.merge(
            ols_results,
            kinase_gene_diff_oe,
            how="left",
            left_on="gene",
            right_index=True,
        )
        # Next with the LOF genes diff
        kinase_gene_diff_lof = diff_gene_df.loc[diff_gene_df["Mutant"] == "LOF"]
        kinase_gene_diff_lof = (
            kinase_gene_diff_lof[["DEG", "Fold-change (log2)"]]
            .rename({"DEG": "gene", "Fold-change (log2)": "log2fc"}, axis=1)
            .set_index("gene")
            .add_prefix(f"{kinase}_LOF_")
        )
        ols_results = pd.merge(
            ols_results,
            kinase_gene_diff_lof,
            how="left",
            left_on="gene",
            right_index=True,
        )

    # Add column for which phosphosites are represented by the
    # representative phosphosite
    representative_phosphosite_df = pd.Series(
        {k: ", ".join(v) for k, v in tf_psite_group_dict.items()}
    ).to_frame(name="represented_psites")
    ols_results = pd.merge(
        ols_results,
        representative_phosphosite_df,
        how="left",
        left_on="phosphosite",
        right_index=True,
    )

    # Add columns with additional information about the genes
    # from mycobrowser
    mycobrowser_df = pd.read_csv(
        DATA_PATH / "Mycobacterium_tuberculosis_H37Rv_txt_v5.txt", sep="\t"
    )
    mycobrowser_df = mycobrowser_df.set_index("Locus")
    ols_results = pd.merge(
        ols_results,
        mycobrowser_df,
        how="left",
        left_on="gene",
        right_index=True,
    )

    logger.info("Saving the final results")
    ols_results.to_csv(
        RESULTS_PATH / "phosphosite_linear_regression_results.csv",
        index=False,
    )
