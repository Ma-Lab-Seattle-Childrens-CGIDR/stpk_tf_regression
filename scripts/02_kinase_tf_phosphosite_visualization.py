"""
Script for generating plots representing the results of the phosphosite
regression for all of the
"""
# DEPENDENCIES: regression/01_kinase_tf_phosphosite_regression

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
from typing import cast, Optional

# External Imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Local Imports

# Path Setup
try:
    BASE_PATH = pathlib.Path(__file__).parent.parent
except NameError:
    BASE_PATH = pathlib.Path(".").absolute()
DATA_PATH = BASE_PATH / "data"
RESULTS_PATH = BASE_PATH / "results"
VISUALIZATION_OUT_PATH = RESULTS_PATH / "visualization"
LOG_PATH = BASE_PATH / "logs"

# Create directories if needed
RESULTS_PATH.mkdir(parents=True, exist_ok=True)
VISUALIZATION_OUT_PATH.mkdir(parents=True, exist_ok=True)
LOG_PATH.mkdir(parents=True, exist_ok=True)

# Logging setup
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "02_kinase_tf_phosphosite_visualization.log",
    filemode="w",
    level=logging.INFO,
)

# Stop Pandas printing erroneous errors (the chained assignment
# is false positive for this file)
pd.options.mode.chained_assignment = None

# Read in the mycobrowser data, and create a dict to rename the Rv numbers to
# gene names
locus_to_name_dict = cast(
    dict[str, str],
    pd.read_csv(
        DATA_PATH / "Mycobacterium_tuberculosis_H37Rv_txt_v5.txt",
        sep="\t",
    )
    .set_index("Locus")["Name"]
    .to_dict(),
)


# Read in the regression results dataframe
logger.info("Reading in the regression results")
regression_results_df = pd.read_csv(
    RESULTS_PATH / "phosphosite_linear_regression_results.csv"
)

# Get a list of te transcription factors
tf_list = list(np.unique(regression_results_df["TF"]))

KINASE_LIST = [
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


# Create a function for generating the output plot
def plot_regression_results(
    tf_id: str,
    regression_results_df: pd.DataFrame,
    outdir: pathlib.Path,
    pval_col,
    rename_dict: Optional[dict[str, str]] = None,
    height_ratios: Optional[list[float]] = None,
    width_ratios: Optional[list[float]] = None,
    figsize: Optional[tuple[float, float]] = None,
    vmin_coef: Optional[float] = None,
    vmax_coef: Optional[float] = None,
    vmin_sig: Optional[float] = None,
    vmax_sig: Optional[float] = None,
    robust: bool = True,
    global_pval_adj: bool = True,
    format: str = "png",
):
    """
    Plot the results of the phosphosite regression for a single transcription factor

    Parameters
    ----------
    tf_id : str
        The Rv number of the transcription factor to plot the results for
    regression_results_df : pd.DataFrame
        The results of the phosphosite regression used for plotting
    outdir : pathlib.Path
        The directory to save the output figure
    pval_col : str
        Column to use for significance
    rename_dict : dict of str to str, optional
    height_ratios : list of float, optional
        Height ratios to be passed to matplotlib.pyplot.subplots
    width_ratios : list of float, optional
        Width ratios to be passed to matplotlib.pyplot.subplots
    figsize : tuple of float
        Figure size to be used
    vmin_coef,vmax_coef : float, optional
        Values to anchor the coefficient colormap, if not provided they
        are inferred from the data
    vmin_sig,vmax_sig : float, optional
        Values to anchor the significance colormap, if not provided they
        are inferred from the data
    robust : bool, default=True
        If vmin, or vmax are missing for a heatmap compute the
        colormap range based on robust quantiles rather than extreme values
    format : str, default='png'
        Format to use when saving the figure
    """
    if not rename_dict:
        rename_dict = {}
    if not figsize:
        figsize = (6.4, 4.8)
    tf_res_df: pd.DataFrame = regression_results_df[
        regression_results_df["TF"] == tf_id
    ]
    tf_res_df["phosphosite"] = tf_res_df["phosphosite"].str.replace(
        f"^{tf_id}_", "", regex=True
    )
    # Construct a dataframe for plotting the phosphosite information
    phosphosite_plot_df: pd.DataFrame = tf_res_df.copy()[
        ["gene", "phosphosite", "interaction_coef", pval_col]
    ]
    # Rename the columns to reflect what they will be in the plot
    phosphosite_plot_df = phosphosite_plot_df.rename(
        {
            "phosphosite": "term",
            "interaction_coef": "coefficient",
            pval_col: "adjusted p-value",
        },
        axis=1,
    )
    # Calculate the significance/confidence of the predictions
    phosphosite_plot_df["significance"] = -np.log10(
        phosphosite_plot_df["adjusted p-value"]
    )
    phosphosite_plot_df = phosphosite_plot_df.drop("adjusted p-value", axis=1)
    phosphosite_plot_df = phosphosite_plot_df.drop_duplicates()
    # Construct a dataframe for plotting the TF information
    tf_plot_df: pd.DataFrame = tf_res_df.copy()[["gene", "TF", "tf_coef"]]
    # Rename the columns to reflect their meaning in the plot
    tf_plot_df = tf_plot_df.rename(
        {"TF": "term", "tf_coef": "coefficient"}, axis=1
    )
    tf_plot_df["significance"] = (
        0.0  # The tf-coefficient doesn't have it's p-value recorded
    )
    tf_plot_df = tf_plot_df.drop_duplicates()
    # Combine the two dataframes into the plot_df
    plot_df = pd.concat([phosphosite_plot_df, tf_plot_df], axis=0)
    # Create the matplotlib figure which will be used for plotting
    fig, axs = plt.subplot_mosaic(
        [
            [
                "coefficients",
                "coefficients_cmap",
                "significance",
                "significance_cmap",
            ],
            ["coef_phosphorylation", ".", "sig_phosphorylation", "."],
        ],
        height_ratios=height_ratios,
        width_ratios=width_ratios,
        layout="tight",
        figsize=figsize,
    )
    # Setup axis sharing
    # Sharing x vertically
    axs["coefficients"].sharex(axs["coef_phosphorylation"])
    axs["significance"].sharex(axs["sig_phosphorylation"])
    # Sharing y horizontally
    axs["coefficients"].sharey(axs["significance"])
    axs["coef_phosphorylation"].sharey(axs["sig_phosphorylation"])
    tf_name = rename_dict[tf_id] if tf_id in rename_dict else tf_id
    fig.suptitle(f"Impact of Phosphorylation for {tf_name}", fontsize=24)
    # Use seaborn to plot the coefficient values in the top left plot
    coefficient_df = plot_df.pivot(
        columns="term", index="gene", values="coefficient"
    ).rename(rename_dict, axis=0)
    coefficient_ax = axs["coefficients"]
    sns.heatmap(
        data=coefficient_df,
        ax=coefficient_ax,
        cmap=sns.diverging_palette(220, 10, as_cmap=True),
        cbar_ax=axs["coefficients_cmap"],
        cbar_kws={
            "shrink": 0.6,
        },
        center=0.0,
        robust=robust,
        vmin=vmin_coef,
        vmax=vmax_coef,
    )
    coefficient_ax.set_title("Coefficient", fontsize=18)
    coefficient_ax.set_ylabel("Genes", fontsize=18)
    # Plot the significance
    significance_df = plot_df.pivot(
        columns="term", index="gene", values="significance"
    ).rename(rename_dict, axis=0)
    significance_ax = axs["significance"]
    sns.heatmap(
        data=significance_df,
        ax=significance_ax,
        cmap="magma",
        cbar_ax=axs["significance_cmap"],
        cbar_kws={
            "shrink": 0.6,
        },
        robust=robust,
        vmin=vmin_sig,
        vmax=vmax_sig,
    )
    significance_ax.set_title("Confidence", fontsize=18)
    significance_ax.set_ylabel("")
    # Common settings for the top two plots
    for ax in [axs["coefficients"], axs["significance"]]:
        ax.set_xlabel("")
        ax.set_yticks(
            list(range(len(coefficient_df.index))),
            list(coefficient_df.index),
            rotation=0,
            ha="right",
        )
        ax.set_xticks(
            list(range(len(coefficient_df.columns))),
            list(coefficient_df.columns),
            rotation=90,
            ha="right",
        )

    # Now create the dataframe that will be used to add information on which kinase targets the phosphosite
    phosphosite_list = list(np.unique(tf_res_df["phosphosite"]))
    phosphosite_list.append(tf_id)
    phosphosite_list = sorted(phosphosite_list)
    xx, yy = np.meshgrid(phosphosite_list, KINASE_LIST)
    xx = xx.ravel()
    yy = yy.ravel()
    phos_plot_df = pd.DataFrame(
        {"term": xx, "kinase": yy, "phosphorylates": 0.0}
    )
    for x, y in zip(xx, yy):
        if (
            tf_res_df.loc[tf_res_df["phosphosite"] == x, y].all()
            and len(tf_res_df.loc[tf_res_df["phosphosite"] == x]) > 0
        ):
            phos_plot_df.loc[
                (phos_plot_df["term"] == x) & (phos_plot_df["kinase"] == y),
                "phosphorylates",
            ] = 1.0
    for ax in [axs["coef_phosphorylation"], axs["sig_phosphorylation"]]:
        sns.scatterplot(
            data=phos_plot_df,
            x="term",
            y="kinase",
            hue="phosphorylates",
            hue_norm=(0.0, 1.0),
            palette="binary",
            ax=ax,
            legend=False,
            s=300,
        )
        ax.set_xticks(
            list(range(len(phosphosite_list))),
            phosphosite_list,
            rotation=90,
            ha="right",
        )
        ax.set_ylabel("STPK", fontsize=18)
        ax.set_xlabel("Phosphorylation Site", fontsize=18)
    # Now that the plot has been generated, save it to the desired output
    format = "png"
    fig.savefig(
        outdir / f"{tf_name}.{format}",
        format=format,
        pad_inches=10.0,
    )
    plt.close(fig)


# Generate the visualization for every tf
for tf in tf_list:
    gene_count = len(
        np.unique(
            regression_results_df[regression_results_df["TF"] == tf]["gene"]
        )
    )
    phosphosite_count = 1 + len(
        np.unique(
            regression_results_df[regression_results_df["TF"] == tf][
                "phosphosite"
            ]
        )
    )
    figsize = (
        (phosphosite_count + 1) * 2 * 1,
        (gene_count + len(KINASE_LIST) + 2) * 0.5,
    )
    width_additional_scale = 10
    gene_additional_scale = 4
    plot_regression_results(
        tf_id=tf,
        regression_results_df=regression_results_df,
        pval_col="interaction_pvalue_global_adj",
        height_ratios=[gene_count + gene_additional_scale, len(KINASE_LIST)],
        width_ratios=[
            phosphosite_count + width_additional_scale,
            1,
            phosphosite_count + width_additional_scale,
            1,
        ],
        figsize=figsize,
        robust=True,
        vmin_sig=0.0,
        vmax_sig=2.0,
        rename_dict=locus_to_name_dict,
        outdir=VISUALIZATION_OUT_PATH,
        format="png",
    )
