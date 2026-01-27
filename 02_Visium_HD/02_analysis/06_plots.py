import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import gc

matplotlib.use('Agg')

# Settings
# Platform
VISIUM_HD = True

# Plot types
# The HI_RES_DPI is required as the quality of the HE image behind the vectorized expression points is based on the DPI.
HI_RES = True
LOW_RES = True
HI_RES_DPI = 3000
LOW_RES_DPI = 1250
GENE_PLOTTING = True
HE_PLOTTING = False
USAGE_PLOTTING = False

BINARIZE = False
GENE_LIST = ["IL1B", "CCL3", "CXCL1", "CXCL2", "CXCL3", "CXCL5",
            "IDO1", "IDO2", "ISG15", "CXCL8", "CXCL9", "CXCL10",
            "ARG1", "MRC1", "CD274", "CX3CR1",
            "MKI67", "CDK1",
            "LYVE1", "HES1", "FOLR2",
            "VEGFA", "SPP1",
            "APOC1", "APOE", "ACP5", "FABP5", "TREM1", "TREM2", "SPP1",
            "GPNMB", "TFE3", "ASPSCR1", "MRC1", "CD3D", "CD3E", "CD3G", "CD4",
            "CD8A", "CD19", "COL6A1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "HIF1A",
            "VWF", "PECAM1", "CD34", "PECAM1", "MS4A1", "TNC", "FN1", "MYH1",
            "PEX7", "GCAR", "CD274", "PDCD1", "ATF1", "HMGB2", "EIF2A", "CYC1",
            "ANXA6", "GUK1", "CD44", "LGALS3BP", "CD68",
            "CTLA4", "CD80", "CD86",
            "TIGIT", "PVR", "NECTIN2",
            "KLRB1", "CLEC2D",
            "FAS", "FASLG"]

def load_metadata(clusters_path, cluster_key=None):
    """
    Load the metadata file and alter it as needed. Kept cluster_key from original code, but unsure of purpose.
    """
    metadata = pd.read_csv(clusters_path)
    if (cluster_key != None):
        # We only need to use this function if we need to add a column to the metadata file that is based on the cluster_key. For example,
        # if we forgot to add a column to the metadata file that specifies what sample each barcode belongs to, we can use this function to
        # split the barcode name and add a column to the metadata file that specifies the sample.

        pass
        # metadata['Sample'] = metadata[cluster_key].astype('category')
        # metadata['Spots'] = metadata['Spots'].apply(lambda x: x.split('_')[1] if '_' in x else x)
        # metadata['Spots'] = metadata['Spots'].apply(lambda x: x.split('-')[2]+'-1' if len(x) > 18 else x)
        # metadata['Spots'] = metadata['Spots'].apply(lambda x: x.split('_')[1]+'')
        # metadata['Spots'] = metadata['Spots'].apply(lambda x: x.split('-')[0]+'-1')
    
    return metadata

def create_figure(ROWS, COLS):
    axes_hist = []
    figure, axes = plt.subplots(nrows=ROWS, ncols=COLS, figsize=(int(5 * COLS), int(5 * ROWS)))

    if ROWS == 1 and COLS == 1:
        axes_hist.append(axes)
    else:
        for row in axes:
            if (ROWS > 1):
                for axe in row:
                    axes_hist.append(axe)
            else:
                axes_hist.append(row)

    [axes_hist[x].legend().set_visible(False) for x in range(0, len(axes_hist))]
    [axes_hist[x].get_xaxis().set_visible(False) for x in range(0, len(axes_hist))]
    [axes_hist[x].get_yaxis().set_visible(False) for x in range(0, len(axes_hist))]

    if GENE_PLOTTING or HE_PLOTTING:
        color_scale = "OrRd"
    else:
        color_scale = LinearSegmentedColormap.from_list("", ['cornsilk', 'magenta', 'navy'])

    return figure, axes_hist, color_scale

def create_rank_plots(usage_norm, metadata_df, samples_df, SAMPLE_DIRECTORY_PATH, metagene, MAX_CUTOFF, ROWS, COLS, plot_key):
    figure, axes_hist, color_scale = create_figure(ROWS, COLS)

    for index, row in samples_df.iterrows():
        sample = row["sample_name"]
        
        # The sample_abbreviated is used when we forget to use the full sample name submitted to SpaceRanger as a prefix for the barcodes. 
        sample_abbreviated = row["sample"]

        sample_path = os.path.join(SAMPLE_DIRECTORY_PATH, sample, "outs/binned_outputs/square_024um/")
        sample_read_data = sc.read_visium(sample_path)
        sample_read_data.var_names_make_unique()
        sc.pp.filter_cells(sample_read_data, min_counts=1)

        sample_read_data.raw = sc.pp.log1p(sample_read_data, copy=True)

        if "sample" in metadata_df.columns:
            metadata_subset = metadata_df[metadata_df["sample"] == sample_abbreviated]
        else:
            metadata_subset = metadata_df

        # Subsetting the usage_norm dataframe to only include the barcodes that are in the metadata file
        usage_norm_sub = usage_norm.loc[metadata_subset.index]

        # If we run several samples through factorization, each barcode repeats, so we make them unique by adding a prefix to each 
        # barcode with the sample name. This means that the barcodes in the usage and metadata file will match, but each sample doesn't
        # contain a sample prefix in the barcode name. This is because each sample has a unique barcode identified, but again, it isn't
        # unique between samples. So, we need to now strip this barcode prefix from the usage and metadata files, so that they match the sample.
        if sample_abbreviated in metadata_subset["Spots"].iloc[0]:
            metadata_subset["Spots"] = metadata_subset["Spots"].apply(lambda x: x.split(sample_abbreviated + SEPERATOR)[1])
        
        # Set the usage_norm_sub index to the metadata_subset Spots column, which should be 1:1.
        usage_norm_sub.index = list(metadata_subset['Spots'])

        sample_read_data.obs = pd.concat([sample_read_data.obs, usage_norm_sub.iloc[:, metagene - 1]], axis=1)

        VMIN = 0
        if (sorted(list(usage_norm_sub.iloc[:, metagene - 1]), reverse=True)[0] == 0):
            VMAX = 0.0001
        else:
            VMAX = sorted(list(usage_norm_sub.iloc[:, metagene - 1]), reverse=True)
            VMAX = [x for x in VMAX if str(x) != "nan"]
            VMAX = VMAX[MAX_CUTOFF]

        cmap = plt.get_cmap(color_scale)
        cbar = plt.colorbar(cm.ScalarMappable(norm = None, cmap = cmap), ax = axes_hist[index], ticks = [0, 1], shrink = 0.8, pad = 0.03)
        cbar.ax.set_yticklabels([str(0.00), '%s' % float('%.1g' % VMAX)], fontsize=15)

        sc.pl.spatial(sample_read_data, img_key=None, color='Usage_' + str(metagene), show=False, title=sample,
                  vmin=VMIN, vmax=VMAX, color_map=color_scale, use_raw=True, ax=axes_hist[index], 
                  legend_loc=None, colorbar_loc=None, na_color="black", na_in_legend=False)

    plt.suptitle(title_key)
    figure.tight_layout()

    if(save):
        if HI_RES:
            figure.savefig(os.path.join(save_dir,"hires_vectors", f"{plot_key}.pdf"), bbox_inches='tight', dpi = HI_RES_DPI)

        if LOW_RES:
            figure.savefig(os.path.join(save_dir,"lowres_jpeg",f"{plot_key}.jpeg"), bbox_inches='tight', dpi = LOW_RES_DPI)

        # Clear the current axes.
        plt.cla()

        # Clear the current figure.
        plt.clf()

        # Closes all the figure windows.
        plt.close('all')
        plt.close(figure)
        gc.collect()

    else:
        plt.show()

def create_gene_plots(metadata_df, samples_df, SAMPLE_DIRECTORY_PATH, gene, MAX_CUTOFF, ROWS, COLS, plot_key):
    figure, axes_hist, color_scale = create_figure(ROWS, COLS)

    for index, row in samples_df.iterrows():
        sample = row["sample_name"]
        sample_abbreviated = row["sample"]

        sample_path = os.path.join(SAMPLE_DIRECTORY_PATH, sample, "outs/binned_outputs/square_024um/")
        sample_read_data = sc.read_visium(sample_path)
        sample_read_data.var_names_make_unique()

        if "Sample" in metadata_df.columns:
            metadata_subset = metadata_df[metadata_df["Sample"] == sample_abbreviated]
        else:
            metadata_subset = metadata_df

        VMIN = 0
        if (sorted(list(sample_read_data[:, gene].X.todense()), reverse = True)[0] == 0):
            VMAX = 0.0001
        else:
            VMAX = sorted(list(sample_read_data[:, gene].X.todense()), reverse = True)[MAX_CUTOFF]

        cmap = plt.get_cmap(color_scale)
        cbar = plt.colorbar(cm.ScalarMappable(norm = None, cmap = cmap), ax = axes_hist[index], ticks = [0, 1], shrink = 0.8, pad = 0.03)
        cbar.ax.set_yticklabels([str(0.00), '%s' % float('%.1g' % VMAX)], fontsize=15)

        sc.pl.spatial(sample_read_data, img_key="hires", color = gene, show=False, title = sample, size = 1.3, 
                      vmin = VMIN, vmax = VMAX, use_raw = False, color_map = color_scale, ax=axes_hist[index], 
                      legend_loc=None, colorbar_loc=None)

    plt.suptitle(title_key)
    figure.tight_layout()

    if (save):
        if HI_RES:
            figure.savefig(os.path.join(save_dir,"hires_vectors", "gene_expression", f"{plot_key}.pdf"), bbox_inches='tight', dpi = HI_RES_DPI)
        if LOW_RES:
            figure.savefig(os.path.join(save_dir,"lowres_jpeg", "gene_expression", f"{plot_key}.jpeg"), bbox_inches='tight', dpi = LOW_RES_DPI)

        # Clear the current axes.
        plt.cla()

        # Clear the current figure.
        plt.clf()

        # Closes all the figure windows.
        plt.close('all')
        plt.close(figure)
        gc.collect()

    else:
        plt.show()

def create_HE_plots(metadata_df, samples_df, SAMPLE_DIRECTORY_PATH, ROWS, COLS, plot_key):
    figure, axes_hist, _  = create_figure(ROWS, COLS)

    for index, row in samples_df.iterrows():
        sample = row["sample_name"]
        sample_abbreviated = row["sample"]

        sample_path = os.path.join(SAMPLE_DIRECTORY_PATH, sample, "outs/binned_outputs/square_024um/")
        sample_read_data = sc.read_visium(sample_path)
        sample_read_data.var_names_make_unique()

        sample_read_data.raw = sc.pp.log1p(sample_read_data, copy=True)

        sc.pp.filter_cells(sample_read_data, min_counts=1)

        if "Sample" in metadata_df.columns:
            metadata_subset = metadata_df[metadata_df["Sample"] == sample_abbreviated]
        else:
            metadata_subset = metadata_df

        # To just show the HE stain, do the following:
        # size = 0; colorbar_loc = None
        sc.pl.spatial(sample_read_data, img_key = "hires", show = False, title = sample,
                      use_raw = False, ax = axes_hist[index], legend_loc = None,
                      colorbar_loc = None, size = 0)

    plt.suptitle("Reference HE Images")

    plt.suptitle(title_key)
    figure.tight_layout()

    if (save):
        figure.savefig(os.path.join(save_dir,f"{plot_key}.pdf"), bbox_inches='tight', dpi = HI_RES_DPI)

        # Clear the current axes.
        plt.cla()

        # Clear the current figure.
        plt.clf()

        # Closes all the figure windows.
        plt.close('all')
        plt.close(figure)
        gc.collect()

    else:
        plt.show()

def check_compatability(usage_df, metadata_df, samples_df, samples_path, intersect):
    for index, row in samples_df.iterrows():

        print("\nThe separators in the following should match, or this script will throw a massive error.")
        print("Folder Name: %s" % row["sample_name"])
        print("SpaceRanger Output Names: %s" % samples_path)
        print("Metadata File: %s" % metadata_df.index[1])
        print("cNMF File Index: %s" % usage_df.index[1])

        sample = row["sample_name"]
        if "Sample" in intersect.columns:
            if sample in intersect.Sample:
                pass
            else:
                print("Sample Name %s not in metadata" % sample)
            print("\n")

if __name__ == "__main__":
    RANK = int(sys.argv[1])
    USAGE_FILE = sys.argv[2]
    METADATA_FILE = sys.argv[3]
    SAMPLES_FILE = sys.argv[4]
    SAMPLE_DIRECTORY_PATH = sys.argv[5]
    PLOT_DIRECTORY_PATH = sys.argv[6]
    ROWS = int(sys.argv[7])
    COLS = int(sys.argv[8])
    MAX_CUTOFF = int(sys.argv[9])
    SEPERATOR = sys.argv[10]

    print(RANK)
    print(USAGE_FILE)
    print(METADATA_FILE)
    print(SAMPLES_FILE)
    print(SAMPLE_DIRECTORY_PATH)
    print(PLOT_DIRECTORY_PATH)
    print(ROWS)
    print(COLS)
    print(MAX_CUTOFF)
    print(SEPERATOR)

    usage_df = pd.read_csv(USAGE_FILE, sep='\t', index_col=0)

    metadata_df = load_metadata(METADATA_FILE)
    metadata_df.index = list(metadata_df['Spots'])

    samples_df = pd.read_csv(SAMPLES_FILE)

    save_dir = str(PLOT_DIRECTORY_PATH)
    save = True

    # The following code finds common rows in the metadata and usages files and removes the unique rows in the metadata and usage dataframes.
    usage_meta_intersect = list(set.intersection(set(usage_df.index), set(metadata_df.index)))
    metadata_intersect = metadata_df.loc[usage_meta_intersect]
    usage_intersect = usage_df.loc[usage_meta_intersect]

    # Add "Usage_" as a prefix to the column names if it isn't already present
    usage_intersect.columns = ["Usage_" + str(col) for col in usage_df.columns if "Usage_" not in usage_df.columns[:]]

    # Normalize the usage
    usage_norm = usage_intersect.div(usage_intersect.sum(axis=1), axis=0)

    check_compatability(usage_df, metadata_df, samples_df, SAMPLE_DIRECTORY_PATH, metadata_intersect)

    if GENE_PLOTTING:
        for gene in GENE_LIST:
            title_key = str(gene)
            create_gene_plots(metadata_df, samples_df, SAMPLE_DIRECTORY_PATH, gene, MAX_CUTOFF, ROWS, COLS, title_key)
    if HE_PLOTTING:
        title_key = "HE_stain"
        create_HE_plots(metadata_df, samples_df, SAMPLE_DIRECTORY_PATH, ROWS, COLS, title_key)
    if USAGE_PLOTTING:
        for metagene in range(1, RANK+1):
            title_key = "K" + str(RANK) + "_metagene" + str(metagene)
            create_rank_plots(usage_norm, metadata_df, samples_df, SAMPLE_DIRECTORY_PATH, metagene, MAX_CUTOFF, ROWS, COLS, title_key)
