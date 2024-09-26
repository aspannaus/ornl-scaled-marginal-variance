# Import libraries
import os
import pickle
import warnings

import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl

import geopandas as gpd
from tqdm import tqdm
import invr

from pysal.lib import weights
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.linalg import solve
from scipy.sparse.linalg import spsolve

# Ignore FutureWarnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Matplotlib default settings
mpl.rcParams.update(mpl.rcParamsDefault)


# Utility functions
def get_folders(location):
    """Get list of folders in a directory."""
    return [
        name
        for name in os.listdir(location)
        if os.path.isdir(os.path.join(location, name))
    ]


def clean_svi_data(dataframe, variables):
    """Clean the SVI data."""
    for var in variables:
        # replace -999 with 0
        dataframe[var] = dataframe[var].replace(-999, 0)

        # drop -999
        # dataframe = dataframe[dataframe[var] != -999]
    return dataframe


def filter_svi_data_by_percentile(dataframe, var, pct):
    """Filter the SVI data based on a given percentile."""

    percentile_value = dataframe[var].quantile(pct)
    filtered_df = dataframe[dataframe[var] < percentile_value]

    return filtered_df


def generate_adjacent_counties(dataframe, var, pct):
    """Generate adjacent counties based on given dataframe and variable."""

    # Only select the counties below the given percentile
    percentile_value = dataframe[var].quantile(pct)
    filtered_df = dataframe[dataframe[var] < percentile_value]
    # filtered_df = dataframe #Uncomment this line to use the entire dataframe

    adjacent_counties = gpd.sjoin(
        filtered_df, filtered_df, predicate="intersects", how="left"
    )
    adjacent_counties = adjacent_counties.query("sortedID_left != sortedID_right")
    adjacent_counties = (
        adjacent_counties.groupby("sortedID_left")["sortedID_right"]
        .apply(list)
        .reset_index()
    )
    adjacent_counties.rename(
        columns={"sortedID_left": "county", "sortedID_right": "adjacent"}, inplace=True
    )
    _adjacencies_list = adjacent_counties["adjacent"].tolist()
    _county_list = adjacent_counties["county"].tolist()
    merged_df = pd.merge(
        adjacent_counties, dataframe, left_on="county", right_on="sortedID", how="left"
    )
    merged_df = gpd.GeoDataFrame(merged_df, geometry="geometry")
    return _adjacencies_list, merged_df, _county_list


def form_simplicial_complex(adjacent_county_list, _county_list):
    """Form a simplicial complex based on adjacent counties."""
    max_dimension = 3
    return invr.incremental_vr([], adjacent_county_list, max_dimension, _county_list)


def plot_simplicial_complex(dataframe, vertices, var, _results_path):
    """Plot the simplicial complex."""

    # city centroids
    # sortedID used in persistence image code for the FIPS
    city_coordinates = {
        city.FIPS: np.array((city.geometry.centroid.x, city.geometry.centroid.y))
        for _, city in dataframe.iterrows()
    }

    # Create a figure and axis
    _, ax = plt.subplots(figsize=(8, 8))
    ax.set_axis_off()

    # Plot the "wyoming_svi" DataFrame
    dataframe.plot(ax=ax, edgecolor="black", linewidth=0.3, color="white")

    # Plot the centroid of the large square with values
    for _, row in dataframe.iterrows():
        centroid = row["geometry"].centroid
        # text_to_display = f"FIPS: {row['FIPS']}\nFilteration: {row['EP_SNGPNT']}"
        plt.text(
            centroid.x,
            centroid.y,
            str(row["FIPS"]),
            fontsize=8,
            ha="center",
            color="black",
        )
        # plt.text(centroid.x, centroid.y, text_to_display, fontsize=10, ha='center', color="black")

    for edge_or_traingle in vertices:

        if len(edge_or_traingle) == 2:
            # Plot an edge
            ax.plot(
                *zip(*[city_coordinates[vertex] for vertex in edge_or_traingle]),
                color="red",
                linewidth=1,
            )
            # img = fig2img(fig)
            # list_gif.append(img)
        elif len(edge_or_traingle) == 3:
            # Plot a triangle
            ax.add_patch(
                plt.Polygon(
                    [city_coordinates[vertex] for vertex in edge_or_traingle],
                    color="green",
                    alpha=0.2,
                )
            )
            # img = fig2img(fig)
            # list_gif.append(img)
    # plt.show()
    # save the plot in the results_path
    plt.savefig(f"{_results_path}/simplicial_complex_{var}.png")
    plt.close()


def create_variable_folders(base_path, variables):
    """Create folders for each variable."""
    for _variable in variables:
        os.makedirs(os.path.join(base_path, _variable), exist_ok=True)
    print("Done creating folders for each variable")


def generate_scaled_marginal_variance(simplices, data_frame, variable_name):

    selected_counties = []

    for _set in simplices:
        if len(_set) == 2 or len(_set) == 3:
            for vertice in _set:
                if vertice not in selected_counties:
                    selected_counties.append(vertice)

    filtered_counties_df = data_frame.loc[
        data_frame["sortedID"].isin(selected_counties)
    ]

    # lattice stored in a geo-table
    wq = weights.contiguity.Queen.from_dataframe(filtered_counties_df)
    neighbors_q = wq.neighbors

    QTemp = pd.DataFrame(*wq.full()).astype(int)
    QTemp = QTemp.multiply(-1)

    QTemp.index = filtered_counties_df["sortedID"].values
    QTemp.columns = filtered_counties_df["sortedID"].values

    # for each row in the fullMatrix dataframe sum the values in the row and take the absolute value and store in the diagonal
    for i in QTemp.index:
        QTemp.loc[i, i] = abs(QTemp.loc[i].sum())

    # Marginal variance code -Multiple clusters

    # transform df to numpy array
    Q = QTemp.to_numpy()

    graph = csr_matrix(QTemp)
    n_components, labels = connected_components(
        csgraph=graph, directed=False, return_labels=True
    )

    # print(f"Number of connected components: {n_components}")

    data_frame[variable_name + "_marginal_variance"] = None

    print(f"Number of connected components: {n_components}")

    for k in range(n_components):
        # print(k)

        # get the length of the labels array where the value is equal to i
        # print(len(labels[labels == k]))

        if len(labels[labels == k]) == 1:

            # get the index of the label
            index = np.where(labels == k)[0][0]
            # print(index)

            # this part is not written becase: does not exists

            # # get the index from Q_df
            # print(Q_df.index[index])

            # print(f"Region {k} is an isolated region")
            # print(f"Marginal Variances with FIPS: {list(zip(Qmatrix[0].index, marginal_variances))}")
        else:
            print(f"Region {k} is a connected region\n")

            # get the location index to an array
            index = np.where(labels == k)
            # print(index)

            # Extract the submatrix
            QQ = Q[np.ix_(index[0], index[0])]

            # print(QQ)

            n = QQ.shape[0]

            Q_jitter = QQ + sp.sparse.diags(np.ones(n)) * max(QQ.diagonal()) * np.sqrt(
                np.finfo(np.float64).eps
            )

            # inverse of precision (Q) is cov

            Q_perturbed = sp.sparse.csc_array(Q_jitter)

            b = sp.sparse.identity(n, format="csc")

            sigma = spsolve(Q_perturbed, b)

            # V \in Null(Q)

            V = np.ones(n)  # from pg. 6

            W = sigma @ V.T  # \Sigma * B in 3.17

            Q_inv = sigma - np.outer(W * solve(V @ W, np.ones(1)), W.T)

            # grabbing diag of cov gives var and

            # arithmetic mean in log-space becomes geometric mean after exp

            scaling = np.exp(np.sum(np.log(np.diag(Q_inv))) / n)

            # scaling_factor.append(scaling)

            # print(f"Scaling/GV: {scaling}")

            marginal_variances = np.diag(Q_inv / scaling)
            # print(f"Marginal Variances: {marginal_variances}")

            # print(f"Marginal Variances with FIPS: {list(zip(Qmatrix[0].index[index[0]], marginal_variances))}")

            # # get the Q_df index
            # print(Q_df.index[index[0]])

            # fill the new column with the marginal variances only matching the (Q_df.index[index[0]]
            for sortedID, marginal_variance in zip(
                QTemp.index[index[0]], marginal_variances
            ):
                data_frame.loc[
                    data_frame["sortedID"] == sortedID,
                    variable_name + "_marginal_variance",
                ] = marginal_variance

            # print(data_frame[['sortedID',variable_name+'_marginal_variance']])

    # turn the column into float
    data_frame[variable_name + "_marginal_variance"] = data_frame[
        variable_name + "_marginal_variance"
    ].astype("float64")

    # print(data_frame)
    return data_frame, Q


# Define the main function
# def main():
if __name__ == "__main__":
    # Main execution

    # Download the SVI data for a required state - county level -
    # ESRI Geodatabase (https://www.atsdr.cdc.gov/placeandhealth/svi/data_documentation_download.html)
    svi_gdb_path = "/Users/Shared/ornldev/projects/GitHub/ornl-scaled-marginal-variance/data/SVI2018_WYOMING_COUNTY.gdb"
    results_path = "/Users/Shared/ornldev/projects/GitHub/ornl-scaled-marginal-variance/results"

    # State abbreviation
    state = "WY"

    # Variables of interest/ EP_PCI not included
    selected_variables = [
        "EP_POV",
        "EP_UNEMP",
        "EP_NOHSDP",
        "EP_UNINSUR",
        "EP_AGE65",
        "EP_AGE17",
        "EP_DISABL",
        "EP_SNGPNT",
        "EP_LIMENG",
        "EP_MINRTY",
        "EP_MUNIT",
        "EP_MOBILE",
        "EP_CROWD",
        "EP_NOVEH",
        "EP_GROUPQ",
    ]
    selected_variables_with_censusinfo = (
        ["FIPS", "STCNTY"] + selected_variables + ["geometry"]
    )

    # Selected percentile
    percentile = 0.90

    # Load the SVI data
    svi = gpd.read_file(svi_gdb_path)

    # Clean the SVI data
    svi_cleaned = clean_svi_data(svi, selected_variables)

    # print(svi_cleaned.head(3))

    for variable in selected_variables:
        # print(variable)
        # Select only the required columns
        df_one_variable = svi_cleaned[["FIPS", variable, "geometry"]]

        # # Sorting the DataFrame based on the 'variable' column
        df_one_variable = df_one_variable.sort_values(by=variable)
        df_one_variable["sortedID"] = range(len(df_one_variable))

        # Convert the DataFrame to a GeoDataFrame
        df_one_variable = gpd.GeoDataFrame(df_one_variable, geometry="geometry")
        df_one_variable.crs = "EPSG:3395"  # This is a commonly used projected CRS

        # Generate adjacent counties # selecting only the counties below the given percentile happens here
        adjacencies_list, adjacent_counties_df, county_list = (
            generate_adjacent_counties(df_one_variable, variable, percentile)
        )

        # Create a dictionary adjacent_counties_df column county as key and column adjacent as value
        # (to avoid NULL adjacencies error)
        adjacent_counties_dict = dict(
            zip(adjacent_counties_df["county"], adjacent_counties_df["adjacent"])
        )

        # This take only counties that have adjacent counties
        county_list = adjacent_counties_df["county"].tolist()

        # returns the simplicial complex connected compontn,edges and triangles (as a list of lists-which
        # cointains the sorted ID of the counties)
        V = form_simplicial_complex(adjacent_counties_dict, county_list)

        # This is a new feature that I added to the code. It creates a new list replace
        # the sorted ID with the FIPS on the V list
        # create a new list replace the sorted ID with the FIPS on the V list
        V_FIPS = [[df_one_variable.iloc[x]["FIPS"] for x in i] for i in V]

        # Plot the simplicial complex # plot will be saved in the results folder
        # change the path to save the plot in the desired location in the plot_simplicial_complex function
        plot_simplicial_complex(df_one_variable, V_FIPS, variable, results_path)

        # Generate scaled marginal variance
        marginal_var_df, adjacency_matrix = generate_scaled_marginal_variance(
            V, df_one_variable, variable
        )

        with open(f"{results_path}/{state}_adj_{variable}.pkl", "wb") as f_out:
            pickle.dump(adjacency_matrix, f_out)

        # print the marginal variance
        print(f"Marginal variances for the variable {variable} \n")
        print(marginal_var_df)

        # plot the scaled marginal variance
        fig, axes = plt.subplots(figsize=(10, 10))

        df_one_variable.plot(
            column=variable + "_marginal_variance",
            cmap="OrRd",
            linewidth=0.5,
            ax=axes,
            edgecolor="black",
            legend=True,
        )
        df_one_variable.boundary.plot(ax=axes, color="black", linewidth=0.5)
        axes.set_title(
            variable + "_marginal_variance", fontsize=12
        )  # Adjust the fontsize here
        axes.set_axis_off()

        plt.tight_layout()
        plt.savefig(f"{results_path}/{variable}_marginal_variance.png")
        plt.close()

        # break

    print(f"All variables processed for the state {state}.")


# if __name__ == "__main__":
#     main()
