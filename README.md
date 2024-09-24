# ornl-scaled-marginal-variance
Generates scaled marginal variance for a particular region based on the variable that you are using.

Here's a basic `README.md` file to help others run your project:

---

# Scaled Marginal Variance Analysis

This repository contains code for analyzing scaled marginal variances for different counties based on Social Vulnerability Index (SVI) data. It also constructs and plots simplicial complexes from the adjacency relationships between counties.

## Getting Started

### Prerequisites

To run the code, you need the following Python libraries:

- `numpy`
- `pandas`
- `geopandas`
- `matplotlib`
- `tqdm`
- `scipy`
- `pysal`
- `invr` (custom library/module included)
  
You can install the necessary libraries using pip:

```bash
pip install numpy pandas geopandas matplotlib tqdm scipy pysal invr
```

Additionally, make sure you have the Social Vulnerability Index (SVI) data in ESRI Geodatabase format for the specific state you are analyzing. This data can be downloaded from the [CDC SVI site](https://www.atsdr.cdc.gov/placeandhealth/svi/data_documentation_download.html).

### Directory Structure

```bash
.
├── data/
│   └── SVI2018_WYOMING_COUNTY.gdb   # The SVI data for Wyoming
├── results/                         # Directory where output plots and results will be saved
├── code/
│   └── main.py                      # Main script for the analysis
└── README.md                        # This file
```

### Running the Code

1. Download and place the SVI geodatabase file in the `data/` directory.
   
2. In the `code/main.py`, update the following paths if necessary:
   - `svi_gdb_path`: Path to the SVI geodatabase file.
   - `results_path`: Path where results and plots will be saved.
   
3. Execute the main script by running:

```bash
python code/main.py
```

### Output

The script will process the SVI data for the specified state and variable, generating:
- **Simplicial Complex Plots**: Visual representations of the adjacency relationships between counties.
- **Marginal Variance Plots**: Plots showing the marginal variance for each variable.

The outputs will be saved in the `results/` folder.

### Customization

You can modify the following variables in the `main.py` script:
- `selected_variables`: The list of SVI variables you want to analyze.
- `percentile`: The threshold to filter the counties for analysis.

### Example

```bash
Processing variables: 100%|██████████| 15/15 [02:15<00:00,  9.03s/it]
Marginal variances for the variable EP_POV:
    FIPS    EP_POV  EP_POV_marginal_variance
0  56001  25.6                0.0213
1  56003  18.9                0.0357
```

### License

This project is licensed under the MIT License.

---

This should provide clear instructions for running your code and understanding its functionality. Let me know if you need further adjustments!
