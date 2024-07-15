# Bayesian Dynamic Borrowing Project Repository

## Project Description

Our collaborative project with Procter & Gamble (P&G) explores the application of Bayesian Dynamic Borrowing (BDB) methodologies to optimize Randomized Controlled Trials (RCTs) in oral care and face cream datasets using historical data integration. The primary aim is to reduce sample size requirements for control groups while maintaining statistical rigor. We focus on developing adaptive strategies that incorporate historical data into current trials, thereby enhancing efficiency and accelerating decision-making in clinical studies.

We evaluate various BDB techniques such as power priors, commensurate priors, Elastic Priors, Bayesian hierarchical models, and meta-analytic-predictive priors. Through simulations and analysis, we assess these methods’ performance in oral care and face cream datasets, particularly examining bias, variance, and statistical power across different scenarios of treatment and historical data alignment.

Initial findings suggest that methods like the power prior and commensurate prior effectively control bias, while the Elastic Prior provides a balanced approach considering bias, variance, and computational efficiency, especially in scenarios of partial congruence between historical and current trial data.

Furthermore, we demonstrate the application of these methodologies through sample analyses, such as varying the means of treatment and control groups, and historical versus control groups, exploring the impact of individual noise, and assessing how historical data influences treatment distributions in a Round-robin model. This research not only highlights the practical advantages of advanced Bayesian techniques in oral care and cosmetic product trials but also contributes to advancing statistical methodology in clinical trial design and analysis.

## Repository Structure

- `Face_Cream` Directory
    - `commensurate.R`: Script implementing commensurate prior methodology for the face cream dataset.
    - `elastic.R`: Script implementing elastic prior methodology for the face cream dataset.
    - `normalized.R`: Script implementing normalized prior methodology for the face cream dataset.
    - `power.R`: Script implementing power prior methodology for the face cream dataset.
    - `results` Directory
        - `elastic_results_nc_fc_temp.csv`: Contains results of the elastic prior simulations for the face cream dataset.

- `Oral_Care` Directory
    - `Plot_Result.R`: Script for plotting results of the analyses.
    - `commensurate.R`: Script implementing commensurate prior methodology for the oral care dataset.
    - `elastic.R`: Script implementing elastic prior methodology for the oral care dataset.
    - `elastic_power.R`: Script implementing elastic power prior methodology for the oral care dataset.
    - `normalized.R`: Script implementing normalized prior methodology for the oral care dataset.
    - `robust_map.R`: Script implementing robust meta-analytic-predictive prior methodology for the oral care dataset.
    - `results` Directory
        - `commensurate_results_nc.csv`: Contains results of the commensurate prior simulations for the oral care dataset.
        - `elastic_power_results_nc.csv`: Contains results of the elastic power prior simulations for the oral care dataset.
        - `elastic_results_nc.csv`: Contains results of the elastic prior simulations for the oral care dataset.
        - `normalized_results_nc.csv`: Contains results of the normalized prior simulations for the oral care dataset.
        - `rMAP_results_nc.csv`: Contains results of the robust meta-analytic-predictive prior simulations for the oral care dataset.

- `app` Directory
    - `app.R`: Main script to run the Shiny app.
    - `helper.R`: Helper functions for the Shiny app.
    - `helper_fc.R`: Helper functions specific to the face cream dataset for the Shiny app.
    - `helper_oc.R`: Helper functions specific to the oral care dataset for the Shiny app.
    - `results_fc` Directory
        - `elastic_results_nc_fc.csv`: Contains elastic prior results for the face cream dataset used in the Shiny app.
        - `elastic_results_tau1_updated.csv`: Contains updated elastic prior results for the face cream dataset with tau1 parameter.
        - `elastic_results_tau2_updated.csv`: Contains updated elastic prior results for the face cream dataset with tau2 parameter.
    - `results_oc` Directory
        - `commensurate_results_nc.csv`: Contains commensurate prior results for the oral care dataset used in the Shiny app.
        - `elastic_power_results_nc.csv`: Contains elastic power prior results for the oral care dataset used in the Shiny app.
        - `elastic_results_nc.csv`: Contains elastic prior results for the oral care dataset used in the Shiny app.
        - `normalized_results_nc.csv`: Contains normalized prior results for the oral care dataset used in the Shiny app.
        - `rMAP_results_nc.csv`: Contains robust meta-analytic-predictive prior results for the oral care dataset used in the Shiny app.
    - `rsconnect/shinyapps.io/3vxc23-naman-agrawal/RIPS_app.dcf`: Configuration file for deploying the Shiny app to shinyapps.io.
    - `sample_generated_data.xlsx`: Sample generated data file for testing the Shiny app.

- `old` Directory
    - Contains older versions of files which may be referred to for historical context or past attempts.

- `requirements.txt`
    - Contains the list of required R packages and their versions to ensure the code runs smoothly. Install the dependencies by running:
        ```bash
        Rscript -e 'install.packages(readLines("requirements.txt"))'
        ```

- `README.md`
    - This file. Provides documentation for the repository.

## R Version
The project was developed using R version 4.1.0. Ensure that you have this version or a compatible one installed to avoid any compatibility issues.

## Running the Shiny App
To run the Shiny app, follow these steps:

1. Install Required Packages: Ensure all required packages are installed. Run the following command in your R console:
     ```bash
     Rscript -e 'install.packages(readLines("requirements.txt"))'
     ```

2. Navigate to the App Directory: Change your working directory to the app directory where the Shiny app files are located.
     ```bash
     setwd("path/to/app")
     ```

3. Run the App: Execute the following command in your R console to start the Shiny app.
     ```bash
     shiny::runApp()
     ```

4. Access the App: Open a web browser and navigate to the provided local URL (typically http://127.0.0.1:XXXX) to access the Shiny app interface.


## Conclusion

This repository contains a comprehensive set of scripts and data for analyzing Bayesian Dynamic Borrowing methodologies applied to oral care and face cream datasets. By following the provided documentation, you can replicate the analyses and explore the results through the Shiny app.


This markdown file provides a clear and detailed documentation for the repository, including the project description, file descriptions, instructions for running the Shiny app, and other relevant details. Let us know if you need any modifications or additional information!

