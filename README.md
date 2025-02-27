# Elastic-Net-Cross-Validation
Functions that perform elastic net modeling using given features and covariates

**R Functions**   
en_cv_test_fast(clin_df, protein_list, cat_control, cont_control, trait_list, alpha, num_ensp, composite=FALSE, direction="up")
en_repeat_fast(clin_df, protein_list, cat_control, cont_control, trait_list, alpha, heatmap=FALSE)


## Overview  
Elastic net is a form of penalized linear regression that performs feature selection by shrinking beta coefficients of predictor variables using a penalty term, 
the stength of which is determined using a hyperperamater lambda, with a smaller lambda selecting more predictor terms, and larger selecting fewer. A range of values are used for lambda for cross validation, and ultimately the model using the lambda value and selected features that produces the lowest mean squared error (MSE) are found. In this way, Elastic Net can be used to identify important features, for example gene expression levels, implicated in clinical parameters (age, bmi, etc.).

One downside of elastic net modeling using common R functions such as cv.glmnet is the variability of the output, as folds for cross validation are selected randomly. This function ensures reproducibility during elastic net modeling by setting explicit fold IDs, performs cross validation reprodicibly using the caret package, auto-scales features, uses parallel processing for faster output, and produces a table indicating out-of-sample predictive ability of each set of features and for each parameter, as well as a list of selected features per parameter.


## Usage  

**Arguments**  
| Parameter       | Type        | Description                                                                                                    |
|----------------|-------------|-----------------------------------------------------------------------------------------------------------------|
| `clin_df `     | `dataframe` | Table of microarray data containing gene expression values and covariates by columns, sample IDs by row         |
| `protein_list` | `vector `   | Vector of genes from which to perform feature selection, and any continuous covariates                          |
| `cat_control`  | `vector `   | Vector of categorical covariates                                                                                |
| `cont_control` | `vector `   | Vector of continuous covariates                                                                                 |
| `trait_list`   | `list `     | Vector of continuous or categorical traits for which to select implicated features                              |
| `alpha`        | `numeric `  | Number indicating hyperparameter alpha (0 for ridge, 1 for lasso, in-between for Elastic Net)                   |
| `num_ensp`     | `numeric `  | Number indicating a limit of features selected per paramter                                                     |
| `composite`    | `boolean `  | Boolean value determining selected features will be condensed into a mean composite score                       |
| `direction`    | `string `   | String "up" or "down" indicating if features positively or negatively associated with parameters will be used   |

## Return Values  
- Coef is a table showing the penalized effect size of each feature for each supplied trait.
- Lambda is a table showing the supplied alpha and optimized lambda across many runs for each parameter of the trait_list.
- Heatmap is a graphical representation of the penalized effect size of each feature for each supplied trait.

<p align="center">
  <img src="images/Example_heatmap2.JPG" alt="Example Image of Selected Features" width="500">
</p>

- IVSum is a bar graph respresenting the number of selected features for each supplied parameter. Example image shows 100 possible features and 2 covariates.

<p align="center">
  <img src="images/Example_ivsum.JPG" alt="Example Image of Selected Features" width="500">
</p>

## Dependencies  
glmnet  
ggplot2  
heatmaply  
doParallel

## Notes  
This function will produce an error if the heatmap functionality is set to true and any of the parameters have 0 selected features. If this error occurs, set heatmap to false.

## Author  
Bradley Olinger, PhD  
b.a.olinger@gmail.com
