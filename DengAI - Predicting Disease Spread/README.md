# DengAI - Predicting Disease Spread using XGRegressor

This project is a challenge from [Driven Data](https://www.drivendata.org/competitions/44/dengai-predicting-disease-spread/page/80/).

The current status: **Finished**

## Overview

### Can you predict local epidemics of dengue fever?

Dengue fever is a mosquito-borne disease that occurs in tropical and sub-tropical parts of the world. In mild cases, symptoms are similar to the flu: fever, rash, and muscle and joint pain. In severe cases, dengue fever can cause severe bleeding, low blood pressure, and even death.

Because it is carried by mosquitoes, the transmission dynamics of dengue are related to climate variables such as temperature and precipitation. Although the relationship to climate is complex, a growing number of scientists argue that climate change is likely to produce distributional shifts that will have significant public health implications worldwide.

Using environmental data collected by various U.S. Federal Government agencies—from the Centers for Disease Control and Prevention to the National Oceanic and Atmospheric Administration in the U.S. Department of Commerce— *can you predict the number of dengue fever cases reported each week in San Juan, Puerto Rico and Iquitos, Peru?*

### Goal
The project's goal is to predict the **total_cases** label for each (**city, year, weekofyear**) in the test set. There are two cities, San Juan and Iquitos, with test data for each city spanning 5 and 3 years respectively. You will make one submission that contains predictions for both cities. The data for each city have been concatenated along with a city column indicating the source: **sj** for San Juan and **iq** for Iquitos. The test set is a pure future hold-out, meaning the test data are sequential and non-overlapping with any of the training data. Throughout, missing values have been filled as *NaNs*.

Click [here](https://github.com/bfroz/Projects/blob/master/DengAI%20-%20Predicting%20Disease%20Spread/DengAI%20-%20Predicting%20Disease%20Spread.ipynb) to check the full project.

### Results

The Mean Absolute Error (MAE) for the whole data is:

|CITY|	MAE		|
|-|-----------------------|
| SJ| 20.017 |
| IQ  | 5.468 |