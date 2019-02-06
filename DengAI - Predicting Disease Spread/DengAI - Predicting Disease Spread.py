
# coding: utf-8

# # DengAI - Predicting Disease Spread
# 
# This project is a challenge from [Driven Data](https://www.drivendata.org/competitions/44/dengai-predicting-disease-spread/page/80/).
# 
# The current status: **Ongoing...**
# 
# ## Overview
# 
# ### Can you predict local epidemics of dengue fever?
# 
# Dengue fever is a mosquito-borne disease that occurs in tropical and sub-tropical parts of the world. In mild cases, symptoms are similar to the flu: fever, rash, and muscle and joint pain. In severe cases, dengue fever can cause severe bleeding, low blood pressure, and even death.
# 
# Because it is carried by mosquitoes, the transmission dynamics of dengue are related to climate variables such as temperature and precipitation. Although the relationship to climate is complex, a growing number of scientists argue that climate change is likely to produce distributional shifts that will have significant public health implications worldwide.
# 
# Using environmental data collected by various U.S. Federal Government agencies—from the Centers for Disease Control and Prevention to the National Oceanic and Atmospheric Administration in the U.S. Department of Commerce— *can you predict the number of dengue fever cases reported each week in San Juan, Puerto Rico and Iquitos, Peru?*
# 
# ### Goal
# The project's goal is to predict the **total_cases** label for each (**city, year, weekofyear**) in the test set. There are two cities, San Juan and Iquitos, with test data for each city spanning 5 and 3 years respectively. You will make one submission that contains predictions for both cities. The data for each city have been concatenated along with a city column indicating the source: **sj** for San Juan and **iq** for Iquitos. The test set is a pure future hold-out, meaning the test data are sequential and non-overlapping with any of the training data. Throughout, missing values have been filled as *NaNs*.

# ## Features description
# 
# The description of every feature is following:
# ### City and date indicators
# - **city** – City abbreviations: sj for San Juan and iq for Iquitos
# - **week_start_date** – Date given in yyyy-mm-dd format
# 
# ### NOAA's GHCN daily climate data weather station measurements
# - **station_max_temp_c** – Maximum temperature
# - **station_min_temp_c** – Minimum temperature
# - **station_avg_temp_c** – Average temperature
# - **station_precip_mm** – Total precipitation
# - **station_diur_temp_rng_c** – Diurnal temperature range
# 
# ### PERSIANN satellite precipitation measurements (0.25x0.25 degree scale)
# - **precipitation_amt_mm** – Total precipitation
# 
# ### NOAA's NCEP Climate Forecast System Reanalysis measurements (0.5x0.5 degree scale)
# - **reanalysis_sat_precip_amt_mm** – Total precipitation
# - **reanalysis_dew_point_temp_k** – Mean dew point temperature
# - **reanalysis_air_temp_k** – Mean air temperature
# - **reanalysis_relative_humidity_percent** – Mean relative humidity
# - **reanalysis_specific_humidity_g_per_kg** – Mean specific humidity
# - **reanalysis_precip_amt_kg_per_m2** – Total precipitation
# - **reanalysis_max_air_temp_k** – Maximum air temperature
# - **reanalysis_min_air_temp_k** – Minimum air temperature
# - **reanalysis_avg_temp_k** – Average air temperature
# - **reanalysis_tdtr_k** – Diurnal temperature range
# 
# ### Satellite vegetation - Normalized difference vegetation index (NDVI) - NOAA's CDR Normalized Difference Vegetation Index (0.5x0.5 degree scale) measurements
# - **ndvi_se** – Pixel southeast of city centroid
# - **ndvi_sw** – Pixel southwest of city centroid
# - **ndvi_ne** – Pixel northeast of city centroid
# - **ndvi_nw** – Pixel northwest of city centroid

# ## Understanding the problem, understanding the features

# Before any other step, it's really important to understand the disease and the given features, then we can get insights from it.

# ### What is Dengue and how is it proliferated?
# 
# Dengue is a mosquito transmitted disease. The mosquito Aedes aegypti, a primary vector of dengue, yellow fever, and chikungunya viruses, is widely distributed in the **subtropics and tropics** [[1]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3516267/). According to WHO (*World Health Organization*), **climate change may affect transmission**, as dengue mosquitoes reproduce more quickly and bite more frequently at **higher temperatures** [[2]](https://www.who.int/heli/risks/vectors/denguecontrol/en/). Also, in proximity to human settlements, Aedes aegypti mosquitoes breed primarily in **artificial water containers**, and the mosquito’s life-cycle is closely associated with **human activities**. Larval habitats are increasing rapidly in **urban areas**.
# 
# The temperature variations and **rainfall intensity** affect the reproductive cycle and survival of the vector, which cause changes in its distribution and density, since mosquitoes need **humidity** and temperatures ranging between **15°C and 35°C** to survive and reproduce [[3]](http://www.scielo.br/pdf/ramb/v63n11/0104-4230-ramb-63-11-0957.pdf). Rainfall provides plenty of breeding sites for mosquito vectors such as puddles, while humidity affects the adult mosquitoes’ survival and **biting frequency** [[4]](https://bmcpublichealth.biomedcentral.com/articles/10.1186/s12889-018-5532-4). Finally, there is a **negative correlation** between Diurnal Temperature Range (**DTR**) and dengue incident. [[5]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5002035/)

# ### Features explained
# 
# The driven features in the dataset is **city**, **year**, **weekofyear**. Also, **week_start_date** can be helpful to understand trends and patterns in data.
# 
# Some explanatory variables are very intuitive. The precipitation features represents the amount (in mm) of water downfall on that current week (**station_precip_mm**, **precipitation_amt_mm**, **reanalysis_sat_precip_amt_mm**, **reanalysis_precip_amt_kg_per_m2**). The same precipitation information are given from different sources and different scales.
# 
# The temperature features represents the weekly climate conditions. This information are given in a tricky way. For example, **station_max_temp_c** is the weekly maximum temperature average, and **station_min_temp_c** weekly minimum average. But **station_avg_temp_c**, which is the average temperature, is not the mean between the weekly maximum average and minimum average, but the mean of the daily average from that week. The same goes to the air temperature features (**reanalysis_air_temp_k**, **reanalysis_max_air_temp_k**, **reanalysis_min_air_temp_k**, **reanalysis_avg_temp_k**). Its important to notice that some temperature features are expressed in **Kelvin** and others in **Celsius**.
# 
# The humidity are measured in several ways in this dataset. There is a Dew Point temperature (**eanalysis_dew_point_temp_k**), relative humidity (**reanalysis_relative_humidity_percent**) and specific humidity (**reanalysis_specific_humidity_g_per_kg**). Some studies confirms that dew point and humidity have a nearly linear relationship [[6]](https://journals.ametsoc.org/doi/10.1175/BAMS-86-2-225). It is possible to measure the humidity percentage using the dew point and air temperature [[7]](http://www.reahvac.com/tools/humidity-formulas/).
# 
# Finally, the **NDVI** is a measure of how much green vegetation is presented at a determined location [[8]](https://gisgeography.com/ndvi-normalized-difference-vegetation-index/). This number ranges from -1 (none or less green) to +1 (very green).

# We can note that the dengue occurrence may be caused by a series of factors:
# - High temperatures
# - Rainfall intensity
# - Humidity
# - Urban areas and human proximity
# 
# Doing a quick recap on the features we have, we can group them related to the possible dengue spread causes.
# 
# |	     Temperature	  |	           Rain             |               Humidity              | Urbanization |
# |-------------------------|-----------------------------|-------------------------------------|--------------|
# |    station_max_temp_c   |   station_precip_mm         |    reanalysis_dew_point_temp_k      |   ndvi_se    |
# |    station_min_temp_c   |   precipitation_amt_mm      |reanalysis_relative_humidity_percent |   ndvi_sw    |
# |    station_avg_temp_c   |reanalysis_sat_precip_amt_mm |reanalysis_specific_humidity_g_per_kg|   ndvi_ne    |
# | station_diur_temp_rng_c |reanalysis_precip_amt_kg_per_m2|        reanalysis_air_temp_k        |   ndvi_nw    |
# |    reanalysis_tdtr_k    |             -               |      reanalysis_max_air_temp_k      |      -       |
# |   reanalysis_avg_temp_k |             -               |      reanalysis_min_air_temp_k      |      -       |
# 
# Important to remember: the target variable is **'total_cases'**. We have to predict the quantity of cases in **two different cities**. In a given moment, we will separate the data in two dataframes and observing the behavior separated.

# ## Importing data
# 
# Let's begin the data manipulation based on what we learning so far.

# In[2]:


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error


# The data is separated in two datasets: features_train and labels.

# In[3]:


df_features_train = pd.read_csv('data/dengue_features_train.csv')
df_labels_train = pd.read_csv('data/dengue_labels_train.csv')


# In[4]:


df_features_train.info()


# In[5]:


df_labels_train.info()


# We can observe that the labels contain the **dependent variable**, our target, **total_cases**. To facilitate the data manipulation, let's put all data into a single *dataframe*.

# In[6]:


df_train = pd.merge(df_features_train, df_labels_train, how='left', left_on=['city','year','weekofyear'], right_on=['city','year','weekofyear']) 


# In[7]:


df_train.head()


# ### Treating missing values
# 
# The first thing to do with the data is checking missing values. We can use a heatmap from seaborn to visualize the missing value relevance before treating data.

# In[8]:


fig, ax = plt.subplots(figsize=(20,11))
sns.heatmap(df_train.isnull(),yticklabels=False,cbar=False,cmap='viridis', ax=ax)


# It is very clear some missing values in this dataset, represented by the yellow color in the heatmap. Some data are **completely missing**, except for the labels. Also is notable the missing values in the **ndvi** features, some **temperature** features from *NOAA's GHCN* source, and some in the **precipitation** features.
# 
# There are some ways to overcome this issue.

# #### Precipitation features
# 
# This dataset contains three precipitation amount features representing similar the same information, but from different sources. Let's see how similar those features are.

# In[9]:


cols_precipitation = ['station_precip_mm', 'precipitation_amt_mm', 'reanalysis_sat_precip_amt_mm', 'reanalysis_precip_amt_kg_per_m2']


# In[10]:


pd.plotting.scatter_matrix(df_train[cols_precipitation], alpha = 0.3, figsize = (14,9), diagonal = 'kde');


# There is a **multicollinearity** between the **precipitation_amt_mm**, **reanalysis_sat_precip_amt_mm** features. Also, we can observe that the **station_precip_mm** and **reanalysis_precip_amt_kg_per_m2** features don't present a linear relation with the other features. This is *plausible*, since those features were generated with different techniques and from different sources.
# > Preciptation measured in kg_per_m2 and mm are equivalent.

# In[11]:


# Filtering just missing values from precipitation_amt_mm and reanalysis_sat_precip_amt_mm, where exist value on station_precip_mm column
df_train[cols_precipitation][~df_train['station_precip_mm'].isnull() & (df_train['precipitation_amt_mm'].isnull() | df_train['reanalysis_sat_precip_amt_mm'].isnull())]


# Observe that all missing values from **precipitation_amt_mm** are presented in the column **reanalysis_sat_precip_amt_mm**. So, it is safe to **remove *any*** of those two features. Just for convention, the **reanalysis_sat_precip_amt_mm** will be removed later, in the appropriate moment.
# 
# Notice that some missing values from those two features are presented on **station_precip_mm**. It is possible to fulfill those missing values just copying from the **station_precip_mm** or **reanalysis_precip_amt_kg_per_m2** features. The inverse is also possible. It is also possible to use the *mean* between the others features to fulfill those values.

# In[12]:


cols_precipitation = ['station_precip_mm', 'precipitation_amt_mm', 'reanalysis_precip_amt_kg_per_m2']


# In[13]:


# Showing just the number of missing values from station_precip_mm that contains value on precipitation_amt_mm or reanalysis_precip_amt_kg_per_m2
df_train[cols_precipitation][(~df_train['precipitation_amt_mm'].isnull() | ~df_train['reanalysis_precip_amt_kg_per_m2'].isnull()) & df_train['station_precip_mm'].isnull()].shape                            


# In[14]:


#Checking one of the NaN feature
df_train[cols_precipitation].iloc[[94]]


# In[15]:


def fill_nan_average(row, cols):
    return row[cols].mean()


# In[16]:


#fill NaN values of precipitation_amt_mm
df_train['precipitation_amt_mm'] = df_train.apply(lambda row: 
                                               fill_nan_average(row, ['station_precip_mm', 'reanalysis_precip_amt_kg_per_m2']) 
                                               if np.isnan(row['precipitation_amt_mm']) & (~np.isnan(row['station_precip_mm']) | ~np.isnan(row['reanalysis_precip_amt_kg_per_m2']))
                                               else row['precipitation_amt_mm'], axis=1)


# In[17]:


#Checking one of the NaN again, to see the change
df_train[cols_precipitation].iloc[[94]]


# In[18]:


#Let's fill the others precipitation features
df_train['station_precip_mm'] = df_train.apply(lambda row: 
                                               fill_nan_average(row, ['precipitation_amt_mm', 'reanalysis_precip_amt_kg_per_m2']) 
                                               if np.isnan(row['station_precip_mm']) & (~np.isnan(row['precipitation_amt_mm']) | ~np.isnan(row['reanalysis_precip_amt_kg_per_m2']))
                                               else row['station_precip_mm'], axis=1)
df_train['reanalysis_precip_amt_kg_per_m2'] = df_train.apply(lambda row: 
                                               fill_nan_average(row, ['station_precip_mm', 'precipitation_amt_mm']) 
                                               if np.isnan(row['reanalysis_precip_amt_kg_per_m2']) & (~np.isnan(row['station_precip_mm']) | ~np.isnan(row['precipitation_amt_mm']))
                                               else row['reanalysis_precip_amt_kg_per_m2'], axis=1)


# #### NDVI features
# 
# The NDVI features represents the vegetation around the city corners. Instead of use those features separated, we can derivate them, conserving the relevance and reducing dimensionality. For example, we can remove the four corners features and substitute for a *mean* value.

# In[19]:


cols_ndvi = ['ndvi_ne', 'ndvi_nw', 'ndvi_se', 'ndvi_sw']


# In[20]:


def ndvi_mean(row):
    return row[cols_ndvi].mean()


# In[21]:


#Creating new features
df_train['ndvi_mean'] = df_train.apply(lambda row: ndvi_mean(row), axis=1)


# In[22]:


#Removing the old NDVI features
df_train.drop(cols_ndvi, axis=1, inplace=True)


# In[23]:


cols_ndvi = ['ndvi_mean']


# #### Temperature features
# 
# Before check any temperature feature, even the ones regard humidity, it is notable that some are given in *Kelvin*. It is necessary to normalize to a unique measure. I will foster *Celsius* as the standard unit.
# 
# > The only exception is the **reanalysis_tdtr_k** feature, since diurnal temperature is independent of unit.

# In[24]:


def convertKelvinToCelsius(k):
    return k - 273.15


# In[25]:


cols_kelvin = ['reanalysis_avg_temp_k', 'reanalysis_dew_point_temp_k', 'reanalysis_air_temp_k', 'reanalysis_max_air_temp_k', 'reanalysis_min_air_temp_k']


# In[26]:


#Creating new Celsius features
df_train['reanalysis_avg_temp_c'] = df_train.apply(lambda row: convertKelvinToCelsius(row['reanalysis_avg_temp_k']), axis=1)
df_train['reanalysis_dew_point_temp_c'] = df_train.apply(lambda row: convertKelvinToCelsius(row['reanalysis_dew_point_temp_k']), axis=1)
df_train['reanalysis_air_temp_c'] = df_train.apply(lambda row: convertKelvinToCelsius(row['reanalysis_air_temp_k']), axis=1)
df_train['reanalysis_max_air_temp_c'] = df_train.apply(lambda row: convertKelvinToCelsius(row['reanalysis_max_air_temp_k']), axis=1)
df_train['reanalysis_min_air_temp_c'] = df_train.apply(lambda row: convertKelvinToCelsius(row['reanalysis_min_air_temp_k']), axis=1)


# In[27]:


#Removing Kelvin features
df_train.drop(cols_kelvin, axis=1, inplace=True)


# In[28]:


cols_temperature = ['station_max_temp_c', 'station_min_temp_c', 'station_avg_temp_c', 'station_diur_temp_rng_c', 'reanalysis_tdtr_k', 'reanalysis_avg_temp_c']


# Let's check again the missing values of the temperature features. Then, we check the correlation.

# In[29]:


fig, ax = plt.subplots(figsize=(12,5))
sns.heatmap(df_train[cols_temperature].isnull(),yticklabels=False,cbar=False,cmap='viridis', ax=ax)


# In[30]:


pd.plotting.scatter_matrix(df_train[cols_temperature], alpha = 0.3, figsize = (14,12), diagonal = 'kde');


# It's notable that there is a high linear correlation between the **station_avg_temp_c** and **reanalysis_avg_temp_c**. Also, the **station_avg_temp_c** contains some missing values. So, it is *safe* to use the **reanalysis_avg_temp_c** to fulfill the missing values from **station_avg_temp_c**. The same applies between the **station_diur_temp_rng_c** and **reanalysis_tdtr_k**.

# In[31]:


df_train['station_avg_temp_c'] = df_train.apply(lambda row: 
                                               row['reanalysis_avg_temp_c'] if np.isnan(row['station_avg_temp_c'])
                                               else row['station_avg_temp_c'], axis=1)
df_train['station_diur_temp_rng_c'] = df_train.apply(lambda row: 
                                               row['reanalysis_tdtr_k'] if np.isnan(row['station_diur_temp_rng_c'])
                                               else row['station_diur_temp_rng_c'], axis=1)


# > For now, we do not fulfill the min and max temperature feature with missing values

# #### The other features
# 
# Let's check now how the dataset is after those fulfilling and creating features.

# In[32]:


fig, ax = plt.subplots(figsize=(20,11))
sns.heatmap(df_train.isnull(),yticklabels=False,cbar=False,cmap='viridis', ax=ax)


# We can see that still have some features to fulfill. Those columns are *harder* to work with, because they don't have dependency with others. However, we have a dataset with an advantage: a **time line**.
# 
# With this powerful resource, we can **interpolate** the missing values, in order to remove NaN values.

# #### Interpolating missing values
# 
# Before starting the interpolation, we have to remember: there are two different cities. That means, there are two different *timelines*. So, to interpolate values, we should first split the dataset in two.

# In[33]:


#Creating week_start_date as index
df_train['week_start_date'] = pd.to_datetime(df_train['week_start_date'], format='%Y-%m-%d')
df_train.set_index('week_start_date', inplace=True)
#sorting index
df_train.sort_index(inplace=True)


# In[34]:


#Checking the cities values
df_train.city.unique()


# In[35]:


df_train_sj = df_train[df_train['city'] == 'sj']
df_train_iq = df_train[df_train['city'] == 'iq']
#Dropping city column
df_train_sj = df_train_sj.drop(['city'], axis=1)
df_train_iq = df_train_iq.drop(['city'], axis=1)


# To interpolate values to fulfill NaN, we use interpolate() from pandas. First we have to index and sort the dataset, using the **week_start_date** feature.
# 
# To be more clear, give a look at a previous known missing value and their neighbors in the sj dataframe.

# In[36]:


df_train_sj[['station_precip_mm']].iloc[77:97].plot()


# There is a missing value in this example. With interpolate method, we can find a value between all missing values from this dataset.

# In[37]:


df_train_sj_inter = df_train_sj.interpolate()
df_train_iq_inter = df_train_iq.interpolate()


# With the new dataframe, with interpolated values, it is possible to check the fulfilled values.

# In[38]:


df_train_sj_inter[['station_precip_mm']].iloc[77:97].plot()


# So, all missing values are fulfilled, without lose much information.

# In[39]:


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(20, 6))
sns.heatmap(df_train_sj_inter.isnull(),yticklabels=False,cbar=False,cmap='viridis', ax=ax[0])
sns.heatmap(df_train_iq_inter.isnull(),yticklabels=False,cbar=False,cmap='viridis', ax=ax[1])


# TODO:
# - Add cycle features and remove year and weekofyear
# - Test regression with presented data
# - Chose most relevant features
# - Time Series those features
# - Predict total_cases
# - Validate model performance

# #### Removing Multicollinearity

# In[40]:


df_train_sj_inter['week_sin'] = np.sin((df_train_sj_inter.weekofyear-1)*(2.*np.pi/53))
df_train_sj_inter['week_cos'] = np.cos((df_train_sj_inter.weekofyear-1)*(2.*np.pi/53))
df_train_iq_inter['week_sin'] = np.sin((df_train_iq_inter.weekofyear-1)*(2.*np.pi/53))
df_train_iq_inter['week_cos'] = np.cos((df_train_iq_inter.weekofyear-1)*(2.*np.pi/53))

df_train_sj_inter = df_train_sj_inter.drop(['weekofyear'], axis=1)
df_train_iq_inter = df_train_iq_inter.drop(['weekofyear'], axis=1)


# In[41]:


df_train_sj_inter = df_train_sj_inter.drop(['year'], axis=1)
df_train_iq_inter = df_train_iq_inter.drop(['year'], axis=1)


# In[42]:


new_cols = ['precipitation_amt_mm','reanalysis_precip_amt_kg_per_m2','reanalysis_relative_humidity_percent','reanalysis_specific_humidity_g_per_kg','reanalysis_tdtr_k','station_avg_temp_c','station_diur_temp_rng_c','station_max_temp_c','station_min_temp_c','station_precip_mm','ndvi_mean','reanalysis_avg_temp_c','reanalysis_dew_point_temp_c','reanalysis_air_temp_c','reanalysis_max_air_temp_c','reanalysis_min_air_temp_c']
new_cols_with_cyclical_features = ['week_sin','week_cos','precipitation_amt_mm','reanalysis_precip_amt_kg_per_m2','reanalysis_relative_humidity_percent','reanalysis_specific_humidity_g_per_kg','reanalysis_tdtr_k','station_avg_temp_c','station_diur_temp_rng_c','station_max_temp_c','station_min_temp_c','station_precip_mm','ndvi_mean','reanalysis_avg_temp_c','reanalysis_dew_point_temp_c','reanalysis_air_temp_c','reanalysis_max_air_temp_c','reanalysis_min_air_temp_c']


# In[43]:


from statsmodels.stats.outliers_influence import variance_inflation_factor

def calculate_vif_(X, thresh=100):
    cols = X.columns
    variables = np.arange(X.shape[1])
    dropped=True
    while dropped:
        dropped=False
        c = X[cols[variables]].values
        vif = [variance_inflation_factor(c, ix) for ix in np.arange(c.shape[1])]

        maxloc = vif.index(max(vif))
        if max(vif) > thresh:
            print('dropping \'' + X[cols[variables]].columns[maxloc] + '\' at index: ' + str(maxloc))
            variables = np.delete(variables, maxloc)
            dropped=True

    print('Remaining variables:')
    print(X.columns[variables])
    return X.columns[variables]


# > Code extracted from [this](https://stats.stackexchange.com/questions/155028/how-to-systematically-remove-collinear-variables-in-python) StackOverflow question.

# In[44]:


remaining_features = calculate_vif_(df_train_sj_inter[new_cols_with_cyclical_features])


# In[45]:


df_train_sj_inter_no_multi = df_train_sj_inter[remaining_features]
df_train_iq_inter_no_multi = df_train_iq_inter[remaining_features]


# In[46]:


df_train_sj_inter.info()


# In[47]:


df_train_sj_inter_no_multi.info()


# LINEAR REGRESSION TESTS BELOW:

# In[48]:


from sklearn.linear_model import LinearRegression
from sklearn import preprocessing


# In[49]:


#Round predictions and thresh to 0 and positives
def correct_prediction(pred):
    arr = np.round(pred).astype(int)
    arr[arr < 0] = 0
    return arr


# In[50]:


X = df_train_sj_inter_no_multi[['precipitation_amt_mm','reanalysis_precip_amt_kg_per_m2','reanalysis_tdtr_k','station_diur_temp_rng_c','station_min_temp_c','station_precip_mm','ndvi_mean']]
y = df_train_sj_inter[['total_cases']]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
regression_model = LinearRegression()
regression_model.fit(X_train, y_train)

y_predict = correct_prediction(regression_model.predict(X_test))
regression_model_mae = mean_absolute_error(y_predict, y_test)
print('sj_no_cycle_no_multi: ' + str(regression_model_mae))


# In[51]:


X = df_train_iq_inter_no_multi[['precipitation_amt_mm','reanalysis_precip_amt_kg_per_m2','reanalysis_tdtr_k','station_diur_temp_rng_c','station_min_temp_c','station_precip_mm','ndvi_mean']]
y = df_train_iq_inter[['total_cases']]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
regression_model = LinearRegression()
regression_model.fit(X_train, y_train)

y_predict = correct_prediction(regression_model.predict(X_test))
regression_model_mae = mean_absolute_error(y_predict, y_test)
print('iq_no_cycle_no_multi: ' + str(regression_model_mae))


# In[52]:


X = df_train_sj_inter[new_cols]
y = df_train_sj_inter[['total_cases']]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
regression_model = LinearRegression()
regression_model.fit(X_train, y_train)

y_predict = correct_prediction(regression_model.predict(X_test))
regression_model_mae = mean_absolute_error(y_predict, y_test)
print('sj_no_cycle: ' + str(regression_model_mae))


# In[53]:


X = df_train_iq_inter[new_cols]
y = df_train_iq_inter[['total_cases']]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
regression_model = LinearRegression()
regression_model.fit(X_train, y_train)

y_predict = correct_prediction(regression_model.predict(X_test))
regression_model_mae = mean_absolute_error(y_predict, y_test)
print('iq_no_cycle: ' + str(regression_model_mae))


# In[54]:


X = df_train_sj_inter[new_cols_with_cyclical_features]
y = df_train_sj_inter[['total_cases']]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
regression_model = LinearRegression()
regression_model.fit(X_train, y_train)

y_predict = correct_prediction(regression_model.predict(X_test))
regression_model_mae = mean_absolute_error(y_predict, y_test)
print('sj_cycle: ' + str(regression_model_mae))


# In[55]:


X = df_train_iq_inter[new_cols_with_cyclical_features]
y = df_train_iq_inter[['total_cases']]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
regression_model = LinearRegression()
regression_model.fit(X_train, y_train)

y_predict = correct_prediction(regression_model.predict(X_test))
regression_model_mae = mean_absolute_error(y_predict, y_test)
print('iq_cycle: ' + str(regression_model_mae))


# In[56]:


X = df_train_sj_inter_no_multi
y = df_train_sj_inter[['total_cases']]

min_max_scaler = preprocessing.MinMaxScaler()
np_scaled = min_max_scaler.fit_transform(X)
X = pd.DataFrame(np_scaled)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
regression_model = LinearRegression()
regression_model.fit(X_train, y_train)

y_predict = correct_prediction(regression_model.predict(X_test))
regression_model_mae = mean_absolute_error(y_predict, y_test)
print('sj_cycle_no_multi_scaled: ' + str(regression_model_mae))


# In[57]:


X = df_train_iq_inter_no_multi
y = df_train_iq_inter[['total_cases']]

min_max_scaler = preprocessing.MinMaxScaler()
np_scaled = min_max_scaler.fit_transform(X)
X = pd.DataFrame(np_scaled)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
regression_model = LinearRegression()
regression_model.fit(X_train, y_train)

y_predict = correct_prediction(regression_model.predict(X_test))
regression_model_mae = mean_absolute_error(y_predict, y_test)
print('iq_cycle_no_multi_scaled: ' + str(regression_model_mae))


# In[58]:


X = df_train_sj_inter_no_multi
y = df_train_sj_inter[['total_cases']]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
regression_model = LinearRegression()
regression_model.fit(X_train, y_train)

y_predict = correct_prediction(regression_model.predict(X_test))
regression_model_mae = mean_absolute_error(y_predict, y_test)
print('sj_cycle_no_multi: ' + str(regression_model_mae))


# In[59]:


X = df_train_iq_inter_no_multi
y = df_train_iq_inter[['total_cases']]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
regression_model = LinearRegression()
regression_model.fit(X_train, y_train)

y_predict = correct_prediction(regression_model.predict(X_test))
regression_model_mae = mean_absolute_error(y_predict, y_test)
print('iq_cycle_no_multi: ' + str(regression_model_mae))


# XGBOOST TESTES BELOW:

# In[60]:


import xgboost as xgb


# In[61]:


X = df_train_sj_inter_no_multi
y = df_train_sj_inter[['total_cases']]
#data_dmatrix = xgb.DMatrix(data=X_train,label=y_train)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=100)
xg_reg = xgb.XGBRegressor()
#xg_reg = xgb.XGBRegressor(objective ='reg:linear', colsample_bytree = 0.3, learning_rate = 0.1, max_depth = 5, alpha = 10, n_estimators = 5)


# In[ ]:


xg_reg.fit(X_train,y_train)


# In[ ]:


xg_reg.fit(X_train,y_train)
y_predict = correct_prediction(xg_reg.predict(X_test))

regression_model_mae = mean_absolute_error(y_predict, y_test)
print('xg_sj_cycle_no_multi: ' + str(regression_model_mae))


