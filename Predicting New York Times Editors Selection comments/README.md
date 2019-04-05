# Predicting New York Times "Editor's Selection" comments using NLP and Gradient Boosting

The current status: **Finished**

### Overview

The following notebook has the goal to predict the "Editor's Selection" comments of the open articles from the [New York Times](https://www.nytimes.com/) webpage.

The imported dataset contains data from January 2017 until May 2017 and from January 2018 until May 2018, related to articles and comments from NY Times. After importing, some features are pre-selected. The columns are treated and some new features are created. The **editorsSelection** represents if the comment is an "Editor's Selection" or not, and is the *target*.

The feature **commentBody** represents the comment itself. I separated it and processed the text using **NLP** technique **Bag-of-words**.

For the classification, the **Gradient Boosting** Machine Learning algorithm is compared with **SVM** and **Naive Bayes** and tunned with **Grid Search**. Three measures are analyzed: **Accuracy, AUROC, and F1-Score**.

Finally, the results are discussed, and some future works are recommended.

>Some codes and ideas were inspired by [Aashita Kesarwani Predicting NYT's pick notebook](https://www.kaggle.com/aashita/predicting-nyt-s-pick/notebook)

Click [here](https://github.com/bfroz/Projects/blob/master/Predicting%20New%20York%20Times%20Editors%20Selection%20comments/Predicting%20New%20York%20Times%20Editors%20Selection%20comments%20using%20NLP%20and%20Gradient%20Boost.ipynb) to check the full project.

### Results

The XGBoost results:

|Raw|-			|
|-|-----------------------|
| Mean Accuracy			| 0.91 |
| Mean AUROC  | 0.86 |
| Mean F1-Score			| 0.81 |

|Grid Search |	-			|
|-|-----------------------|
| Mean Accuracy| 0.91 |
| Mean AUROC  | 0.87 |
| Mean F1-Score			| 0.82 |

The results now are way slightly better after Grid Search the parameters. Improving the f1-score in Grid Search helped the XGBoost to penalize misclassification of the class = 1.

When dealing with such unbalanced dataset, the F1-Score is the appropriate measurement to see how the model handles overfitting. The value of 0.82 of F1-Score is good and *acceptable* to use in a production level since the **New York Times** need to look at most comments to make a good choice. If the model is used, it's possible to **reduce** the number of analyzed comments and optimize the editor's work. Also, it can be used to automatize the editor's choice.