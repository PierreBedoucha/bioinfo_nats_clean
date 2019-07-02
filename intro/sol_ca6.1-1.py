#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 10:52:37 2017

@author: Yaxin Xue, yaxin.xue@uib.no
"""

# import all the packages
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split

# four machine learning models
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier

# evaluation methods
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score

# task1: read the file
def load_dataset(in_file):
    colnames = ['sepal-length', 'sepal-width', 'petal-length', 'petal-width', 'class']
    iris_data = pd.read_csv(in_file, names = colnames)
    return (iris_data)


# task2: summarize the dataset
def summarize_dataset(dataset):
    # print dimension
     n_samples, n_features = dataset.shape
     print("Dataset has", n_samples, "rows and", n_features,"columns.")
     # print top 5 lines
     top_5 = dataset.head(5)
     print(top_5)
     # print mean and standard variation of sepal-width 
     all_mean = dataset.mean()
     sw_mean = all_mean['sepal-width']
     all_sd = dataset.std()
     sw_sd = all_sd['sepal-width']
     print("Mean of sepal-width is", sw_mean)
     print("Standard variation of sapal-width is", sw_sd)
     # draw box plot of each attribute
     dataset.plot(kind='box', subplots=True, layout=(2,2), sharex=False, sharey=False, 
                  title = 'Boxplot of each feature')
     plt.show()
    
# task3: run one model
def run_model(dataset):
    # step1: split the columns of dataset into features and target classes
    array = dataset.values
    X = array[:, 0:4]
    Y = array[:, 4]
    # step2: split the dataset into train and validation , and set parameters
    validation_size = 0.50
    seed = 9
    scoring = 'accuracy'
    X_train, X_val, Y_train, Y_val = train_test_split(X, Y, test_size = validation_size, random_state = seed)
    # step 3: fit the model, here we use LogisticRegression()    
    lr_fit = LogisticRegression().fit(X_train,Y_train)
    # step4: test model with validation dataset
    lr_pred = lr_fit.predict(X_val)
    # step5: evaluate the performance
    print("The classification result of LR is:")
    print(classification_report(Y_val, lr_pred))
        
# task4: compare models and let user input fraction for validation
def cmp_models(dataset):
    # step1: split the dataset into train and validation: ask user to set up validation size
    array = dataset.values
    X = array[:, 0:4]
    Y = array[:, 4]
    validation_size = eval(input("input the percentage of validation dataset(0.0 - 1.0):"))
    seed = 10
    scoring = 'accuracy'
    X_train, X_val, Y_train, Y_val = train_test_split(X, Y, test_size = validation_size, random_state = seed)
    # step2: select models, here we compare four models
    models = []
    models.append(('LR', LogisticRegression()))
    models.append(('NB', GaussianNB()))
    models.append(('SVM', SVC()))
    models.append(('KNN', KNeighborsClassifier()))
    accuracy = dict()
   # step3: build a for loop to run and compare each model
    for name, model in models:
        # fit new_model with training dataset
        new_model = model.fit(X_train, Y_train)
        # test new_model with validation dataset
        new_pred = new_model.predict(X_val)
        # evaluate the performance 
        acc = accuracy_score(Y_val, new_pred)
        accuracy[name] = acc
        print("The acccuracy of", name, "is:", acc)
   # step4: answer which one is the winner based on the accuracy score
    print("The best machine learning method for my dataset is: ", max(accuracy, key=accuracy.get))

# run the program    
iris_data = load_dataset('iris_data.csv')

summarize_dataset(iris_data)
run_models(iris_data)
cmp_models(iris_data)

### TASK 5 ###
# explain which one is best for high fraction, which for low fraction.
