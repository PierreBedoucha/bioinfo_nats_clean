# Statistical learning: the setting and the estimator
# object in scikit-learn

# Datasets

from sklearn import datasets

iris = datasets.load_iris()
data = iris.data
print(data.shape)

# Reshaping data with the digits dataset
digits = datasets.load_digits()
print(digits.images.shape)
# (1797, 8, 8)
import matplotlib.pyplot as plt

plt.imshow(digits.images[-1], cmap=plt.cm.cmap_d)

# To use this dataset with scikit-learn, we transform
# each 8x8 image into a feature vector of length 64
data = digits.images.reshape((digits.images.shape[0], -1))

# Estimators objects

# An estimator is any object that learns from data
estimator = Estimator(param1=1, param2=2)
estimator.fit(data)
print(estimator.param1)

# Supervised learning: predicting an output variable
# from high-dimensional observations

# Nearest neighbor and the curse of dimensionality
import numpy as np
from sklearn import datasets

iris = datasets.load_iris()
iris_X = iris.data
iris_y = iris.target
np.unique(iris_y)

# k-Nearest neighbors classifier
# Split iris data in train and test data
# A random permutation, to split the data randomly
np.random.seed(0)
indices = np.random.permutation(len(iris_X))
iris_X_train = iris_X[indices[:-10]]
iris_y_train = iris_y[indices[:-10]]
iris_X_test = iris_X[indices[-10:]]
iris_y_test = iris_y[indices[-10:]]
# Create and fit a nearest-neighbor classifier
from sklearn.neighbors import KNeighborsClassifier

knn = KNeighborsClassifier()
knn.fit(iris_X_train, iris_y_train)

# Lnear model: from regression to sparsity

diabetes = datasets.load_diabetes()
diabetes_X_train = diabetes.data[:-20]
diabetes_X_test  = diabetes.data[-20:]
diabetes_y_train = diabetes.target[:-20]
diabetes_y_test  = diabetes.target[-20:]

#Linear regression
from sklearn import linear_model
 regr = linear_model.LinearRegression()
 regr.fit(diabetes_X_train, diabetes_y_train)
# LinearRegression(copy_X=True, fit_intercept=True, n_jobs=None,
#                  normalize=False)
print(regr.coef_)

# The mean square error
 np.mean((regr.predict(diabetes_X_test) - diabetes_y_test)**2)
 
# Explained variance score: 1 is perfect prediction
 # and 0 means that there is no linear relationship
 # between X and y.
 regr.score(diabetes_X_test, diabetes_y_test)

# Shrinkage
X = np.c_[ .5, 1].T
 y = [.5, 1]
 test = np.c_[ 0, 2].T
 regr = linear_model.LinearRegression()

import matplotlib.pyplot as plt
 plt.figure()

np.random.seed(0)
 for _ in range(6):
    this_X = .1 * np.random.normal(size=(2, 1)) + X
    regr.fit(this_X, y)
    plt.plot(test, regr.predict(test))
    plt.scatter(this_X, y, s=3)
    
alphas = np.logspace(-4, -1, 6)
 print([regr.set_params(alpha=alpha)
            .fit(diabetes_X_train, diabetes_y_train)
            .score(diabetes_X_test, diabetes_y_test)
        for alpha in alphas])
 
# Sparsity
regr = linear_model.Lasso()
 scores = [regr.set_params(alpha=alpha)
               .fit(diabetes_X_train, diabetes_y_train)
               .score(diabetes_X_test, diabetes_y_test)
           for alpha in alphas]
 best_alpha = alphas[scores.index(max(scores))]
 regr.alpha = best_alpha
 regr.fit(diabetes_X_train, diabetes_y_train)

# Classification
log = linear_model.LogisticRegression(solver='lbfgs', C=1e5,
                                       multi_class='multinomial')
 log.fit(iris_X_train, iris_y_train)


# Support Vector Machine (SVM)
from sklearn import svm
>>> svc = svm.SVC(kernel='linear')
>>> svc.fit(iris_X_train, iris_y_train)

# Using kernel
svc = svm.SVC(kernel='linear')
svc = svm.SVC(kernel='poly', degree=3)
svc = svm.SVC(kernel='rbf')

