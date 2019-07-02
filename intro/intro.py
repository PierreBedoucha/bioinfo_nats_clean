from sklearn import datasets
from sklearn import svm
from sklearn.svm import SVC

# Loading an example dataset

iris = datasets.load_iris()
digits = datasets.load_digits()

# A dataset is a dictionary like object, holding the data
# and some metadata.

# .data member is a n-samples, n-features array.
print(digits.data)
# .target member stores the response variable.
print(digits.target)
# Each original sample is an image of shape (8,8)
print(digits.images[0])

# Learning and predicting

# Here, we are given samples of each of the 10 possible classes
# (the digits zero through nine) on which we fit an estimator
# to be able to predict the classes of new samples

# An estimator is a Python object that implements the methods
# fit(x, y) and predict(T)
# Ex: Support Vector Classification
# Here gamma is set manually. But to find good values for
# these parameters, we can use tools like grid search
# or cross validation
clf = svm.SVC(gamma=0.001, C=100.)

# The fitting, or 'learning' from the model
# Here, the training set: all the images but the last one for
# the target (prediction)
clf.fit(digits.data[:-1], digits.target[:-1])

# Predicting from the last image
clf.predict(digits.data[-1:])

# Model persistence

# Save a model in scikit-learn with using 'pickle' Python's persistence
clf = svm.SVC(gamma='scale')
X, y = iris.data, iris.target
clf.fit(X, y)

import pickle

s = pickle.dumps(clf)
clf2 = pickle.loads(s)
clf2.predict(X[0:1])

# Using joblib replacement for pickle
from joblib import dump, load

dump(clf, 'filename.joblib')
clf = load('filename.joblib')

# Conventions

# Unless otherwise specified, input will be cast to float64
# Classification targets are maintained
clf = SVC(gamma='scale')
clf.fit(iris.data, iris.target)
list(clf.predict(iris.data[:3]))
# [0, 0, 0]
clf.fit(iris.data, iris.target_names[iris.target])
list(clf.predict(iris.data[:3]))
# ['setosa', 'setosa', 'setosa']

# Refitting and updating parameters

from sklearn.datasets import load_iris

X, y = load_iris(return_X_y=True)

# Hyper-parameters of an estimator can be updated
# after it has been constructed via the set_params() method.
# Calling fit() more than once will overwrite what was learned by any previous fit()

clf = SVC()
clf.set_params(kernel='linear').fit(X, y)
clf.predict(X[:5])
clf.set_params(kernel='rbf', gamma='scale').fit(X, y)
clf.predict(X[:5])

# Multiclass vs. multilabel fitting

# using multiclass classifiers, the learning and prediction task that is performed
# is dependent on the format of the target data fit upon:
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import LabelBinarizer

X = [[1, 2], [2, 4], [4, 5], [3, 2], [3, 1]]
y = [0, 0, 1, 1, 2]

classif = OneVsRestClassifier(estimator=SVC(gamma='scale',
                                            random_state=0))

classif.fit(X, y).predict(X)

#  It is also possible to fit upon a 2d array of binary label indicators
y = LabelBinarizer().fit_transform(y)
classif.fit(X, y).predict(X)

# With multilabel outputs, it is similarly possible for an instance
# to be assigned multiple labels
from sklearn.preprocessing import MultiLabelBinarizer
y = [[0, 1], [0, 2], [1, 3], [0, 2, 3], [2, 4]]
y = MultiLabelBinarizer().fit_transform(y)
classif.fit(X, y).predict(X)

# In this case, the classifier is fit upon instances each assigned
# multiple labels. The MultiLabelBinarizer is used to binarize
# the 2d array of multilabels to fit upon.
# As a result, predict() returns a 2d array with multiple
# predicted labels for each instance.