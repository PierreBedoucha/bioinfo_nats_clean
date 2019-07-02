import matplotlib.pyplot as plt
import pandas as pd

# Read in the data to a pandas DataFrame using the read_csv method.
train = pd.read_csv('titanic_data.csv')

# We are using the train data as placeholder for the test data,
# so we can implement code to run predictions on the test pandas DataFrame
test = pd.read_csv('titanic_data.csv')

# Uncomment this line to read in the true test, it will be revealed in due time....
# test=pd.read_csv('titanic_test.csv')
# train

# EXAMPLE
# This loop over both the train and test dataframes
for f in [train, test]:
    # This will add a new column with the name 'child' containing all 0
    f['child'] = 0
    # This will change the values of rows that statisfy the train['age']<10 condition
    f['child'].loc[f['age'] < 10] = 1
    f['adult_male'] = 0
    # This will change the values of rows that statisfy the train['age']> 16 and train['sex']=='male' condition,
    # the paranteses around the expression are important.
    f['adult_male'].loc[(f['age'] > 16) & (f['sex'] == 'male')] = 1
    ########
    # Now try adding the other examples

    f['adult'] = 0
    f['adult'].loc[f['age'] > 16] = 1
    f['adult_female'] = 0
    f['adult_female'].loc[(f['age'] > 16) & (f['sex'] == 'female')] = 1
    f['female'] = 0
    f['male'] = 0
    f['female'].loc[(f['sex'] == 'female')] = 1
    f['male'].loc[(f['sex'] == 'male')] = 1

### And maybe some more features


# train
# print test[['embarked','embarked_S','embarked_C','embarked_Q']]

# make cross_val sets by grouping ticket number
ticket_number_cv = {}
cv = 1
train['cv'] = None
cv_count = {}
for i, row in train.iterrows():
    if cv not in cv_count:
        cv_count[cv] = 0
    if row['ticket'] not in ticket_number_cv:
        ticket_number_cv[row['ticket']] = cv
        cv = cv + 1
        if cv > 5:
            cv = 1
    train.loc[i, 'cv'] = ticket_number_cv[row['ticket']]
    cv_count[ticket_number_cv[row['ticket']]] = cv_count[ticket_number_cv[row['ticket']]] + 1

print(train.shape)
print(cv_count)

## SCALING

from sklearn import preprocessing

min_max_scaler = preprocessing.MinMaxScaler()

# Pick the columns that have numbers that can be used for training:
trainable_cols = ["age", "fare", "pclass", "has_cabin_number", "male", "female", "sibsp", "parch", "child", "adult",
                  "adult_female", "adult_male"]

# Example with three features, do try other features.
# trainable_cols=["age","fare","male"]

# For the training we include the 'cv' column
train_columns = trainable_cols + ["survived", "cv"]
test_columns = trainable_cols + ["survived"]

# Make train_data and test_data to be used for training and testing.
train_data = train[train_columns].dropna()
test_data = test[test_columns].dropna()
# df=train_target[trainable_cols]

scaling = True
if scaling:
    columns_to_scale = trainable_cols
    # Fit the scaler on the training data
    min_max_scaler.fit(train_data[columns_to_scale].values)
    # Transform the scaling to the train_data
    train_data.loc[:, columns_to_scale] = min_max_scaler.transform(train_data[columns_to_scale].values)
    # Transform the scaling to the test_data
    test_data.loc[:, columns_to_scale] = min_max_scaler.transform(test_data[columns_to_scale].values)

## Setting up the data for training

from sklearn.model_selection import PredefinedSplit

(size_x, size_y) = train_data.shape
target_index = size_y - 2
cv_index = size_y - 1

# Put the training data in X the .values method returns a numpy matrix of the numbers in the DataFrame.
X = train_data[trainable_cols].values  # ,0:target_index]
print(X)

# Put the target value in Y
Y = train_data['survived'].values

# Use the PredefinedSplit class to define the cross-validation sets from before.
cv = PredefinedSplit(train_data['cv'].values)

# Training Machine Learning methods

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score

# This is one example, change here to try other classifiers, what are the default hyperparameters?
# clf = DecisionTreeClassifier(max_depth=None)
clf = RandomForestClassifier(max_depth=None)
# clf = RandomForestClassifier(n_estimators=100, max_depth=None, min_samples_split=2, n_jobs=4,random_state=None,verbose=0)


# Some dictionaries to store cross-validated predictions


pred_save = []
true_save = []
pred_prob_save = []

for i, (train_index, val_index) in enumerate(cv.split(), 1):
    print("Set: ", i)
    print("Training on", len(train_index), "examples")
    print("Testing on", len(val_index), "examples")
    (X_train, X_val) = X[train_index, :], X[val_index, :]
    (Y_train, Y_val) = Y[train_index], Y[val_index]
    # print X_train.shape
    # print Y_train.shape
    # train_pred=clf.predict(X_train)
    # acc_train=sklearn.metrics.accuracy_score(train_pred,Y_train)
    # print acc_train
    clf = clf.fit(X_train, Y_train)
    #   continue

    # Predict on the training data
    pred = clf.predict(X_train)
    # Calculate performance measures on the validation data
    acc_train = accuracy_score(pred, Y_train)
    mcc_train = matthews_corrcoef(pred, Y_train)
    f1_train = f1_score(pred, Y_train)

    # Predict on the validation data
    val_pred = clf.predict(X_val)
    # print val_pred
    # Predict the probability (to use the roc-plot later)
    val_pred_prob = val_pred

    # val_pred_prob=clf.predict_proba(X_val)

    # Save the values to have predictions for all folds.
    pred_save.append(val_pred)
    pred_prob_save.append(val_pred_prob)
    true_save.append(Y_val)
    # Calculate performance measures on the validation data
    acc = accuracy_score(val_pred, Y_val)
    mcc = matthews_corrcoef(val_pred, Y_val)
    f1 = f1_score(val_pred, Y_val)

    print("Training performance", "f1", f1_train, "acc", acc_train, "mcc", mcc_train)
    print("Validation performance", "f1", f1, "acc", acc, "mcc", mcc)
    print("==============")

# Calculate overall validation performance
predictions = np.concatenate(pred_save)
correct = np.concatenate(true_save)
predicted_prob = np.concatenate(pred_prob_save)
acc = accuracy_score(predictions, correct)
mcc = matthews_corrcoef(predictions, correct)
f1 = f1_score(predictions, correct)
print("==============")
print("Overall Validation Performance", "f1", f1, "acc", acc, "mcc", mcc)
print("==============")

## Check the performance

import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_curve

# Define a classifier
clf = RandomForestClassifier(n_estimators=100, max_depth=None, min_samples_split=2, n_jobs=4, random_state=None,
                             verbose=0)

# Use the cv.split() method to generate an iterator of the cross-validation sets.


# Some dictionaries to store cross-validated predictions
predictions = {}
correct = {}
predicted_prob = {}
pred_sorted = []
importances = []
pred_save = []
true_save = []
pred_prob_save = []
legend_text = []

# trainable_cols=["age","fare","pclass","has_cabin_number","male","female","sibsp","parch","child","adult","adult_female","adult_male"]
# trainable_cols=["pclass","male"]

for feat_stop in range(0, X.shape[1] + 1):
    name = "-".join(trainable_cols[0:feat_stop + 1])
    # name
    legend_text.append(name)  # trainable_cols[feat_stop-1]) #len(feat_name))
    pred_save = []
    true_save = []
    pred_prob_save = []
    for i, (train_index, val_index) in enumerate(cv.split(), 1):
        # print "Set: ",i, name
        # print "Training on",len(train_index),"examples"
        # print "Testing on",len(val_index),"examples"
        (X_train, X_val) = X[train_index, 0:feat_stop + 1], X[val_index, 0:feat_stop + 1]
        (Y_train, Y_val) = Y[train_index], Y[val_index]
        # print X_train.shape
        # print Y_train.shape
        # train_pred=clf.predict(X_train)
        # acc_train=sklearn.metrics.accuracy_score(train_pred,Y_train)
        # print acc_train
        clf = clf.fit(X_train, Y_train)
        #   continue

        # Predict on the training data
        pred = clf.predict(X_train)
        # Calculate performance measures on the validation data
        acc_train = accuracy_score(pred, Y_train)
        mcc_train = matthews_corrcoef(pred, Y_train)
        f1_train = f1_score(pred, Y_train)

        # Predict on the validation data
        val_pred = clf.predict(X_val)
        # Predict the probability (to use the roc-plot later)
        val_pred_prob = clf.predict_proba(X_val)
        # Save the values to have predictions for all folds.
        pred_save.append(val_pred)
        pred_prob_save.append(val_pred_prob)
        true_save.append(Y_val)
        # Calculate performance measures on the validation data
        acc = accuracy_score(val_pred, Y_val)
        mcc = matthews_corrcoef(val_pred, Y_val)
        f1 = f1_score(val_pred, Y_val)

        # print "Training performance","f1",f1_train,"acc",acc_train,"mcc",mcc_train
        # print "Validation performance","f1",f1,"acc",acc,"mcc",mcc

    # Calculate overall validation performance
    predictions[name] = np.concatenate(pred_save)
    correct[name] = np.concatenate(true_save)
    predicted_prob[name] = np.concatenate(pred_prob_save)
    acc = accuracy_score(predictions[name], correct[name])
    mcc = matthews_corrcoef(predictions[name], correct[name])
    f1 = f1_score(predictions[name], correct[name])
    print("==============")
    print("Training on", name)
    print("Overall Validation Performance", "f1", f1, "acc", acc, "mcc", mcc)
    print("==============")

    plt.clf()

    legend_text = []
    pred_sorted = sorted(predictions.items(), key=lambda kv: (len(kv[1]), kv[0]))
    print(pred_sorted)

    print(clf.feature_importances_)

    importances = clf.feature_importances_
    std = np.std([tree.feature_importances_ for tree in clf.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]

    # Print the feature ranking
    print("Feature ranking:")

    for f in range(X_train.shape[1]):
        print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

for (name, value) in pred_sorted:
    # print key, value
    # continue
    legend_text.append(name)
    acc = accuracy_score(predictions[name], correct[name])
    mcc = matthews_corrcoef(predictions[name], correct[name])
    f1 = f1_score(predictions[name], correct[name])
    # (prec,recall,thres)=precision_recall_curve(true_save,pred_prob_save[:,1])
    (fpr, tpr, thres_roc) = roc_curve(correct[name], predicted_prob[name][:, 1])
    plt.plot(fpr, tpr)
    plt.title('ROC curve')
    plt.xlabel('fpr')
    plt.ylabel('tpr')

    plt.savefig('RF.png', dpi=300)
    plt.legend(legend_text, bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=4, borderaxespad=0.)  # mode="expand"
    plt.show()

# Plot the feature importances of the forest
plt.figure()
plt.title("Feature importances")
plt.bar(range(X_train.shape[1]), importances[indices],
        color="r", yerr=std[indices], align="center")
plt.xticks(range(X_train.shape[1]), indices)
plt.xlim([-1, X_train.shape[1]])

plt.show()

# from sklearn.metrics import matthews_corrcoef
# trainable_cols=["age","fare","pclass","has_cabin_number","male","female","sibsp","parch","passenger_on_ticket","kids_on_ticket","teens_on_ticket","adults_on_ticket","child"]
# trainable_cols.reverse()
trainable_cols.append("survived")
trainable_cols.append("cv")
# train[trainable_cols].dropna().to_csv('titanic_data.csv')
# print trainable_cols
# trainable_cols=["age","fare","pclass","survived","cv"]
# tmp=train[["cv","survived"]]
# corr=np.corrcoef(tmp)
# df_clean=train[trainable_cols].dropna()
# df_clean["survived"].unique
# df_clean.to_excel('tmp.xls')
# g = sns.factorplot("cv", "survived", "male", data=df_clean, kind="bar", palette="muted", legend=True)
# plt.show()

print(trainable_cols)

data = train[trainable_cols].dropna().values.astype(float)
print(data.shape[1])
(size_x, size_y) = data.shape
X = data[:, 0:size_y - 2]
Y = data[:, size_y - 2]
CV_splits = data[:, size_y - 1]
print(X)
# print data.shape
corr = np.corrcoef(data)
# corr=np.corrcoef(data, rowvar=False)
# print corr.shape

# print sklearn.metrics
# np.set_printoptions(precision=3)
# print(corr)
for i, name in enumerate(trainable_cols):
    #    c=matthews_corrcoef(data[:,i],data[:,size_y-1])
    print(i, name, corr[i, size_y - 1])

plt.clf()
legend_text = []
measures = {}

measures['MC_train'] = {}
measures['MC_test'] = {}
measures['acc_train'] = {}
measures['acc_test'] = {}
measures['F1_train'] = {}
measures['F1_test'] = {}
for feat_stop in range(X.shape[1], X.shape[1] + 1):
    feat_name = "-".join(trainable_cols[0:feat_stop + 1])
    legend_text.append(trainable_cols[feat_stop - 1])  # len(feat_name))
    print(feat_name)

    for trees in (1, 10, 100, 1000):
        clf = RandomForestClassifier(n_estimators=trees, max_depth=None, min_samples_split=2, n_jobs=4,
                                     random_state=None, verbose=0)
        fold = 1
        pred_save = []
        true_save = []
        pred_prob_save = []
        for train_index, test_index in cv.split():
            print(trees, fold, len(train_index), len(test_index))
            # continue
            print(feat_name, feat_stop)

            X_train, X_test = X[train_index, 0:feat_stop + 1], X[test_index, 0:feat_stop + 1]
            Y_train, Y_test = Y[train_index], Y[test_index]
            # continue
            # print feat_name,feat_stop,X_test.shape
            # continue
            clf = clf.fit(X_train, Y_train)
            #   continue
            train_pred = clf.predict(X_train)
            test_pred = clf.predict(X_test)
            test_pred_prob = clf.predict_proba(X_test)
            acc_train = sklearn.metrics.accuracy_score(train_pred, Y_train)
            MC_train = sklearn.metrics.matthews_corrcoef(train_pred, Y_train)
            acc = sklearn.metrics.accuracy_score(train_pred, Y_train)
            MC = sklearn.metrics.matthews_corrcoef(train_pred, Y_train)

            name = feat_name + '-trees-' + str(trees) + '-fold-' + str(fold)

            measures['MC_train'][name] = sklearn.metrics.matthews_corrcoef(train_pred, Y_train)
            measures['MC_test'][name] = sklearn.metrics.matthews_corrcoef(test_pred, Y_test)
            measures['acc_train'][name] = sklearn.metrics.accuracy_score(train_pred, Y_train)
            measures['acc_test'][name] = sklearn.metrics.accuracy_score(test_pred, Y_test)
            measures['F1_train'][name] = sklearn.metrics.f1_score(train_pred, Y_train)
            measures['F1_test'][name] = sklearn.metrics.f1_score(test_pred, Y_test)
            pred_save.append(test_pred)
            pred_prob_save.append(test_pred_prob)
            true_save.append(Y_test)
            fold = fold + 1
        # continue
        name = feat_name + '-trees-' + str(trees) + '-overall'
        # print Y_test
        # print np.concatenate(true_save)
        pred_save = np.concatenate(pred_save)
        true_save = np.concatenate(true_save)
        pred_prob_save = np.concatenate(pred_prob_save)
        measures['MC_test'][name] = sklearn.metrics.matthews_corrcoef(pred_save, true_save)
        measures['acc_test'][name] = sklearn.metrics.accuracy_score(pred_save, true_save)
        measures['F1_test'][name] = sklearn.metrics.f1_score(pred_save, true_save)
        (prec, recall, thres) = sklearn.metrics.precision_recall_curve(true_save, pred_prob_save[:, 1])
        (fpr, tpr, thres_roc) = sklearn.metrics.roc_curve(true_save, pred_prob_save[:, 1])
        plt.plot(fpr, tpr)
        plt.title('ROC curve')
        plt.xlabel('fpr')
        plt.ylabel('tpr')
        # plt.xlabel('Prec')
        # plt.ylabel('Recall')
print(legend_text)
plt.legend(legend_text)
plt.show()
