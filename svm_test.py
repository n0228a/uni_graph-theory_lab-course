
from sklearn import svm
from sklearn.model_selection import train_test_split
import pandas as pd
import glob
import random

all_x_lists = []

for file in glob.glob("data/*.xlsx"):
    data = pd.read_excel(file)
    print(f"File: {file}, Rows: {len(data)}")
    all_x_lists.append(data['DRF Nodes'].values)

X = np.concatenate(all_x_lists)
df=pd.read_csv('schneider50k.tsv', sep='\t')
Y = df['rxn_class'].values

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

def my_kernel(X1, X2):
    """
    Computes the intersection kernel between two arrays of sets.
    K(x, y) = |x intersection y|
    Args:
        X1: array-like of shape (n_samples_X,) containing python sets.
        X2: array-like of shape (n_samples_Y,) containing python sets.

    Returns:
        K: Kernel matrix of shape (n_samples_X, n_samples_Y).
    """

    n_x = len(X1)
    n_y = len(X2)

    K = np.zeros((n_x, n_y))

    for i in range(n_x):
        for j in range(n_y):
            # Compute intersection size
            K[i, j] = len(X1[i].intersection(X2[j]))

    return K


# we create an instance of SVM and fit out data.
clf = svm.SVC(kernel=my_kernel)

clf.fit(X_train, Y_train)
score = clf.score(X_test, Y_test)
