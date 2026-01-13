
from sklearn import svm
from sklearn.model_selection import train_test_split
import numpy as np
from scripts import create_varied_set
import glob
import pandas as pd

def my_kernel(X1, X2):
    """
    Computes the intersection kernel between two arrays of sets.
    K(x, y) = |x intersection y|
    Args:
        X: array-like of shape (n_samples_X,) containing python sets.
        Y: array-like of shape (n_samples_Y,) containing python sets.

    Returns:
        K: Kernel matrix of shape (n_samples_X, n_samples_Y).
    """
    # Convert to list if pandas Series
    if hasattr(X1, 'tolist'):
        X1 = X1.tolist()
    if hasattr(X2, 'tolist'):
        X2 = X2.tolist()
    
    n_x1 = len(X1)
    n_x2 = len(X2)

    K = np.zeros((n_x1, n_x2))

    for i in range(n_x1):
        for j in range(n_x2):
            # Compute intersection size
            K[i, j] = len(X1[i].intersection(X2[j]))

    return K

# START SCRIPT

dataset = pd.DataFrame()

# START CONFIGURATION VARIABLES
CLASSES = 2
REACTIONS_PER_CLASS = 50
FEATURE_SET = 'DRF Nodes'
# END CONFIGURATION VARIABLES

# 1st step: Read files and load to dataset DataFrame
with open("data/combined_data.xlsx", "rb") as f:
    data = pd.read_excel(f)

# 2nd step: Create a varied dataset according to configuration variables
data = create_varied_set(data, classes=CLASSES, reactions_per_class=REACTIONS_PER_CLASS)

# 3rd step: Prepare data for SVM
X = data[f'{FEATURE_SET}']
Y = data['rxn_class']

# 4th step: Preprocess strings from the excel file into sets if possible (not so for e.g. NaN)
for idx in X.index:
    print("------------------------")
    x = X[idx]
    print(x)
    try:
        X[idx] = set(map(str, x.strip('{}').split(', ')))
    except AttributeError:
        X[idx] = set()

# 5th step: Split data into training and test set
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

#6th step: Set up SVM classifier with custom kernel
clf = svm.SVC(kernel=my_kernel) # Here, one may try different kernels, see documentation

# 7th step: Train model on training set
clf.fit(X_train, Y_train)

# 8th step: Evaluate model on test set
score = clf.score(X_test, Y_test)

# 9th step: Store setting and metric

# TODO:!

# TODO: May use other learning techniques