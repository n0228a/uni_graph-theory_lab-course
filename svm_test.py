
from sklearn import svm
from sklearn.model_selection import train_test_split
import pandas as pd

df=pd.read_csv('schneider50k_clean.tsv', sep='\t')
df2=pd.read_csv()

X = # Feature vectors of reactions
Y = df['rxn_class'].values
print(Y)

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
