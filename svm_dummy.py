
from sklearn import svm
from sklearn.model_selection import train_test_split

X = # Feature vectors of reactions
Y = # Vector of Target Labels

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

def my_kernel(X, Y):
    """
    Computes the intersection kernel between two arrays of sets.
    K(x, y) = |x intersection y|
    Args:
        X: array-like of shape (n_samples_X,) containing python sets.
        Y: array-like of shape (n_samples_Y,) containing python sets.

    Returns:
        K: Kernel matrix of shape (n_samples_X, n_samples_Y).
    """

    n_x = len(X)
    n_y = len(Y)

    K = np.zeros((n_x, n_y))

    for i in range(n_x):
        for j in range(n_y):
            # Compute intersection size
            K[i, j] = len(X[i].intersection(Y[j]))

    return K


# we create an instance of SVM and fit out data.
clf = svm.SVC(kernel=my_kernel)

clf.fit(X_train, Y_train)
score = clf.score(X_test, Y_test)
