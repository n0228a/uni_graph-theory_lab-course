
from sklearn import svm
from sklearn.model_selection import train_test_split
import numpy as np
from scripts import create_varied_set
import glob
import pandas as pd
from math import log, ceil

### INTRODUCTION
# There are 50 different classes with 1000 reactions each in the dataset.

# The repetitions depends on the reactions per class chosen. As 1000 is the maximum it is cumbersome to do many repetitions as the dataset won't vary much.
### END INTRODUCTION

# START CONFIGURATION VARIABLES
CLASSES = 5
REACTIONS_PER_CLASS = 200
FEATURE_SET = 'DRF Shortest Paths'
# END CONFIGURATION VARIABLES

dataset = pd.DataFrame()

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

def run_single_experiment(feature_set:str, chosen_classes:list, reactions_per_class:int):
    # 2nd step: Create a varied dataset according to configuration variables
    data = create_varied_set(dataset, chosen_classes=chosen_classes, reactions_per_class=reactions_per_class)

    # 3rd step: Prepare data for SVM
    X = data[f'{feature_set}']
    Y = data['rxn_class']

    # 4th step: Preprocess strings from the excel file into sets if possible (not so for e.g. NaN) # TODO: Add a counter for empty sets or some other sort of tracking
    for idx in X.index:
        x = X[idx]
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
    return clf.score(X_test, Y_test)


def run_experiments(feature_set:str = FEATURE_SET, used_classes:int = CLASSES, reactions_per_class:int = REACTIONS_PER_CLASS):
    scores = []
    chosen_classes = dataset["rxn_class"].drop_duplicates().sample(n=used_classes, random_state=42)
    repetitions = required_repetitions(reactions_per_class)

    # TODO: Introduce for loop over different class sets to evaluate dependence on chosen classes

    for i in range(repetitions):
        score = run_single_experiment(feature_set, chosen_classes, reactions_per_class)
        scores.append(score)
    
    # Summary
    avg_score = sum(scores) / len(scores)
    # We need: Feature Set, used classes, number of used classes, reactions per class, repetitions, average score, single scores
    print(f"Feature Set: {feature_set}, Used Classes: {chosen_classes.tolist()}, Number of Used Classes: {used_classes}, Reactions per Class: {reactions_per_class}, Repetitions: {repetitions}, Average Score: {avg_score}, Scores: {scores}")
    
def required_repetitions(sample_size:int, target_coverage:float = 0.95) -> int:
    """
    Calculate the number of repetitions required to achieve a target coverage
    of unique reactions in the dataset.

    Args:
        sample_size: Number of reactions sampled in each experiment.
        target_coverage: Desired coverage of unique reactions (between 0 and 1).

    Returns:
        Number of repetitions needed to achieve the target coverage.
    """

    # Using the formula: repetitions = log(1 - target_coverage) / log(1 - (sample_size / total_reactions))
    p = sample_size / 1000 # Total reactions per class is 1000
    if p >= 1.0:
        return 1
    
    repetitions = log(1 - target_coverage) / log(1 - p)

    return ceil(repetitions)


# 1st step: Read files and load to dataset DataFrame
with open("data/combined_data.xlsx", "rb") as f:
    dataset = pd.read_excel(f)

run_experiments()



# TODO: Add used classes to compare if it depends on the classes chosen - Carefully! Needs to be evaluated on all combinations
# TODO: May use other learning techniques