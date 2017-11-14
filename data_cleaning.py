import numpy as np
import random
import collections


# presuming import format as provided by e.g. genome2simplesttrain.py
# prep to import data (this time with some transparency)
def split_with_axe(wood, width_x=784):
    x = wood[:, :width_x]
    y = wood[:, width_x:]
    return x, y


def csv2dataset(filein, width_x=784):
    Dataset = collections.namedtuple('Dataset', ['data', 'target'])
    wood = np.genfromtxt(filein, delimiter=',')
    x, y = split_with_axe(wood, width_x)
    out = Dataset(data=x, target=y)
    return out


def score_row(rowy):
    if np.mean(rowy) >= 0.5:
        out = 1
    else:
        out = 0
    return out


def balance_mats2(x, y):
    score_cols = np.apply_along_axis(lambda w: (score_row(w), not score_row(w)), 1, y)
    score_cols = score_cols.astype(bool)
    # split both x and y by trend in y scoring (e.g. less/more than 50% cds)
    x_one = x[score_cols[:, 0], :]
    y_one = y[score_cols[:, 0], :]
    x_zero = x[score_cols[:, 1], :]
    y_zero = y[score_cols[:, 1], :]
    # shuffle everything
    x_one, y_one = shuff_mat_together(x_one, y_one)
    x_zero, y_zero = shuff_mat_together(x_zero, y_zero)
    # and shorten everything to minimum length
    target_nrows = min(x_one.shape[0], x_zero.shape[0])
    x_one = x_one[:target_nrows, :]
    y_one = y_one[:target_nrows, :]
    x_zero = x_zero[:target_nrows, :]
    y_zero = y_zero[:target_nrows, :]
    # re concatenate
    out_x, out_y = shuff_mat_together(np.vstack((x_one, x_zero)), np.vstack((y_one, y_zero)))
    return out_x, out_y


def shuff_mat_together(x, y):
    indexes = list(range(x.shape[0]))
    random.shuffle(indexes)
    x = np.array([x[i, :] for i in indexes])
    y = np.array([y[i, :] for i in indexes])
    return x, y
