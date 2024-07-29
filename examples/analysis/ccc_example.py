"""Created on Wed Apr  3 11:25:43 2024.

@author: dpm42
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams['figure.dpi'] = 400


def concordance_correlation_coefficient(y_true, y_pred):
    """Concordance correlation coefficient."""
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    # Raw data
    dct = {"y_true": y_true, "y_pred": y_pred}
    df = pd.DataFrame(dct)
    # Remove NaNs
    df = df.dropna()
    # Pearson product-moment correlation coefficients
    y_true = df["y_true"]
    y_pred = df["y_pred"]
    cor = np.corrcoef(y_true, y_pred)[0][1]
    # Means
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)
    # Population variances
    var_true = np.var(y_true)
    var_pred = np.var(y_pred)
    # Population standard deviations
    sd_true = np.std(y_true)
    sd_pred = np.std(y_pred)
    # Calculate CCC
    numerator = 2 * cor * sd_true * sd_pred
    denominator = var_true + var_pred + (mean_true - mean_pred) ** 2

    return numerator / denominator


import matplotlib.pyplot as plt
import numpy as np


def generate_gaussian_points(n_points=100, mean=0, spread_along=1, spread_perpendicular=0.1):
    """Generate points following a 2D Gaussian distribution along y=x line,
    with independent control over the spread along and perpendicular to the line.

    Parameters:
    - n_points: Number of points to generate.
    - mean: Mean of the Gaussian distribution.
    - spread_along: Spread of the points along the y=x line.
    - spread_perpendicular: Spread of the points perpendicular to the y=x line.

    Returns:
    - A tuple of (x_points, y_points)
    """
    # Adjust the angle for rotation to align the major axis with y=x
    theta = np.pi / 4  # 45 degrees for y=x line
    rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    # Original covariance before rotation, with independent spread along and perpendicular
    original_covariance = np.array([[spread_along, 0], [0, spread_perpendicular]])

    # Rotate covariance matrix to align with y=x
    rotated_covariance = rotation_matrix @ original_covariance @ rotation_matrix.T

    # Generate the points
    points = np.random.multivariate_normal([mean, mean], rotated_covariance, n_points)

    x_points, y_points = points[:, 0], points[:, 1]
    return x_points, y_points


import seaborn as sns

sns.set(context='poster', style='white')
# Example usage
n_points = 100
mean = 5
spread_along = 3  # Spread along y=x
spread_perpendicular = 1  # Spread perpendicular to y=x
for spread_perpendicular in [0.1, 2]:
    x_points, y_points = generate_gaussian_points(n_points, mean, spread_along, spread_perpendicular)
    plt.figure(figsize=(4, 4))
    plt.scatter(x_points, y_points, s=25, color='gray')
    plt.plot([0, 10], [0, 10], 'k--')
    plt.ylim(0, 10)
    plt.xlim(0, 10)
    ccc = concordance_correlation_coefficient(x_points, y_points)
    plt.title(f'CCC={ccc:.3f}')
    plt.xticks([0, 10])
    plt.yticks([0, 10])
    plt.xlabel('Measure A')
    plt.ylabel('Measure B')
    plt.savefig(f'ccc {spread_perpendicular}.svg')
    plt.show()
