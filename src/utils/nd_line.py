# -*- coding: utf-8 -*-
"""Created on Tue Sep  7 14:07:13 2021.

@author: dpm42
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splev, splprep
from scipy.spatial import distance


class nd_line:
    def __init__(self, points, inplace=False):
        self.points = np.array([tuple(x) for x in points])
        alldist = self._lengths(self.points)
        self.length = sum(alldist)
        self.cumul = np.cumsum([0] + alldist)
        self.type = 'linear'

    def _lengths(self, points):
        """calculate the length (sum of the euclidean distance between
        points)"""
        return [self.e_dist(points[i], points[i + 1]) for i in range(len(points) - 1)]

    def _length(self, points):
        """calculate the length (sum of the euclidean distance between
        points)"""
        return sum([self.e_dist(points[i], points[i + 1]) for i in range(len(points) - 1)])

    def interp(self, dist):
        """return a point a specified distance along the line."""
        assert dist <= self.length, 'length cannot be greater than line length'
        assert dist >= 0, 'length cannot be less than zero'
        if dist == 0:
            return self.points[0]
        if dist == self.length:
            return self.points[-1]
        index = np.where(self.cumul < dist)[0][-1]
        d = self.cumul[index]
        vector = (self.points[index + 1] - self.points[index]) / self.e_dist(self.points[index], self.points[index + 1])
        remdist = dist - d
        final_point = remdist * vector + self.points[index]
        return final_point

    def interp_rat(self, ratio):
        assert ratio <= 1, "Ratio for interp_rat() must be a value from 0 to 1"
        return self.interp(ratio * self.length)

    def splineify(self, samples=None, s=0):
        """Turn line into a spline approximation, currently occurs in place."""
        if samples is None:
            samples = len(self.points)
        tck, u = splprep([self.points[:, i] for i in range(self.points.shape[1])], s=s)
        self.points = np.transpose(splev(np.linspace(0, 1, num=samples), tck))
        self.length = self._length(self.points)
        self.type = 'spline'

    def plot2d(self):
        fig = plt.figure()
        plt.scatter(self.points[:, 0], self.points[:, 1])

    def plot3d(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(self.points[:, 0], self.points[:, 1], self.points[:, 2])

    def e_dist(self, a, b):
        return np.sqrt(np.sum((a - b) ** 2, axis=0))
