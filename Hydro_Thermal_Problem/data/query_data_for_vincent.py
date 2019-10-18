#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: lingquan
"""
import pandas
import numpy

T = 12
size = 50
hydro_ = pandas.read_csv("./data/hydro.csv", index_col=0)
demand = pandas.read_csv("./data/demand.csv", index_col=0)
deficit_ = pandas.read_csv("./data/deficit.csv", index_col=0)
exchange_ub = pandas.read_csv("./data/exchange.csv", index_col=0)
thermal_ = [
    pandas.read_csv("./data/thermal_{}.csv".format(i), index_col=0)
    for i in range(4)
]
hist = [
    pandas.read_csv("./data/hist_{}.csv".format(i), sep=";")
    for i in range(4)
]
hist = pandas.concat(hist, axis=1)
hist.dropna(inplace=True)
hist.drop(columns='YEAR', inplace=True)
scenarios = numpy.array([
    hist.iloc[:,12*i:12*(i+1)].transpose().values for i in range(4)
])
log_scenarios = numpy.log(scenarios)
scenarios_resampled = numpy.empty((T-1,size,4))
for t in range(1,T):
    mean = numpy.mean(log_scenarios[:,t%12,:], axis=1)
    cov = numpy.cov(log_scenarios[:,t%12,:])
    scenarios_resampled[t-1] = numpy.exp(
        numpy.random.RandomState(0).multivariate_normal(
            mean=mean,
            cov=cov,
            size=size
        )
    )
for i in range(4):
    pandas.DataFrame(numpy.round(scenarios_resampled[:,:,i],1)).to_csv(
        "/Users/lingquan/Desktop/thesis/vincent/reservior_{}.csv".format(i)
    )
