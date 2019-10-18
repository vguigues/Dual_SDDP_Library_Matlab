#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: lingquan
"""
import pandas
import numpy
from msppy.msp import MSLP
from msppy.solver import SDDP, Extensive
import gurobipy

T = 4
size = 25
hydro_ = pandas.read_csv("./data/hydro.csv", index_col=0)
demand = pandas.read_csv("./data/demand.csv", index_col=0)
deficit_ = pandas.read_csv("./data/deficit.csv", index_col=0)
exchange_ub = pandas.read_csv("./data/exchange.csv", index_col=0)
thermal_ = [
    pandas.read_csv("./data/thermal_{}.csv".format(i), index_col=0)
    for i in range(4)
]
scenarios = numpy.empty((T-1,size,4))
for t in range(1,T):
    scenarios[t-1] = pandas.read_csv("/Users/lingquan/Desktop/thesis/vincent/stage_{}.csv".format(t+1), index_col=0)        
    
HydroThermal = MSLP(T=T, discount=0.9906, bound=0)
for t in range(T):
    m = HydroThermal[t]
    stored_now,stored_past = m.addStateVars(4, ub=hydro_['UB'][:4], name="stored")
    spill = m.addVars(4, name="spill")
    hydro = m.addVars(4, ub=hydro_['UB'][-4:], name="hydro")
    deficit = m.addVars(
        [(i,j) for i in range(4) for j in range(4)],
        ub = [
            demand.iloc[t][i] * deficit_['DEPTH'][j]
            for i in range(4) for j in range(4)
        ],
        obj = [
            deficit_['OBJ'][j]
            for i in range(4) for j in range(4)
        ],
        name = "deficit")
    thermal = [None] * 4
    for i in range(4):
        thermal[i] = m.addVars(
            len(thermal_[i]),
            ub=thermal_[i]['UB'],
            lb=thermal_[i]['LB'],
            obj=thermal_[i]['OBJ'],
            name="thermal_{}".format(i)
        )
    exchange = m.addVars(5,5, ub=exchange_ub.values.flatten(), name="exchange")
    thermal_sum = m.addVars(4, name="thermal_sum")
    m.addConstrs(thermal_sum[i] == gurobipy.quicksum(thermal[i].values()) for i in range(4))
    for i in range(4):
        m.addConstr(
            thermal_sum[i]
            + gurobipy.quicksum(deficit[(i,j)] for j in range(4))
            + hydro[i]
            - gurobipy.quicksum(exchange[(i,j)] for j in range(5))
            + gurobipy.quicksum(exchange[(j,i)] for j in range(5))
            == demand.iloc[t][i]
        )
    m.addConstr(
        gurobipy.quicksum(exchange[(j,4)] for j in range(5))
        - gurobipy.quicksum(exchange[(4,j)] for j in range(5))
        == 0
    )
    for i in range(4):
        if t == 0:
            m.addConstr(
                stored_now[i] + spill[i] + hydro[i] - stored_past[i]
                == hydro_['INITIAL'][4:8][i]
            )
        else:
            m.addConstr(
                stored_now[i] + spill[i] + hydro[i] - stored_past[i] == 0,
                uncertainty = {'rhs': scenarios[t-1,:,i]}
            )
    if t == 0:
        m.addConstrs(stored_past[i] == hydro_['INITIAL'][:4][i] for i in range(4))
HydroThermal_ext = Extensive(HydroThermal)
HydroThermal_ext.solve()
HydroThermal_SDDP = SDDP(HydroThermal)
HydroThermal_SDDP.solve(
    max_time=200,
    n_processes=3,
    n_steps=3,
)
