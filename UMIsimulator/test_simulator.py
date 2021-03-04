import UMIsimulator.simulator as sim
import numpy as np

def test_mp_iter():
	res = sim.mp_iter(100, 7, 5, 0, 0, 0.7, 1, seed=123)
	assert sum(res.umis_amplified.values()) == 2290
	assert len(res.umis_amplified.keys()) == 100

def test_CountsPerMethod():
	res = sim.CountsPerMethod(10, 100, 7, 5, 0, 0, 0.7, 1, ["naive", "unique"], seed=123)
	assert sum(res["naive"]) == 990
	assert sum(res["unique"]) == 1000

def test_simulationCycle():
	res = sim.simulationCycle(10, 100, 7, 5, 0, 0, 0.7, 1, ["naive", "unique"], seed=123)
	(rows, cols) = res.shape
	assert rows == 2
	assert cols == 4