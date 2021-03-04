import UMIsimulator.ampcycle as ampcycle
import UMIsimulator.dedup as dedup
from collections import Counter
import numpy as np


def generate_UMIs():
	np.random.seed(123)
	pool = ampcycle.Pool()
	pool.createUMIsWithBias(100, 0.7, 1, 0.001)
	pool.SimulateCycleWithErrors()
	pool.PCRcyclesWithErrorsAndBias(5)
	return pool

def test_dedup_naive():
	pool = generate_UMIs()
	count = dedup.dedup_naive(pool.umis_amplified)
	assert count == 111

def test_dedup_hierarchical():
	pool = generate_UMIs()
	count = dedup.dedup_hierarchical(pool.umis_amplified, method="single")
	assert count == 98
	count = dedup.dedup_hierarchical(pool.umis_amplified, method="complete")
	assert count == 101
	count = dedup.dedup_hierarchical(pool.umis_amplified, method="ward")
	assert count == 101

def test_dedup_unique():
	pool = generate_UMIs()
	count = dedup.dedup_unique(pool.umis_amplified)
	assert count == 126

def test_dedup_percentile():
	pool = generate_UMIs()
	count = dedup.dedup_percentile(pool.umis_amplified)
	assert count == 126

def test_dedup_graph():
	pool = generate_UMIs()
	count = dedup.dedup_graph(pool.umis_amplified)
	assert count == 98

def test_dedup_adj():
	pool = generate_UMIs()
	count = dedup.dedup_adj(pool.umis_amplified)
	assert count == 98

def test_dedup_dir_adj():
	pool = generate_UMIs()
	count = dedup.dedup_dir_adj(pool.umis_amplified)
	assert count == 99