import UMIsimulator.ampcycle as ampcycle

def test_createUMIsWithBias():
	pool = ampcycle.Pool()
	pool.createUMIsWithBias(10)
	assert len(pool.umis.keys()) == 10

def test_SimulateCycleWithErrors():
	pool = ampcycle.Pool()
	pool.createUMIsWithBias(10, 1, 1, 0)
	pool.SimulateCycleWithErrors()
	assert len(pool.umis_amplified.keys()) == 10
	assert sum(pool.umis_amplified.values()) == 20

def test_PCRcyclesWithErrorsAndBias():
	pool = ampcycle.Pool()
	pool.createUMIsWithBias(10, 1, 1, 0)
	pool.SimulateCycleWithErrors()
	pool.PCRcyclesWithErrorsAndBias(2)
	assert sum(pool.umis_amplified.values()) == 80

def test_downsampleUMIs():
	pool = ampcycle.Pool()
	pool.createUMIsWithBias(10, 1, 1, 0)
	pool.SimulateCycleWithErrors()
	pool.downsampleUMIs(0.1)
	assert len(pool.umis_amplified.keys()) <= 10

def test_addSequencingErrors():
	pool = ampcycle.Pool()
	pool.createUMIsWithBias(10, 1, 1, 0)
	pool.SimulateCycleWithErrors()
	pool.addSequencingErrors()
	assert len(pool.umis_amplified.keys()) >= 10