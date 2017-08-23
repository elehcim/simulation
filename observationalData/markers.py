#markers = {22000: '*', 25000: '*', 24000: 's', 14000: 'p', 27000: '^', 26000: 'v', 21000: 'H', 29000: 'o', 30000: 'D', 33000: '<', 40000: '>'}

def getMarker(sim):
	basesim = sim // 1000 * 1000
	markers = {33000: 'o', 25000: '^', 24000: 'v', 29000: 's', 14000: 'D', 30000: 'p', 21000: 'H', 40000: '8'}
	simsNP3 = (27000, 34000, 32000, 31000)

	if basesim in markers:
		return markers[basesim]
	elif basesim in simsNP3
		return '*'
	else
		return 'x'
