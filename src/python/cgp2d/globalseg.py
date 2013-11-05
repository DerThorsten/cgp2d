


"""
X  are region  variables


INPUT:
	E_G(X) :
		- is a very high order function which should be minimized

	DIFFERENT_REDUCED_SEG  = { watershed, mc on gradient , etc }

WORKING DATA:
	P_L(X) :
		- is the probability desity  which will be modified 

		- it can be a graphical model of any structure

		- the structure of the graphical model
		  should resamble somehow E_G(X)



ALGORITHM:
	
Initialization : 
	- evaluate all seg in DIFFERENT_REDUCED_SEG

	- get elite samples from evaluation

	- update probability from elite samples

Iteration :
	


"""


# defined by user

def globalEnergyFunction(state):
    pass


def sampleStateFromProbability(stateMean,stateStd):

    nVar  = len(stateMean)

    # probablity of beeing on  to maximize
    edgeProbability = numpy.ones(nvar)

    # energy of beeing on
    edgeEnergy 
    



def optimizer(nVar,initStd=0.5,damping=0.5,globalEnergyFunction,stateToStructure,sampleStateFromProbability):

    nSamples      = 1000
    nEliteSamples = 20 
    maxIterations = 100

    
    stateMean =  numpy.zeros(nVar)+0.5
    stateStd  =  numpy.zeros(nVar)+initStd

    bestState = numpy.zeros(nVar,dtype=uint32)
    bestValue = float('inf')

    for iteration in range(maxIterations):

        # get samples states from probablity desnisty
        sampleStates = [sampleStateFromProbability(stateMean,stateStd) for n in range(nSamples)]

        # evaluate samples
        sampleEnergy = numpy.array([ globalEnergyFunction(s) for s in sampleStates ])

        # sort samples by energy
        sortedIndex  = numpy.argsort(sampleEnergy)

        # get elitem samples
        eliteSamples  = numpy.array(sampleStates[sortedIndex[0:nEliteSamples]])

        if(sampleEnergy[sortedIndex[0]]<bestValue):
            bestState[:] = eliteSamples[0][:]


        # shift probability density 
        stateMean = stateMean*damping + (1.0-damping)*numpy.mean(eliteSamples,axis=1)
        stateStd  = stateStd*damping  + (1.0-damping)*numpy.stddev(eliteSamples,axis=1)



