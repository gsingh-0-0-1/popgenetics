import numpy as np
import random
import matplotlib.pyplot as plt

#parameters for the simulation

#set up our alleles
DOMINANT_ALLELE = "A"
RECESSIVE_ALLELE = "a"

#probabilities of surviving given any genotype (0 - 1)
#never make the probability entire 1! Not only is this
#a bit unrealistic, it'll cause the program to hang when 
#dealing with the carrying capacity
HOM_DOM_SURV_PROB = 0.99
HET_SURV_PROB = 0.99
HOM_REC_SURV_PROB = 0.99

#number of organisms that we want to start with
N_INITIAL_ORGANISMS = 100

#number of generations to run the simulation for
N_GENERATIONS = 20

#mean and standard deviation for the number of offspring each parent produces
#we'll assume that number of offspring is normally distributed - we could
#theoretically hold it constant, but doing so would make an even larger assumption
#about the behavior of any given population
OFFSPRING_N_MEAN = 2.5
OFFSPRING_N_STD = 0.5

#carrying capacity - really, this is a computational limit, since we don't want to
#have to run list comprehensions on lists of tens of thousands of items, but we can
#also use this limit to simulate a biological or ecological carrying capacity
CARRYING_CAPACITY = 1000


#convenient functions to check what genotype an organism has
def isHomDom(org):
	if org == DOMINANT_ALLELE + DOMINANT_ALLELE:
		return True
	return False

def isHomRec(org):
	if org == RECESSIVE_ALLELE + RECESSIVE_ALLELE:
		return True
	return False

def isHet(org):
	if org == DOMINANT_ALLELE + RECESSIVE_ALLELE or org == RECESSIVE_ALLELE + DOMINANT_ALLELE:
		return True
	return False

#pick a random allele
def getRandomAllele():
	return random.choice([DOMINANT_ALLELE, RECESSIVE_ALLELE])

#generate a random genotype
def generateGenotype():
	return getRandomAllele() + getRandomAllele()

#see if an organism survives based on its genotype
def checkSurvival(genotype):
	if isHomDom(genotype):
		probability = HOM_DOM_SURV_PROB
	if isHet(genotype):
		probability = HET_SURV_PROB
	if isHomRec(genotype):
		probability = HOM_REC_SURV_PROB
	#First, we generate a random number from 0 to 1. If it's less than the
	#probability of survival, then the organism lives. If not, the organism
	#dies.
	if np.random.uniform() < probability:
		return True
	return False

#return the number of offspring that a pair of organisms will have
def howManyOffspring():
	#make sure to round the result to the nearest integer
	return round(np.random.normal(OFFSPRING_N_MEAN, OFFSPRING_N_STD))

#create offspring given two parent genotypes
def generateOffspring(p1, p2):
	#get one allele from both parents
	a1 = random.choice(list(p1))
	a2 = random.choice(list(p2))

	return a1 + a2


#create an empty list to store generational data
GENERATIONAL_DATA = []

#generate initial population of organisms - we want them to all start as heterozygotes
CURRENT_POPULATION = [DOMINANT_ALLELE + RECESSIVE_ALLELE for organism in range(N_INITIAL_ORGANISMS)]


#main simulation loop
for generation_num in range(N_GENERATIONS):
	#apply the survival conditions to the current population
	CURRENT_POPULATION = [genotype for genotype in CURRENT_POPULATION if checkSurvival(genotype)]

	#enforce the carrying capacity
	#we can re-apply the survival conditions here, so only the fittest will survive when the population
	#grows too large for the environment
	while len(CURRENT_POPULATION) > CARRYING_CAPACITY:
		CURRENT_POPULATION = [genotype for genotype in CURRENT_POPULATION if checkSurvival(genotype)]


	#make sure to update our generational data with our current population
	GENERATIONAL_DATA.append(CURRENT_POPULATION)

	#create a list where we'll create our new population
	NEW_POPULATION = []


	#generate mating pairs for organisms

	#first we generate a list of indices corresponding
	#to each organism in the current population
	l = list(range(len(CURRENT_POPULATION)))
	#now we mix them up
	random.shuffle(l)
	#then we can split the shuffled indices into two groups
	#that will "mate" with each other
	l1 = l[::2]
	l2 = l[1::2]
	#if there's ever an odd number of organisms, we'll 
	#leave one out from the process

	#loop through the lists together - we can just take the minimum of
	#the two lengths in order to keep in line with the rule of leaving out
	#the last organism in case of an odd total number
	for pairnum in range(min(len(l1), len(l2))):
		#we want to create a number of offspring in line with the parameters
		#we defined earlier
		for offspringnum in range(howManyOffspring()):
			#generate new offspring and add it to the new population
			NEW_POPULATION.append(generateOffspring(CURRENT_POPULATION[l1[pairnum]], CURRENT_POPULATION[l2[pairnum]]))

	CURRENT_POPULATION = NEW_POPULATION
	del NEW_POPULATION

#calculate some statistics from the data we've generated / collected

#create lists to keep track of how many of each genotype are present
#over each generation
HOM_DOM_N_LIST = []
HET_N_LIST = []
HOM_REC_N_LIST = []

#create lists to keep track of allele frequency over each generation
DOM_ALLELE_N = []
REC_ALLELE_N = []

#keep track of total organisms over time
TOTAL_N = []

for generation in GENERATIONAL_DATA:
	#tracker variables for each genotype and allele
	n_hom_dom = 0
	n_het = 0
	n_hom_rec = 0
	n_dom_allele = 0
	n_rec_allele = 0

	#for each organism, we check its phenotype and use that
	#to increment our counter variables
	for organism in generation:
		if isHomDom(organism):
			n_hom_dom += 1
			n_dom_allele += 2
			continue
		if isHomRec(organism):
			n_hom_rec += 1
			n_rec_allele += 2
			continue
		if isHet(organism):
			n_het += 1
			n_dom_allele += 1
			n_rec_allele += 1

	#make sure to append percentage values to our lists, not the raw ones

	#get the totals
	genotype_total = n_hom_dom + n_het + n_hom_rec
	allele_total = n_dom_allele + n_rec_allele

	#if any of the totals are 0, the population has died out - we can just
	#set the totals to 1 to avoid a division by zero error, and the graphs will display
	#that there are no organisms have died out
	if genotype_total == 0 or allele_total == 0:
		genotype_total = 1
		allele_total = 1

	HOM_DOM_N_LIST.append(100 * n_hom_dom / genotype_total)
	HOM_REC_N_LIST.append(100 * n_hom_rec / genotype_total)
	HET_N_LIST.append(100 * n_het / genotype_total)

	DOM_ALLELE_N.append(100 * n_dom_allele / allele_total)
	REC_ALLELE_N.append(100 *n_rec_allele / allele_total)

	TOTAL_N.append(genotype_total)

#set up the title string for the plot to display the survival probabilities
#and the total number of generations
titlestr = DOMINANT_ALLELE + DOMINANT_ALLELE + " = " + str(HOM_DOM_SURV_PROB)
titlestr += ", " + DOMINANT_ALLELE + RECESSIVE_ALLELE + " = " + str(HET_SURV_PROB)
titlestr += ", " + RECESSIVE_ALLELE + RECESSIVE_ALLELE + " = " + str(HOM_REC_SURV_PROB)
titlestr += ", " + str(N_GENERATIONS) + " Generations"

plt.suptitle("Simulation Results: " + titlestr)

plt.subplot(2, 2, 1)
plt.title("Genotype Presence Over Time (%)")
plt.plot(HOM_DOM_N_LIST, color = "blue", label = "Homozygous Dominant")
plt.plot(HET_N_LIST, color = "green", label = "Heterozygous")
plt.plot(HOM_REC_N_LIST, color = "red", label = "Homozygous Recessive")
plt.legend()

plt.subplot(2, 2, 2)
plt.title("Allele Distribution Over Time (%)")
plt.plot(DOM_ALLELE_N, color = "blue", label = "Dominant Allele (p)")
plt.plot(REC_ALLELE_N, color = "red", label = "Recessive Allele (q)")
plt.legend()

plt.subplot(2, 2, 3)
plt.title("Total Population Over Time")
plt.plot(TOTAL_N, color = "blue")


plt.show()


