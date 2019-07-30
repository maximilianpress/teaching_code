#!/usr/bin/R
# A SIMPLE EVOLUTION SIMULATOR ("EvoSim")
#============================================================
# written for matlab by max press (university of washington, 
# dept. of genome sciences) 1/2011, reimplemented in R 7/2012
# for SEP demonstration with some modifications to allow more
# sophisticated fitness models and to allow some randomness in 
# selection
#============================================================
# If for some reason distribution comes up, this script is 
# free software under the GNU General Public License (GPL). That 
# is, you can do whatever you want with it,
# so long as you don't sell or copyright it.  
#============================================================


require(gplots)

# PARAMETERS, to be altered at will for different models of evolution

select_coef = 1		# strength of selection
generations = 100	# number of generations of evolution to run
mu_rate = 0.01		# rate of mutation
dev_noise = 0.01		# random noise in phenotypes from development
pop_size = 100		# size of population; MUST BE EVEN!!

fitness_funct = 'simple'	# use this parameter to change fitness regimes
				# options are 'stable' (stabilizing selection),
				# 'simple' (directional selection to be all '1's), and 'epi' (epistasis).

# initialize population matrix and data storage vectors
# generate a population of N individuals, each a binary vector(10), or 1 row 
# of an Nx10 matrix
pop = round(matrix(runif(10*pop_size),pop_size))

fitness = rep(0,pop_size)
fit = rep(0,pop_size)

Wavg = rep(0,generations)
Wstd = rep(0,generations)
Wmin = rep(0,generations)
Wmax = rep(0,generations)
gene_add = rep(0,generations)

# function FITNESS (simple, directional selection)
calcFitness = function(pop) {
    #step through population row by row 
    for (n in 1:pop_size) {
        individual = pop[n,]
        # all 1s is the ideal genotype 
        fit[n] = mean(individual) + dev_noise*rnorm(1)
	  }
	return(fit)
	}

# function FITNESS (stable)
calcStableFitness = function(pop) {
    #step through population row by row 
    for (n in 1:pop_size) {
        individual = pop[n,]
        # calculate the distance from the ideal genotype (which sums to 5)
        fit[n] = 1-abs((mean(individual)-.5)) + dev_noise*rnorm(1)
	  }
	return(fit)
	}

# function to simulate EPISTASIS
calcEpiFitness = function(pop) {
	# for each ind
	for (n in 1:pop_size) {
		individual = pop[n,]
        # make gene 1 (first entry in genotype vector) epistatic to genes 3-6 and 7-10, 
	# but in opposite directions, with only gene 2 being additive.
        fit[n] = (1-abs(individual[2]+sum(individual[1]*individual[3:6]) - sum(abs(individual[1]-1)*individual[7:10])) / 10) + dev_noise*rnorm(1)
	  }
	return(fit)
	}



# function SELECTION
select = function(pop, fitness) {

    # apply selection to whole population- sample from population,
    # while weighting the sample probability by the fitness of each 
    # individual, subject to the selection coefficient.
    selectedpop = pop[sample(1:pop_size,pop_size/2,prob=abs(fitness/sum(fitness))^select_coef),]
    return(selectedpop)    
	}

# function MUTATION
mutate = function(popul) {

# work through each element of the matrix, give each a random chance of mutating
# INEFFICIENT- but demonstrates the wackiness of this process
	for (k in 1:pop_size) {
   		for (j in 1:10) {
      		if (runif(1) < mu_rate) {
          		popul[k,j] = abs(popul[k,j]-1)
          		}
    		}
		}
		return(popul)
   }

# for each generation of evolution, calculate fitness of pop,
# apply selection and mutation:
for (i in 1:generations) {
	if (i %% 10 == 0) {
        cat(i,' GENERATIONS\n')
        }
    
    # calculate fitness for each individual - in simple case, (1-s)^x where x is the #0s in
    # vector individual
    
    if (fitness_funct == 'simple') {
    	fitness = calcFitness(pop)
    	} else if (fitness_funct == 'epi') {
    	fitness = calcEpiFitness(pop)
    	} else if (fitness_funct == 'stable') {
    	fitness = calcStableFitness(pop)
    	} 
    
    Wavg[i] = mean(fitness)
    Wstd[i] = sd(fitness)
    Wmax[i] = max(fitness)
    Wmin[i] = min(fitness)
    gene_add[i] = mean(apply(pop,2,mean))
    
    # apply selection to the population with a selection function
    selectedpop = select(pop, fitness)
   
    # selected individuals (1/2 pop) duplicate themselves 
    pop = rbind(selectedpop,selectedpop)
    
    # each offspring has a chance of mutation at each position, determined by mu_rate
	pop = mutate(pop)
    }


# plot fitness over time - print to a .pdf, and also plot to screen
pdf(gsub(':','_',paste(date(),'EvoSim.pdf')))
plotCI(1:generations,Wavg, Wstd,type='l',xlab='Generations',ylim=c(0,1),ylab='Fitness +/- SD (black, min = red, max = green), Proportion of \"1\"s (blue)')
lines(1:generations,Wmin,col = 'red')
lines(1:generations,Wmax,col='green')
lines(1:generations,gene_add,col='blue')
dev.off()
plotCI(1:generations,Wavg, Wstd,type='l',xlab='Generations',ylim=c(0,1),ylab='Fitness +/- SD (black, min = red, max = green), Proportion of \"1\"s (blue)')
lines(1:generations,Wmin,col = 'red')
lines(1:generations,Wmax,col='green')
lines(1:generations,gene_add,col='blue')



