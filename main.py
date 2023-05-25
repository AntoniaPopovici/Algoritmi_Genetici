import random
import math


f = open("input.txt", "r")
g = open("evolution.txt", "w")


# relevant functions

def createChrom(chromLength):
    chromosome = ""
    for i in range(chromLength):
        chromosome += str(random.randint(0, 1))
    return chromosome

def chromosomeToX(chromosome, precision, xMin):
    intreg = int(chromosome, 2)
    x = xMin + intreg * pow(0.1, precision)
    return round(x, precision)

def fitness(x, a, b, c):
    f = a * (x ** 2) + b * x + c
    return f

def elitism(chromosomeArray, fitnessArray):
    maxFitness = max(fitnessArray)
    maxIndex = fitnessArray.index(maxFitness)
    maxChromosome = chromosomeArray[maxIndex]
    return maxFitness, maxChromosome


def binarySearch(intervaleProb, u, stanga, dreapta):
    if stanga > dreapta:
        return dreapta
    else:
        mijloc = int((stanga + dreapta) / 2)
        if u == intervaleProb[mijloc]:
            return mijloc
        if u < intervaleProb[mijloc]:
            return binarySearch(intervaleProb, u, stanga, mijloc - 1)
        return binarySearch(intervaleProb, u, mijloc + 1, dreapta)

def crossover(chrom1, chrom2):
    breaking_point = random.randint(0, len(chrom1))
    crossed_chrom1 = chrom1[:breaking_point] + chrom2[breaking_point:]
    crossed_chrom2 = chrom2[:breaking_point] + chrom1[breaking_point:]
    return crossed_chrom1, crossed_chrom2, breaking_point

def mutation(chromosome):
    mutatedChromosome = list(chromosome)
    geneIndex = random.randint(0, len(chromosome) - 1)  # Select a random gene index
    mutatedChromosome[geneIndex] = '1' if chromosome[geneIndex] == '0' else '0'  # Flip the selected gene
    return ''.join(mutatedChromosome)

def main():

    # variables

    sir = f.readline().strip()
    popSize = int(sir)
    sir = f.readline().strip().split(" ")
    xMin = float(sir[0])
    xMax = float(sir[1])
    sir = f.readline().strip().split(" ")
    #    a, b, c = map(int, input().split())
    coeficients = [float(x) for x in sir]
    a = coeficients[0]
    b = coeficients[1]
    c = coeficients[2]
    precision = int(f.readline().strip())
    probCrossover = float(f.readline().strip())
    probMutation = float(f.readline().strip())
    iterations = int(f.readline().strip())
    chromLength = round(math.log((xMax - xMin) * pow(10, precision), 2))

    # initialize population
    chromosomeArray = []
    fitnessArray = []
    fitnessSum = 0
    g.write("Initial population\n")

    while len(chromosomeArray) < popSize:
        chromosome = createChrom(chromLength)
        x = chromosomeToX(chromosome, precision, xMin)
        if xMin <= x <= xMax:
            g.write(str(len(chromosomeArray) + 1) + ": " + chromosome + " x = " + str(x) + " f = " + str(
                fitness(x, a, b, c)) + "\n")
            chromosomeArray.append(chromosome)
            fitnessSum += fitness(x, a, b, c)
            fitnessArray.append(fitness(x, a, b, c))

    # elitism
    maxFitness, maxChromosome = elitism(chromosomeArray[0], fitnessArray)

    # selection
    selectionProbArray = []
    g.write("\n\nSelection probabilities\n")
    g.write("-----------------------\n")

    i = 1
    for chromosome in chromosomeArray:
        x = chromosomeToX(chromosome, precision, xMin)
        formula = fitness(x, a, b, c)
        selectionProb = formula / fitnessSum
        selectionProbArray.append(selectionProb)
        g.write(f"Chromosome number {i} has the probability {selectionProb}\n")
        i += 1

    selectionProbArray.sort()

    probabilitiesIntervals = [0]
    for i in range(len(selectionProbArray)):
        probabilitiesIntervals.append(selectionProbArray[i] + probabilitiesIntervals[-1])

    g.write("\nProbabilities intervals selected:\n")
    for i in range(len(probabilitiesIntervals)):
        g.write(str(probabilitiesIntervals[i]) + "\n")

    chromModified = []
    for idx in range(popSize):
        u = random.random()
        index = binarySearch(probabilitiesIntervals, u, 0, len(probabilitiesIntervals))
        g.write("\nFor u = " + str(u) + " we select the chromosome " + str(index + 1))
        chromModified.append(chromosomeArray[index])

    g.write("\nAfter selection: \n")
    for i in range(len(chromModified)):
        chromosome = chromModified[i]
        x = chromosomeToX(chromosome, precision, xMin)
        fit = fitness(x, a, b, c)
        g.write(str(i + 1) + ": " + chromModified[i] + " x = " + str(x) + " f = " + str(fit) + "\n")

    g.write(f"\nCrossover probability: {probCrossover}\n")
    recombinedChromosomes = []

    for i in range(len(chromModified)):
        chromosome = chromModified[i]
        u = random.random()
        if u < probCrossover:
            recombinedChromosomes.append(chromosome)
            g.write(f"index: {i + 1} chromosome: {chromosome} u = {u} < {probCrossover} participates\n")
        else:
            g.write(f"index: {i + 1} chromosome: {chromosome} u = {u}\n")

    g.write("\n")

    crossover_pairs = []

    while len(recombinedChromosomes) >= 2:
        x = random.sample(recombinedChromosomes, 2)
        recombinedChromosomes.remove(x[0])
        recombinedChromosomes.remove(x[1])
        chromModified.remove(x[0])
        chromModified.remove(x[1])
        new_chrom1, new_chrom2, break_point = crossover(x[0], x[1])
        crossover_pairs.append((x[0], x[1], new_chrom1, new_chrom2))
        g.write(f"Recombination between chromosome {i + 1} and chromosome {i + 2}:\n")
        g.write(f"{x[0]} {x[1]} break_point {break_point}\n")
        g.write(f"Result:    {new_chrom1} {new_chrom2}\n")
        chromModified.append(new_chrom1)
        chromModified.append(new_chrom2)

    g.write("\nAfter crossover:\n")
    i = 1
    for chromosome in chromModified:
        x = chromosomeToX(chromosome, precision, xMin)
        fit = fitness(x, a, b, c)
        g.write(f"{i}: {chromosome} x = {x} f = {fit}\n")
        i = i + 1

    # mutation
    g.write(f"\nMutation probability for each gene: {probMutation}\n")
    g.write("Chromosomes modified:\n")
    mutatedIndices = []

    for i in range(len(chromModified)):
        u = random.random()
        if u < probMutation:
            mutatedIndices.append(i)
            g.write(str(i + 1) + "\n")

    chromModifiedAgain=[]
    g.write("After mutation:\n")
    for i in range(len(chromModified)):
        chromosome = chromModified[i]
        if i in mutatedIndices:
            chromosome = mutation(chromosome)
        x = chromosomeToX(chromosome, precision, xMin)
        chromModifiedAgain.append(chromosome)
        fit = fitness(x, a, b, c)
        g.write(f"{i + 1}: {chromosome} x = {x} f = {fit}\n")

    print(chromModified)
    print(chromModifiedAgain)
    chromModified=chromModifiedAgain
    print(chromModified)


    for i in range(len(chromModified)):
        cr = chromModified[i]
        x = chromosomeToX(cr,precision,xMin)
        ok = 0
        if fitness(x,a,b,c) > maxFitness:
            ok = 1
            break
    if ok == 0:
        chromModified[0] = maxChromosome

    def findMax(population):
        f_max = 0
        for i in range(len(population)):
            cr = population[i]
            x = chromosomeToX(cr,precision,xMin)
            fit = fitness(x,a,b,c)
            if fit > f_max:
                f_max = fit
        return f_max

    g.write("\n\nMaximum evolution:\n")
    g.write(str(findMax(chromModified)) + "\n")

    def sumFit(population):
        suma = 0
        for idx in range(len(population)):
            cr = population[idx]
            x = chromosomeToX(cr,precision,xMin)
            fit = fitness(x,a,b,c)
            suma += fit
        return suma

    def createNewPop(lista_crom):
        sumaF = sumFit(lista_crom)
        fit_max = 0
        selectionProbArray = []

        i = 1
        for chromosome in chromosomeArray:
            x = chromosomeToX(chromosome, precision, xMin)
            formula = fitness(x, a, b, c)
            selectionProb = formula / fitnessSum
            selectionProbArray.append(selectionProb)
            i += 1

        selectionProbArray.sort()
        probabilitiesIntervals = [0]
        for i in range(len(selectionProbArray)):
            probabilitiesIntervals.append(selectionProbArray[i] + probabilitiesIntervals[-1])


        chromModified = []
        for idx in range(popSize):
            u = random.random()
            index = binarySearch(probabilitiesIntervals, u, 0, len(probabilitiesIntervals))
            chromModified.append(chromosomeArray[index])

        for i in range(len(chromModified)):
            chromosome = chromModified[i]
            x = chromosomeToX(chromosome, precision, xMin)
            fit = fitness(x, a, b, c)

        recombinedChromosomes = []

        for i in range(len(chromModified)):
            chromosome = chromModified[i]
            u = random.random()
            if u < probCrossover:
                recombinedChromosomes.append(chromosome)

        crossover_pairs = []

        while len(recombinedChromosomes) >= 2:
            x = random.sample(recombinedChromosomes, 2)
            recombinedChromosomes.remove(x[0])
            recombinedChromosomes.remove(x[1])
            chromModified.remove(x[0])
            chromModified.remove(x[1])
            new_chrom1, new_chrom2, break_point = crossover(x[0], x[1])
            crossover_pairs.append((x[0], x[1], new_chrom1, new_chrom2))
            chromModified.append(new_chrom1)
            chromModified.append(new_chrom2)

        i = 1
        for chromosome in chromModified:
            x = chromosomeToX(chromosome, precision, xMin)
            fit = fitness(x, a, b, c)
            i = i + 1

        # mutation
        mutatedIndices = []

        for i in range(len(chromModified)):
            u = random.random()
            if u < probMutation:
                mutatedIndices.append(i)

        chromModifiedAgain=[]
        for i in range(len(chromModified)):
            chromosome = chromModified[i]
            if i in mutatedIndices:
                chromosome = mutation(chromosome)
            x = chromosomeToX(chromosome, precision, xMin)
            chromModifiedAgain.append(chromosome)
            fit = fitness(x, a, b, c)

        chromModified = chromModifiedAgain

        if fit_max > findMax(chromModified):
            chromModified = maxChromosome
        return chromModified

    idx = 0
    pop = chromModified
    while idx != iterations:
        pop = createNewPop(pop)
        g.write(str(findMax(pop)) + "\n")
        idx += 1
if __name__ == '__main__':
    main()