import matplotlib.pyplot as plt
import numpy as np
import readmap as rm
import timeit
import random



 


def main():

    #Search
    genome = ''.join(random.choice("acgt") for i in range(500))
    readList = []
    for r in range(20):
        readList.append(("".join(random.choice("acgt") for i in range(2*r))))
    
    runtimes = runtimeSearch(genome, readList)
    x = np.array([i[0] for i in runtimes])
    y = np.array([i[1] for i in runtimes])
    plt.title("Approximate string matching algorithm")
    plt.xlabel("pattern length")
    plt.ylabel("Runtime")
    plt.plot(x,y)
    plt.show()

    genomes = []
    for j in range(20):
        genomes.append(["", ''.join(random.choice("acgt") for i in range(j*10))])
    read = "".join(random.choice("acgt") for i in range(5))
    
    runtimes = runtimeSearchGen(genomes, read)
    x = np.array([i[0] for i in runtimes])
    y = np.array([i[1] for i in runtimes])
    plt.title("Approximate string matching algorithm")
    plt.xlabel("string length")
    plt.ylabel("Runtime")
    plt.plot(x,y)
    plt.show()





def runtimeSearch(string: str, patterns: list[str]):
    """
    measures runtimes for approximate fm-index search
    """
    runtimes = []


    bwtList = rm.preprocess_genomes([["",string]])

    


    for i in range(len(patterns)):
        start = timeit.default_timer()
        matches = rm.searchPattern(patterns[i],bwtList[0],4)
        for l in matches:
            pass
        stop = timeit.default_timer()
        
        runtimes.append([len(patterns[i]), stop - start])


    return runtimes

def runtimeSearchGen(strings: list[str], pattern: str):
    """
    measures runtimes for approximate fm-index search 
    """
    runtimes = []


    bwtList = rm.preprocess_genomes(strings)

    


    for i in range(len(strings)):
        start = timeit.default_timer()
        matches = rm.searchPattern(pattern,bwtList[i],4)
        for l in matches:
            pass
        stop = timeit.default_timer()
        
        runtimes.append([len(strings[i][1]), stop - start])


    return runtimes




if __name__ == "__main__":
    main()