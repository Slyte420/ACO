import colony

if __name__ == "__main__":
    input_file = "graph1.in"
    evapore_rate = 0.01
    alpha = 1
    beta = 1
    antNumbers = 250
    maxRanking = 25
    startingPheromone = 1
    starting_node = 0
    target_node = 4
    numberRuns = 250
    ACO = colony.antColony(antNumbers, maxRanking, startingPheromone,
                           evapore_rate, alpha, beta, input_file)
    ACO.run_SearchTSP_AS(starting_node, numberRuns)
    ACO.print_bestSolution()
    ACO.run_SearchTSP_EAS(starting_node, numberRuns)
    ACO.print_bestSolution()
    ACO.run_SearchTSP_RAS(starting_node, numberRuns)
    ACO.print_bestSolution()
    ACO.run_SearchTSP_ACS(starting_node,
                          numberRuns)
    ACO.print_bestSolution()
    #ACO.run_SearchShortestPath(starting_node, target_node, numberRuns)
    # ACO.print_bestSolution()
