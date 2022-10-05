from tracemalloc import start
import numpy as np


class antColony:

    def __init__(self, antNumbers, maxRanking, startingPheromone, evapore, alpha, beta, input_file):
        self.evaporeRate = evapore
        self.alpha = alpha
        self.beta = beta
        self.solutions = [[] for i in range(antNumbers)]
        self.adj_mat = self.graph_input(input_file)
        self.antNumbers = antNumbers
        self.node_numbers = np.shape(self.adj_mat)[0]
        self.startingPheromone = startingPheromone
        self.pheromones = np.full_like(self.adj_mat, self.startingPheromone)
        self.bestSolution = []
        self.bestSolutionLength = -1
        self.bestSolutionWeigth = 1.25
        if ((maxRanking-1) > self.antNumbers):
            self.maxRank = self.antNumbers
        else:
            self.maxRank = maxRanking
        self.Q = 1
        self.constantACS = self.startingPheromone

    def __functionF(self, length):
        return self.Q / length

    def __clear_solutions(self):
        for solution in self.solutions:
            solution.clear()

    def __reset(self):
        self.solutions.clear()
        self.solutions = [[] for i in range(self.antNumbers)]
        self.bestSolution.clear()
        self.bestSolutionLength = -1
        self.pheromones.fill(self.startingPheromone)

    def graph_input(self, name):
        with open(name) as file:
            graph = np.loadtxt(file)
            print(graph)
            return graph

    def print_solution(self):
        for antIndex in range(self.antNumbers):
            print(self.solutions[antIndex])
            print(self.calculate_solution_length(antIndex))
        print(self.pheromones)

    def print_bestSolution(self):
        print(self.bestSolution, self.bestSolutionLength)

    def calculate_solution_length(self, antIndex):
        length = 0
        for nodeIndex in range(len(self.solutions[antIndex])-1):
            length = length + self.adj_mat[self.solutions[antIndex][nodeIndex],
                                           self.solutions[antIndex][nodeIndex+1]]
        return length

    def __evaporatePheromons(self):
        for i in range(self.node_numbers):
            for j in range(i, self.node_numbers):
                if (i != j):
                    self.pheromones[i, j] = (
                        1-self.evaporeRate) * self.pheromones[i, j]
                    self.pheromones[j, i] = self.pheromones[i, j]

    def __pheromonsAntSystem(self, antIndex):
        self.__evaporatePheromons()

        currentSolutionLength = self.calculate_solution_length(antIndex)
        for nodeIndex in range(len(self.solutions[antIndex])-1):
            startingNode = self.solutions[antIndex][nodeIndex]
            targetNode = self.solutions[antIndex][nodeIndex+1]
            self.pheromones[startingNode, targetNode] = self.pheromones[startingNode,
                                                                        targetNode] + (self.__functionF(currentSolutionLength))
            self.pheromones[targetNode,
                            startingNode] = self.pheromones[startingNode, targetNode]
        if self.bestSolutionLength == -1:
            self.bestSolution = self.solutions[antIndex].copy()
            self.bestSolutionLength = currentSolutionLength
        else:
            if self.bestSolutionLength > currentSolutionLength:
                self.bestSolution = self.solutions[antIndex].copy()
                self.bestSolutionLength = currentSolutionLength

    def __bestSolutionUpdateEAS(self):
        self.__evaporatePheromons()
        for nodeIndex in range(len(self.bestSolution)-1):
            startingNode = self.bestSolution[nodeIndex]
            targetNode = self.bestSolution[nodeIndex]
            self.pheromones[startingNode, targetNode] = self.pheromones[startingNode,
                                                                        targetNode] + self.bestSolutionWeigth * self.__functionF(self.bestSolutionLength)
            self.pheromones[targetNode,
                            startingNode] = self.pheromones[startingNode, targetNode]

    def __pheromonsElitistAntSystem(self, antIndex):
        self.__evaporatePheromons()

        currentSolutionLength = self.calculate_solution_length(antIndex)
        for nodeIndex in range(len(self.solutions[antIndex])-1):
            startingNode = self.solutions[antIndex][nodeIndex]
            targetNode = self.solutions[antIndex][nodeIndex+1]
            self.pheromones[startingNode, targetNode] = self.pheromones[startingNode,
                                                                        targetNode] + (self.__functionF(currentSolutionLength))
            self.pheromones[targetNode,
                            startingNode] = self.pheromones[startingNode, targetNode]
        if self.bestSolutionLength == -1:
            self.bestSolution = self.solutions[antIndex].copy()
            self.bestSolutionLength = currentSolutionLength
        else:
            if self.bestSolutionLength > currentSolutionLength:
                self.bestSolution = self.solutions[antIndex].copy()
                self.bestSolutionLength = currentSolutionLength
        if (antIndex == (self.antNumbers-1)):
            self.__bestSolutionUpdateEAS()

    def __calculate_probabilityPheOnly(self, available_nodes, i):
        pheromons_AvailabeNodesSum = np.sum([self.pheromones[i, node]
                                            for node in available_nodes])
        probability = [(self.pheromones[i, selected_note] / pheromons_AvailabeNodesSum)
                       for selected_note in available_nodes]
        return probability

    def __calculate_probabilityAlphaBeta(self, available_nodes, i):
        etapheromons_AvailableNodes = [
            np.float_power(self.pheromones[i, node], self.alpha)*np.float_power((1 / self.adj_mat[i, node]), self.beta) for node in available_nodes]
        etapheromons_AvailableNodesSum = np.sum(
            etapheromons_AvailableNodes
        )
        a = [etapheromon / etapheromons_AvailableNodesSum for etapheromon in etapheromons_AvailableNodes]
        a_sum = np.sum(a)
        probabilty = [a_ij / a_sum for a_ij in a]
        return probabilty

    def __calculate_probability(self, available_nodes, antIndex, i):
        return self.__calculate_probabilityAlphaBeta(available_nodes, i)

    def calculate_probabilityShortestPath(self, antIndex, i):
        available_nodes = [node for node in range(
            self.node_numbers)if node not in self.solutions[antIndex] and self.adj_mat[i, node] != 0]
        probability = self.__calculate_probability(
            available_nodes, antIndex, i)
        # print(available_nodes)
        # print(probability)
        return available_nodes, probability

    def __run_AntShortestPath(self, antIndex, startingNode, targetNode):
        for i in range(self.node_numbers):
            self.solutions[antIndex].append(startingNode)
            if (startingNode == targetNode):
                break
            available_nodes, probability = self.calculate_probabilityShortestPath(
                antIndex, startingNode)
            if (len(available_nodes) == 0):
                break
            startingNode = np.random.choice(available_nodes, p=probability)

    def calculate_probabilityTSP(self, antIndex, i):
        available_nodes = [node for node in range(
            self.node_numbers)if (node not in self.solutions[antIndex] and self.adj_mat[i, node] != 0)]
        probability = self.__calculate_probability(
            available_nodes, antIndex, i)
        return available_nodes, probability

    def __run_AntTSP_AS(self, antIndex, startingNode):
        while True:
            self.solutions[antIndex].append(startingNode)
            if (len(self.solutions[antIndex]) == self.node_numbers and self.adj_mat[self.solutions[antIndex][0], self.solutions[antIndex][len(self.solutions[antIndex])-1]] != 0):
                self.solutions[antIndex].append(self.solutions[antIndex][0])
                break
            available_nodes, probability = self.calculate_probabilityTSP(
                antIndex, startingNode)
            if (len(available_nodes) == 0):
                break
            startingNode = np.random.choice(available_nodes, p=probability)

    def run_SearchShortestPath(self, startingNode, targetNode, numberRuns):
        self.__reset()
        for i in range(numberRuns):
            for antIndex in range(self.antNumbers):
                self.__run_AntShortestPath(antIndex, startingNode, targetNode)

            for antIndex in range(self.antNumbers):
                if (self.solutions[antIndex][len(self.solutions[antIndex])-1] == targetNode):
                    self.__update_pheromons(antIndex)
            # self.print_solution()
            self.__clear_solutions()

    def run_SearchTSP_AS(self, startingNode, numberRuns):
        self.__reset()
        for i in range(numberRuns):
            for antIndex in range(self.antNumbers):
                self.__run_AntTSP_AS(antIndex, startingNode)

            for antIndex in range(self.antNumbers):
                if (self.solutions[antIndex][len(self.solutions[antIndex])-1] == self.solutions[antIndex][0]):
                    self.__pheromonsAntSystem(antIndex)
            # self.print_solution()
            self.__clear_solutions()

    def run_SearchTSP_EAS(self, startingNode, numberRuns):
        self.__reset()
        for i in range(numberRuns):
            for antIndex in range(self.antNumbers):
                self.__run_AntTSP_AS(antIndex, startingNode)

            for antIndex in range(self.antNumbers):
                if (self.solutions[antIndex][len(self.solutions[antIndex])-1] == self.solutions[antIndex][0]):
                    self.__pheromonsElitistAntSystem(antIndex)
            # self.print_solution()
            self.__clear_solutions()

    def __sortRankingSolutions(self, antIndex, rankingSolutions, rankingSolutionsLength):
        currentSolutionLength = self.calculate_solution_length(antIndex)
        copied_solution = self.solutions[antIndex].copy()
        if (len(rankingSolutionsLength) == 0):
            rankingSolutions.append(copied_solution)
            rankingSolutionsLength.append(currentSolutionLength)
        else:
            for i in range(len(rankingSolutions)):
                if rankingSolutionsLength[i] > currentSolutionLength:
                    rankingSolutions.insert(i, copied_solution)
                    rankingSolutionsLength.insert(i, currentSolutionLength)
                    break
            if self.solutions[antIndex] not in rankingSolutions:
                rankingSolutions.append(copied_solution)
                rankingSolutionsLength.append(currentSolutionLength)
        while (len(rankingSolutions) > self.maxRank):
            rankingSolutions.pop()
            rankingSolutionsLength.pop()

    def __update_pheromons_RAS(self, rankingSolutions, rankingSolutionsLengths):
        ranking = 0
        for solution in rankingSolutions:
            self.__evaporatePheromons()
            for nodeIndex in range(len(solution) - 1):
                weightRanking = self.maxRank - ranking
                startingNode = solution[nodeIndex]
                targetNode = solution[nodeIndex+1]
                self.pheromones[startingNode,
                                targetNode] = self.pheromones[startingNode, targetNode] + weightRanking * self.__functionF(rankingSolutionsLengths[ranking])
                self.pheromones[targetNode,
                                startingNode] = self.pheromones[startingNode, targetNode]
            ranking += 1

    def __print_RAS(self, rankingSolutions, rankingSolutionsLengths):
        for count, solution in enumerate(rankingSolutions):
            print(solution, rankingSolutionsLengths[count])

    def run_SearchTSP_RAS(self, startingNode, numberRuns):
        self.__reset()
        rankingSolutions = []
        rankingSolutionsLengths = []
        for i in range(numberRuns):
            for antIndex in range(self.antNumbers):
                self.__run_AntTSP_AS(antIndex, startingNode)

            for antIndex in range(self.antNumbers):
                if (self.solutions[antIndex][len(self.solutions[antIndex])-1] == self.solutions[antIndex][0]):
                    self.__sortRankingSolutions(
                        antIndex, rankingSolutions, rankingSolutionsLengths)
            self.__update_pheromons_RAS(
                rankingSolutions, rankingSolutionsLengths)
            self.__clear_solutions()

        self.bestSolution = rankingSolutions[0].copy()
        self.bestSolutionLength = rankingSolutionsLengths[0]

    def __calculate_probabilityACS(self, antIndex, i):
        available_nodes = [node for node in range(
            self.node_numbers)if (node not in self.solutions[antIndex] and self.adj_mat[i, node] != 0)]
        etapheromons_AvailableNodes = [
            np.float_power(self.pheromones[i, node], self.alpha)*np.float_power((1 / self.adj_mat[i, node]), self.beta) for node in available_nodes]
        return available_nodes, etapheromons_AvailableNodes

    def __update_pheromones_ACS(self, bestSolution, bestSolutionLength):
        functionResult = self.__functionF(bestSolutionLength)
        for node_index in range(len(bestSolution)-1):
            startingNode = bestSolution[node_index]
            targetNode = bestSolution[node_index+1]
            self.pheromones[startingNode, targetNode] = (
                1-self.evaporeRate) * self.pheromones[startingNode, targetNode] + self.evaporeRate * (functionResult)
            self.pheromones[targetNode,
                            startingNode] = self.pheromones[startingNode, targetNode]

    def __update_construction_step(self, startingNode, targetNode):
        self.pheromones[startingNode, targetNode] = (
            1-self.evaporeRate) * self.pheromones[startingNode, targetNode] + self.evaporeRate*self.constantACS
        self.pheromones[targetNode,
                        startingNode] = self.pheromones[startingNode, targetNode]

    def __run_AntTSP_ACS(self, antIndex, startingNode):
        while True:
            self.solutions[antIndex].append(startingNode)
            if (len(self.solutions[antIndex]) == self.node_numbers and self.adj_mat[self.solutions[antIndex][0], self.solutions[antIndex][len(self.solutions[antIndex])-1]] != 0):
                self.solutions[antIndex].append(self.solutions[antIndex][0])
                break
            available_nodes, probability = self.__calculate_probabilityACS(
                antIndex, startingNode)
            if (len(available_nodes) == 0):
                break
            maxValue = probability[0]
            targetNode = available_nodes[0]
            for count, node in enumerate(available_nodes):
                if (probability[count] > maxValue):
                    targetNode = node
                    maxValue = probability[count]
            self.__update_construction_step(startingNode, targetNode)
            startingNode = targetNode

    def run_SearchTSP_ACS(self, startingNode, numberRuns):
        self.__reset()
        bestSolution = []
        bestSolutionLength = -1
        for i in range(numberRuns):
            for antIndex in range(self.antNumbers):
                self.__run_AntTSP_ACS(antIndex, startingNode)
                currentLength = self.calculate_solution_length(antIndex)
                if (bestSolutionLength == -1):
                    bestSolution = self.solutions[antIndex].copy()
                    bestSolutionLength = currentLength
                elif bestSolutionLength > currentLength:
                    bestSolution = self.solutions[antIndex].copy()
                    bestSolutionLength = currentLength
            self.__update_pheromones_ACS(bestSolution, bestSolutionLength)
        self.bestSolution = bestSolution.copy()
        self.bestSolutionLength = bestSolutionLength
