#!/usr/bin/env python

import random
import numpy
import sys
import copy
import argparse

class Organism:
  '''A class to contain an organism'''
  def __init__(self, cellID, cheatID = 0, genome=numpy.array([]), parent=None, empty=False, lineage = -1):
    self.age = 0
    self.cheatID = cheatID
    self.empty = empty
    self.ID = cellID
    self.fitness = 0
    self.genome = genome
    if not self.empty:
      # Making assert instead of printing error
      assert len(genome) or parent, "Must have a genome or parent set"
      if len(genome):
        self.genome = genome
      else:
        self.genome = parent.genome[:]
        self.cheatID = parent.cheatID
        self.mutate()
        parent.mutate()
        parent.fitness = 0
        parent.age = 0


        
  def __repr__(self):
    # Using string formatting
    return "empty: {}, ID: {}, genome: {}, fitness: {}\n".format(self.empty, self.ID, self.genome, self.fitness)


  def __eq__(self, other):
    # Boolean Zen
    return numpy.array_equal(self.genome, other.genome) and (self.lineage == other.lineage)


  def __hash__(self):
    full = numpy.array_str(self.genome) + str(self.lineage)
    return hash(full)


  def update(self):
    '''Updates the organism's fitness based on its age'''
    if not self.empty:
      self.age += 1
      cur_gene = self.genome[self.age % len(self.genome)]
      self.fitness += cur_gene
      return True

    return False
      

  def mutate(self):
    newGenome = numpy.copy(self.genome)
    for i in range(len(newGenome)):
      if random.random() < args.genome_mut:
        if newGenome[i] == 0:
          newGenome[i] = 1
        else:
          newGenome[i] = 0
    self.genome = newGenome
    if random.random() < args.cheat_mut:
      self.cheatID = 1

    
      
  def findNeighbors(self):
    cellID = self.ID
    radius = 1
    world_x = deme_x
    world_y = deme_y
    cell_x = cellID % world_x
    cell_y = (cellID - cell_x)/world_x
    neighbor_ids = []
    for i in range(cell_x-radius, cell_x+radius+1):
      for j in range(cell_y-radius, cell_y+radius+1):
        if i == cell_x and j == cell_y:
          continue
        if i<0:
          x = world_x + i
        elif i>=world_x:
          x = i-world_x
        else:
          x = i

        if j<0:
          y = world_y + j
        elif j>=world_y:
          y = j-world_y
        else:
          y = j

        neighbor_ids.append(y*world_x+x)

    return neighbor_ids

class Deme:
  '''A class to contain a deme.'''
  def __init__(self, deme_size, prop_cheaters=0, parent=None):
    self.orgs = []
    self.deme_size = deme_size
    if not parent:
      self.prop_cheaters = prop_cheaters

      for i in range(deme_size):
        if random.random() < self.prop_cheaters:
          newOrg = self.makeOrg(1)
          self.orgs.append(newOrg)
        else:
          newOrg = self.makeOrg(0)
          self.orgs.append(newOrg)

    else:
      copy_parents = parent.orgs[:]
      random.shuffle(copy_parents)
      parent_org = copy_parents[0]
      count = 0
      while parent_org.empty and count < len(copy_parents):
        parent_org = copy_parents[count]
        count += 1
      assert not parent_org.empty, "Empty deme!"

      newSeedOrg = Organism(0, parent = parent_org)

      self.orgs.append(newSeedOrg)
      for i in range(1,deme_size):
        self.orgs.append(Organism(len(self.orgs), empty=1))



  def makeOrg(self, cheatID):
    '''A function to make a new organism randomly'''
    randomBitArray = numpy.random.randint(2, size=(100,))
    newOrg = Organism(len(self.orgs), cheatID = cheatID, genome=randomBitArray)


    return newOrg

  def reproduceOrg(self, org):
    '''A helper function to reproduce a host organism.'''
    ## Time to reproduce
    ## Returns list of neighbor indices
    dead_neighbor = False
    neighbors = org.findNeighbors()
    for ID in neighbors:
      if self.orgs[ID].empty:
        dead_neighbor = ID
        break
    if dead_neighbor:
      position = dead_neighbor
    else:
      position = random.choice(neighbors)
    newOrg = Organism(position, parent = org)


    self.orgs[position] = newOrg

  def update(self):
    '''A function that runs a single update'''
    for org in self.orgs:
      if not org.empty:
        result = org.update()
        #check neighbor cooperating
        neighbors = org.findNeighbors()
        competitorID = random.choice(neighbors)
        competitor = self.orgs[competitorID]
        if org.cheatID > competitor.cheatID:
          org.fitness *= 2
        elif org.cheatID < competitor.cheatID:
          org.fitness *= 0.5
        elif org.cheatID == 1 and competitor.cheatID ==1:
          org.fitness *= 0.2
        elif org.cheatID == 0 and competitor.cheatID == 0:
          org.fitness = org.fitness

    for org in self.orgs:
      #Need to decide reproduction after all fitnesses have been updated to avoid weirdness
      if not org.empty and org.fitness >= len(org.genome):
        self.reproduceOrg(org)
          
  def getFitness(self):
    '''Calculates the average fitness of the deme for the competitions.'''
    sum_fit = 0.0
    for org in self.orgs:
      sum_fit += org.fitness

    if len(self.orgs):
      return sum_fit/len(self.orgs)
    else:
      return 0



        
class Population:
  '''A class to contain the population and do stuff'''
  def __init__(self, deme_size, num_demes, cheat_demes, prop_cheaters):
    self.currentUpdate = 0
    self.sub_pops = []
    self.num_demes = num_demes
    self.cheat_demes = cheat_demes
    self.prop_cheaters = prop_cheaters
    self.deme_size = deme_size

    deme_nums = range(num_demes)
    cheater_demes_list = random.sample(deme_nums, cheat_demes)

    for i in deme_nums:
      if i in cheater_demes_list:
        self.sub_pops.append(Deme(deme_size, self.prop_cheaters))
      else:
        self.sub_pops.append(Deme(deme_size, 0))

  def update(self):
    '''Go through each subpopulation/deme to have them update their orgs and migrate or reproduce demes'''
    self.currentUpdate += 1
    for sub_pop in self.sub_pops:
      sub_pop.update()

    if self.currentUpdate % args.competition_updates == 0:
      #time to compete the demes, there can be only one!
      competitors = random.sample(self.sub_pops, args.comp_deme_num)
      max_deme, min_deme = self.competeDemes(competitors)


      newDeme = Deme(self.deme_size, parent=max_deme)
      self.sub_pops[self.sub_pops.index(min_deme)]=newDeme
        
    

  def competeDemes(self, competitors):
    '''A helper function to compete whichever demes are given.'''
    max_deme_fit = -1
    max_deme = None
    min_deme_fit = 101
    min_deme = None
    for deme in competitors:
      deme_fitness = deme.getFitness()

      if deme_fitness > max_deme_fit:
        max_deme = deme
        max_deme_fit = deme_fitness
      if deme_fitness < min_deme_fit:
        min_deme = deme
        min_deme_fit = deme_fitness

      
    assert (max_deme and min_deme), "Deme competition yielded error"
    return max_deme, min_deme


  def getStats(self):
    '''Calculate the proportion of cheaters and the average fitness of cheaters and cooperators.'''
    stats = {}
    for s in range(len(self.sub_pops)):
      num_cheaters = 0.0
      num_cooperators = 0.0
      sum_cheater_fitness = 0.0
      sum_coop_fitness = 0.0
      
      for org in self.sub_pops[s].orgs:
        if org.cheatID:
          num_cheaters += 1
          sum_cheater_fitness += org.fitness
        else:
          num_cooperators += 1
          sum_coop_fitness += org.fitness
      

      if num_cheaters > 0:
        avg_cheater_fit = sum_cheater_fitness/num_cheaters
      else:
        avg_cheater_fit = 0

      if num_cooperators > 0:
        avg_cooperator_fit = sum_coop_fitness/num_cooperators
      else:
        avg_cooperator_fit = 0

      stats[s] = {"coop_fit":avg_cooperator_fit, "num_cooperators":num_cooperators, "cheater_fit": avg_cheater_fit, "num_cheaters":num_cheaters}

    return stats

def parseArgs(parser):
  ''' A function to parse the command line arguments'''
  parser.add_argument('--num_demes', help="the number of demes that should be created, default is 1",default=1, type=int)
  parser.add_argument('--cheat_demes', help="the number of demes to start off with cheaters",default=0, type=int)
  parser.add_argument("--prop_cheaters", help="the proportion of cheaters to start cheating demes with", default = 0.5, type=float)
  parser.add_argument('deme_x', help="the x value for the deme size", type=int)
  parser.add_argument('deme_y', help="the y value for the deme size", type=int)
  parser.add_argument('--cheat_mut', help="the mutation rate for turning into a cheater, default 0.25",default=0.25, type=float)
  parser.add_argument('filename', help='the name you want data written to')
  parser.add_argument('--genome_mut', help="the mutation rate for point mutations in the genome, default = 0.007", default=0.007, type=float)
  parser.add_argument("--num_updates", help="the number of updates to run for", default=1000, type=float)
  parser.add_argument("random_seed", help="the seed to set the random number generators to", type=int)
  parser.add_argument("--competition_updates", help="the number of updates between deme competitions", default = 100, type=int)
  parser.add_argument("--comp_deme_num", help="the number of demes competing in each competition, default = 2", default = 2, type=int)

  return parser.parse_args()

parser= argparse.ArgumentParser(description="Simulation evolution of prisoner's dilemma in viruses")
args = parseArgs(parser)

seed = args.random_seed
random.seed(seed)
numpy.random.seed(seed)

num_updates = args.num_updates
deme_x = args.deme_x
deme_y = args.deme_y
deme_size = deme_x*deme_y
num_demes = args.num_demes
cheat_demes = args.cheat_demes
prop_cheaters = args.prop_cheaters

data_file = open(args.filename+str(seed)+".dat", 'w')
data_file.write("Update num_cheaters avg_cheat_fit avg_coop_fit\n")


population_orgs = Population(deme_size, num_demes, cheat_demes, prop_cheaters)

for u in range(num_updates):
  population_orgs.update()
  
  if u % 10 == 0:
    stats = population_orgs.getStats()
    num_cheaters_tot = 0.0
    avg_cheater_fit_tot = 0.0
    avg_cooperator_fit_tot = 0.0
    for key in stats:
      num_cheaters_tot += stats[key]["num_cheaters"]
      avg_cheater_fit_tot += stats[key]["cheater_fit"]
      avg_cooperator_fit_tot += stats[key]["coop_fit"]

    num_cheaters = num_cheaters_tot/num_demes
    avg_cheater_fit = avg_cheater_fit_tot/num_demes
    avg_cooperator_fit = avg_cooperator_fit_tot/num_demes
    data_file.write("{} {} {} {}\n".format(u, num_cheaters, avg_cheater_fit, avg_cooperator_fit))
    
data_file.close()






