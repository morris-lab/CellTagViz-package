import numpy as np
import pandas as pd
import itertools
import collections

def findDuplicates(N, L, MOI):
    '''
    This function takes as an input the number of cells (N), the library size (L)
    and the average MOI for the virus and returns the number of duplicate cells.
    '''
    n_tags_per_cell = np.random.poisson(MOI, N)
    cell_tags = np.empty(N, dtype=object)
    rr = range(1, L + 1)
    for i in range(N):
        cell_tag_digits = np.random.randint(1, L + 1, size = n_tags_per_cell[i])
        cell_tags[i] = list(cell_tag_digits)
    
    num_duplicates = 0
    seen = set()
    for x in cell_tags:
        if len(x) > 1:
            srtd = tuple(sorted(set(x)))
            if srtd not in seen:
                seen.add(srtd)
            else:
                num_duplicates += 1

    return num_duplicates

def runDuplicatesExperiment(n_experiments, N, L, MOI):
    '''
    Runs the experiment n_experiments times to initialize the cell population.
    '''
    n_experiments = int(n_experiments)
    N = int(N)
    L = int(L)
    MOI = int(MOI)
    duplicates_experiment = np.empty(n_experiments)
    for i in range(int(n_experiments)):
        duplicates_experiment[i] = findDuplicates(N, L, MOI)
    
    #duplicates_experiment_normalized = duplicates_experiment/N * 100
    duplicates_experiment_normalized = duplicates_experiment/N 
    
    return duplicates_experiment_normalized
    
def initialCellPopulation(N, L, MOI):
    '''
    Create the starting cell population with the CellTags.
    '''
    n_tags_per_cell = np.random.poisson(MOI, N)
    cell_tags = np.empty(N, dtype=object)

    for i in range(N):
        cell_tag_digits = np.random.randint(1, L + 1, size = n_tags_per_cell[i])
        cell_tags[i] = list(cell_tag_digits)
    
    seen = set()
    clones = {}
    for x in cell_tags:
        srtd = tuple(sorted(set(x)))
        if srtd not in seen:
            seen.add(srtd)
            clones[srtd] = 1
        else:
            clones[srtd] += 1
    
    return clones    
    
def runClonalSimulation(N, L, MOI, division_rate, passage_rate, passage_fraction, sequence_time):
    '''
    Returns clonesDistribution: 
        1. First row is the number of cell tags in each clone
        2. Each subsequent row is the number of cells in each clone
        3. The number of columns is the number of different clones
    '''
    clones = initialCellPopulation(N, L, MOI)
    clonesDistribution = np.zeros((int(sequence_time/division_rate) + 2, len(clones)), dtype=int)
    
    x = 0
    for seq in clones:
        clonesDistribution[0][x] = len(seq)
        clonesDistribution[1][x] = clones[seq]
        x += 1
        
    dividing_times = list(range(division_rate, (int(sequence_time/division_rate) + 1) * division_rate, division_rate))
    passage_times = list(range(passage_rate, (int(sequence_time/passage_rate) + 1) * passage_rate, passage_rate))
    
    all_times = dividing_times + passage_times
    all_times = list(set(all_times))
    all_times.sort()
    
    y = 2
    for i in range(len(all_times)):
        if all_times[i] in dividing_times and all_times[i] in passage_times:
            # Double the cells first
            clonesDistribution[y][:] = clonesDistribution[y - 1][:] * 2
            
            # Passage the Cells
            clonesDistribution[y][:] = np.random.binomial(clonesDistribution[y][:], passage_fraction)
            y += 1
        
        else:
            if all_times[i] in dividing_times:
                clonesDistribution[y][:] = clonesDistribution[y - 1][:] * 2
                y += 1
            if all_times[i] in passage_times:
                clonesDistribution[y - 1][:] = np.random.binomial(clonesDistribution[y - 1][:], passage_fraction)
    
    return clonesDistribution
    
def meanClonePopulation(n_simulations, N, L, MOI, division_rate, passage_rate, passage_fraction, sequence_time):
    mean_clones = np.zeros(n_simulations)
    for i in range(n_simulations):
        cDist = runClonalSimulation(N, L, MOI, division_rate, passage_rate, passage_fraction, sequence_time)
        final_clone_distribution = cDist[-1][:]
        mean_clones[i] = np.mean(final_clone_distribution)
    
    return mean_clones    
    
def meanSequencedCloneSize(n_simulations, N, L, MOI, division_rate, passage_rate, passage_fraction, sequence_time,
                           n_cells_sequenced, threshold_clone_size):

    n_simulations = int(n_simulations)
    N = int(N)
    L = int(L)
    MOI = int(MOI)
    division_rate = int(division_rate)
    passage_rate = int(passage_rate)
    sequence_time = int(sequence_time)
    n_cells_sequenced = int(n_cells_sequenced)
    threshold_clone_size = int(threshold_clone_size)
    
    mean_clones = np.zeros(n_simulations)
    mean_clones_above_threshold = np.zeros(n_simulations)
    
    for i in range(n_simulations):
        cDist = runClonalSimulation(N, L, MOI, division_rate, passage_rate, passage_fraction, sequence_time)
        df = pd.DataFrame()
        df['NumTags'] = cDist[0]
        df['Final Clones'] = cDist[-1]
        df['SampledCells'] = np.zeros(len(cDist[0]))
        
        weights = cDist[-1] / np.sum(cDist[-1], dtype = np.float64)
        sampled_cells = np.random.choice(np.arange(0, np.size(cDist, axis=1)), size=n_cells_sequenced, p=weights)
        
        counter=collections.Counter(sampled_cells)
        
        temp = np.zeros(len(cDist[0]))
        temp[list(counter.keys())] = list(counter.values())
        df['SampledCells'] = temp
        
        # Remove Cells (rows) with < 2 cell tags
        df2 = df.drop(df[df.NumTags < 2].index)
        
        mean_clones[i] = np.mean(df2['SampledCells'])
        
        # Only consider clones above the specified threshold
        df3 = df2.drop(df2[df2.SampledCells < threshold_clone_size].index) 
        mean_clones_above_threshold[i] = np.mean(df3['SampledCells'])
    
    return mean_clones, mean_clones_above_threshold    
    #return mean_clones_above_threshold    
    
def cloneSizeSingleExperiment(N, L, MOI, division_rate, passage_rate, passage_fraction, sequence_time,
                              n_cells_sequenced, threshold_clone_size):

    N = int(N)
    L = int(L)
    MOI = int(MOI)
    division_rate = int(division_rate)
    passage_rate = int(passage_rate)
    sequence_time = int(sequence_time)
    n_cells_sequenced = int(n_cells_sequenced)
    threshold_clone_size = int(threshold_clone_size)
    

    cDist = runClonalSimulation(N, L, MOI, division_rate, passage_rate, passage_fraction, sequence_time)
    df = pd.DataFrame()
    df['NumTags'] = cDist[0]
    df['Final Clones'] = cDist[-1]
    df['SampledCells'] = np.zeros(len(cDist[0]))

    weights = cDist[-1] / np.sum(cDist[-1], dtype = np.float64)
    sampled_cells = np.random.choice(np.arange(0, np.size(cDist, axis=1)), size=n_cells_sequenced, p=weights)

    counter=collections.Counter(sampled_cells)

    temp = np.zeros(len(cDist[0]))
    temp[list(counter.keys())] = list(counter.values())
    df['SampledCells'] = temp

    # Remove Cells (rows) with < 2 cell tags
    df2 = df.drop(df[df.NumTags < 2].index)

    mean_clone_size = np.mean(df2['SampledCells'])

    # Only consider clones above the specified threshold
    df3 = df2.drop(df2[df2.SampledCells < threshold_clone_size].index) 
    mean_clone_size_above_threshold = np.mean(df3['SampledCells'])
    
    all_clones = np.array(df2['SampledCells'])
    
    #print(mean_clone_size)
    
    return mean_clone_size, mean_clone_size_above_threshold, all_clones
    #return all_clones
    
    
    
    
    
    
    
