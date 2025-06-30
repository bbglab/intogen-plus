from math import sqrt, floor
from collections import Counter


def center_of_geometry(Coordinates):
    '''
    Give me x, y, z coordinates, get back the cooresponding center of geometry
    '''

    return [sum(Coordinate)/len(Coordinate) for Coordinate in zip(*Coordinates)]


def distance(Coordinates):
    '''
    Euclidean distance between to x, y, z coordinates
    '''
    x, y, z = zip(*Coordinates)

    return sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)


def density(Distances, MinDist=0, MaxDist=10):
    '''
    For a given residues center of geometry, count the number of other residue
    center of geometries between a min and max disatance
    '''
    Counts = []

    for Distance in Distances:
       Counts.append(floor(Distance))

    Density = []
    for Distance, Count in sorted(Counter(Counts).items(), key=lambda Count: Count[0]):
        if Distance in range(MinDist, MaxDist):
            Density.append(Count)

    return sum(Density)


def cumulative_density(Distances, num_res, MinDist=0, MaxDist=10):
    '''
    For a given residues center of geometry, count the number of other residue
    center of geometries between a min and max disatance
    '''
    Counts = []

    for Distance in Distances:
       Counts.append(floor(Distance))

    Density = []
    dists = []
    widths = []
    for i, (Distance, Count) in enumerate(sorted(Counter(Counts).items(),
                                                 key=lambda Count: Count[0])):
        if Distance in range(MinDist, MaxDist):
            #dists.append(Distance)
            if Density == []:
                Density.append(Count)
                dists.append(Distance)
            else:
                Density.append(Density[-1] + Count)
                widths.append(Distance - dists[-1])
                dists.append(Distance)
    widths.append(MaxDist+1-Distance)
    cum_dens = [Density[k]*widths[k] for  k in range(len(Density))]

    return float(sum(cum_dens)) / num_res


def density2(Distances, mutation_cts, MinDist=0, MaxDist=10):
    '''
    For a given residues center of geometry, count the number of other residue
    center of geometries between a min and max disatance
    '''
    #Counts = []

    #for Distance in Distances:
       #Counts.append(floor(Distance))

    #Density = []
    #for Distance, Count in sorted(Counter(Counts).items(), key=lambda Count: Count[0]):
        #if Distance in range(MinDist, MaxDist):
            #Density.append(Count)

    return sum(mutation_cts[d[1]] for d in Distances if MinDist <= int(d[0]) <= MaxDist)
