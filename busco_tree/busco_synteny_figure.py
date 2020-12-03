#! /usr/bin/env python
import itertools
import numpy as np

### BUSCO sampling approach code
threads=80

#parse data, pull out information
infile = open("complete_busco_locations.csv", "r")
lines=infile.readlines()
lines=[x.strip('\n').split(',') for x in lines]

data = lines[1:]
species, buscoIDs, contig, start, end = zip(*data)

#get unique species names
species = set(species)
buscoIDs = set(buscoIDs)

#store genomes and BUSCO sequences in nested dictionary
#top level key is species name
#second level key is contig name
genomes = {}
for sp in list(species):
    data_subset = [x for x in data if x[0] == sp]
    contig_names = list(set(list(zip(*data_subset))[2]))
    contigs = {}
    for contig in contig_names:
        genome_subset = [x for x in data_subset if x[2] == contig]
        genome_subset = [[sp,busco,contig,int(round(np.mean(list(map(int,[start,end])))))]
                         for sp,busco,contig,start,end in genome_subset]
        #sort by busco position on contig
        genome_subset = sorted(genome_subset, key = lambda x: x[3])
        contigs[contig] = [x[1] for x in genome_subset]
    genomes[sp] = contigs

#get flat contigs
flat_contigs = [[genomes[x][y] for y in genomes[x]] for x in genomes]
flat_contigs = [x for y in flat_contigs for x in y]

#get all n-wise neighbor combinations and store as sets
n=2 #size of chunks
k=1 #size of overlap
assert n > k
neighbors=[]
for contig in flat_contigs:
    #ignore stuff without neighbors
    if len(contig) <= (n):
        continue
    obs = [[contig[x:x+n], 
           contig[-(len(contig)-x+1):-(len(contig)-x-n+1)][::-1]] for 
           x in range(1, len(contig)-k, n-k)]
    neighbors.append([x for y in obs for x in y])

#flatten list of lists
neighbors = [x for y in neighbors for x in y]

#now for every unique BUSCO pair count their frequencies
#neighbor_counts=[]
#for neighbor in set(neighbors):
#    neighbor_counts.append([neighbor,neighbors.count(neighbor)])

# parallelize counting op
import multiprocessing as mp
def f(buscoID):
    pair_list = [x for x in neighbors if buscoID == x[0]]
    pair_list = ['{}to{}'.format(x[0],x[1]) for x in pair_list]
    uniq_pair_list = list(set(pair_list))
    n = len(uniq_pair_list)
    uniq_count_list = np.array([pair_list.count(x) for x in uniq_pair_list])
    uniq_prop_list = uniq_count_list/np.sum(uniq_count_list)
    return list(zip([buscoID]*n, uniq_pair_list, uniq_count_list,
                    uniq_prop_list))

p = mp.Pool(threads)
neighbor_counts = p.map(f, list(set(list(zip(*neighbors))[0])))
p.terminate()
p.join()

neighbor_counts = [x for y in neighbor_counts for x in y]

outfile = open('busco_pairs.csv','w')
outfile.write('busco,connection,count,proportion\n')
for busco, conn, count, prop in neighbor_counts:
    outline = '{},{},{},{}\n'.format(busco, conn, count, prop)
    bytes=outfile.write(outline)

outfile.close()

#gephi formatted nodes.csv: Id,Label,Attribute
# do later: get D.melanogaster chromosome element as attribute
outfile=open("nodes.csv","w")
outfile.write('Id,Label\n')
for busco in buscoIDs:
    bytes=outfile.write('{},{}\n'.format(busco,busco))

outfile.close()

#gephi formatted edges.csv: Source,Target,Type,Weight
# Weight = count, Type=Undirected
outfile=open("edges.csv","w")
bytes=outfile.write('Source,Target,Type,Weight\n')

for busco,connection,count,proportion in neighbor_counts:
    p1, p2 = connection.split("to")
    #outputs COUNT!
    bytes=outfile.write('{},{},Undirected,{}\n'.format(p1,p2,count))

outfile.close()

#designate BUSCO nodes with melanogaster chromosome
busco_locations={}
busco_locations_data = open('dmel_complete_buscos.csv','r').readlines()
busco_locations_data = [x.strip('\n').split(',') for x in busco_locations_data]
for key, val in busco_locations_data:
    busco_locations[key] = val

lines=open('nodes.csv','r').readlines()[1:]
lines=[x.strip('\n').split(',') for x in lines]

outfile=open("nodes_with_attrs.csv","w")
bytes=outfile.write('Id,Label,Attribute\n')

for ID,label in lines:
    if (ID in busco_locations):
        chrom="chr{}".format(busco_locations[ID])
    else:
        chrom="not present in D. mel"
    bytes=outfile.write('{},{},{}\n'.format(ID,label,chrom))

outfile.close()