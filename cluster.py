#!/usr/bin/env python
import numpy as np
import subprocess,os,itertools,copy

class Cluster():
    def __init__(self):
        subprocess.check_call('getclus > clusters.tmp',shell=True)
        self.data=np.loadtxt('clusters.tmp')
        os.remove('clusters.tmp')

    def get_cluster_number(self):
        'get the number of all clusters'
        return len(self.data)
    
    def get_order_cluster_number(self,order):
        'get the number of a given order cluster'
        return list(self.data[:,0]).count(order)

    def get_exclude_cluster_number(self,order):
        'get the number of excluded clusters from front'
        exclude=0
        for i in range(order):
            exclude+=self.get_order_cluster_number(i)
        return exclude

    def get_cluster_index(self,order):
        'get the index of clusters for the given order'
        index=[]
        for i in range(len(self.data)):
            if self.data[i,0]==order:
                index.append(i)
        return index

    def get_cluster_multi(self):
        return self.data[:,2]

    def get_cluster_radius(self):
        index=range(self.get_exclude_cluster_number(3))
        radius=[0.0]
        for i in index:
            if abs(self.data[i,1]-radius[-1])>=0.01:
                radius.append(float(self.data[i,1]))
        return radius

    def get_cutoff_clusters(self,radius,order=6):
        'get the index of clusters within certain given cut-off radius'
        index=[]
        clus_num_2=0
        clus_num_3=0
        clus_num_4=0
        for i in range(len(self.data)):
            if self.data[i,0]<=order and self.data[i,1]<=radius+0.01:
                index.append(i)
                if self.data[i,0]==2:
                    clus_num_2+=1
                elif self.data[i,0]==3:
                    clus_num_3+=1
                elif self.data[i,0]==4:
                    clus_num_4+=1
        print clus_num_2,clus_num_3,clus_num_4
        return index

    def get_cluster_functions(self):
        'get cluster functions from every directory without error'
        cluster_function=[]
        clus_multi=self.get_cluster_multi()
        for item in os.listdir(os.getcwd()):
            fullpath=os.path.join(os.getcwd(),item)
            if os.path.isdir(fullpath):
                os.chdir(fullpath)
                if os.path.isfile('str.out') and not os.path.isfile('error'):
                    cluster_function.append(np.multiply(map(float,subprocess.check_output('corrdump -c -l=../lat.in -cf=../clusters.out',shell=True).split()),clus_multi))
                os.chdir('../')
        cluster_function=np.array(cluster_function,dtype='float')
        condition_number=np.linalg.cond(cluster_function)
        if condition_number >= 1E15:
            print "WARNING: The condition number of the matrix of cluster functions is too large: %.6e\n" %(np.linalg.cond(cluster_function))
        return np.array(cluster_function)

    def get_cluster_function_from_file(self):
        'get cluster function from allcorr.out'
        return np.loadtxt('allcorr.out')

    def generate_candidates(self,order,radius_cutoff=float('inf')):
        'use a hierarchical method to generate candidates of cluster models. Clusters of the same radius are added simultaneously.'
        candidates=[]
        index=[]
        yield index #the empty 
        for i in range(len(self.data)):
            if self.data[i,0]==order and self.data[i,1]<=radius_cutoff:
                if i not in index:
                    index.append(i)
                if i+1<len(self.data) and self.data[i,1]==self.data[i+1,1]:
                    index.append(i+1)
                else:
                    yield index

    def construct_candidates(self,max_order,max_cluster_number):
        'construct candidates'
        #for i in itertools.product(self.generate_candidates(2),self.generate_candidates(3)):
        index=range(self.get_exclude_cluster_number(2)) #index of null and point clusters
        for i in self.generate_candidates(2):
            index_2=index+i
            for j in self.generate_candidates(3,self.data[index_2[-1],1]):
                index_3=index_2+j
                for k in self.generate_candidates(4,self.data[index_3[-1],1]):
                    index_4=index_3+k 
                    if len(index_4)<=max_cluster_number:
                        #print index_4
                        yield index_4
            #index=[0,1,2]
            #index=list(np.array(i).flatten())
            #index+=i
            #index+=j
            #index+=k
            #if len(index)<=max_cluster_number:
            #    print index
                #yield index
       
class Eci(Cluster):
    def init(self,name):
        try:
            self.eci = np.loadtxt('%s.eci' %(name))
        except:
            self.eci = None


#cluster=Cluster()
#a=cluster.generate_candidates(2)
#for i in a:
#    print i
#models=cluster.construct_candidates(5,80)
#for j in models:
#    print j
#two_body_2=cluster.generate_candidates(2)
