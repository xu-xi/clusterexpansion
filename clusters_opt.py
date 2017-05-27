#!/usr/bin/env python
import subprocess,numpy,argparse
from myfunc import read_clusters

def clusters_optimizer(max_cluster_number,property_to_expand,average=False):
    'optimal cluster expansion construction'
    pair_distance=0
    triplet_distance=0
    quad_distance=0
    cluster_number_list=[[0,0,0]]
    cv_list=[]

    subprocess.check_call('corrdump -clus -2=0.5',shell=True)
    if average==False:
        cv_list.append(float(subprocess.check_output('clusterexpand -e -cv %s | tail -1' %(property_to_expand),shell=True)))
    else:
        cv_list.append(float(subprocess.check_output('clusterexpand -e -pa -cv %s | tail -1' %(property_to_expand),shell=True)))
    #null_point_cluster_number=int(subprocess.check_output('getclus | wc -l',shell=True))

    while 1:
        pair_distance+=0.5
        for triplet_distance in numpy.arange(0,pair_distance,0.5):
            for quad_distance in numpy.arange(0,triplet_distance,0.5):
                subprocess.check_call('corrdump -clus -2=%s -3=%s -4=%s' %(pair_distance,triplet_distance,quad_distance),shell=True)
                clusters=[]
                for i in [2,3,4]:
                    a,b=read_clusters(i)
                    clusters.append(a)
                cluster_number=int(subprocess.check_output('getclus | wc -l',shell=True))
                if cluster_number>=max_cluster_number:
                    break
                if clusters not in cluster_number_list:
                    cluster_number_list.append(clusters)
                    #print cluster_number
                    #print pair_distance,triplet_distance,quad_distance
                    if average==False:
                        cv=float(subprocess.check_output('clusterexpand -e -cv %s | tail -1' %(property_to_expand),shell=True))
                    else:
                        cv=float(subprocess.check_output('clusterexpand -e -pa -cv %s | tail -1' %(property_to_expand),shell=True))
                    if cv<min(cv_list):
                        optimal_pair_distance=pair_distance
                        optimal_triplet_distance=triplet_distance
                        optimal_quad_distance=quad_distance
                    cv_list.append(cv)
            else:
                continue
            break
        else:
            continue
        #print cluster_number_list,cv_list
        subprocess.check_call('corrdump -clus -2=%s -3=%s -4=%s' %(optimal_pair_distance,optimal_triplet_distance,optimal_quad_distance),shell=True)
        if average==False:
            optimal_cv=float(subprocess.check_output('clusterexpand -e -cv %s | tail -1' %(property_to_expand),shell=True))
        else:
            optimal_cv=float(subprocess.check_output('clusterexpand -e -pa -cv %s | tail -1' %(property_to_expand),shell=True))
        return optimal_cv

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='optimal cluster expansion construction',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n',dest='max_clusters_number',type=int,help="the maximum number of clusters")
    parser.add_argument('-p',default='energy',dest='property',help="the property to expand")
    #parser.add_argument('-o',action='store_true',dest='opt_only',help="do not do cluster expansion, only optimize the choice of clusters")
    #parser.add_argument('--version',action='version',version='2017.2.23',help="output the version of the program")

    args=parser.parse_args()
    
    if args.max_clusters_number!=None:
        clusters_optimizer(args.max_clusters_number,args.property)
    else:
        structure_number=0
        #scan current directory to collect finished calculation
        for item in os.listdir(os.environ['PWD']):
            fullpath=os.path.join(os.environ['PWD'],item)
            if os.path.isdir(fullpath):
                os.chdir(fullpath)
                if os.path.isfile('energy') and not os.path.isfile('error'):
                    structure_number+=1
                    os.chdir('../')
        clusters_optimizer(structure_number,args.property)
