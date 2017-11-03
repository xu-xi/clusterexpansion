#!/usr/bin/env python
import subprocess,numpy,argparse,os
from myfunc import read_clusters

def clusters_optimizer(property_to_expand,max_cluster_number,scan_step,average=False):
    'optimal cluster expansion construction'
    pair_distance=0.5
    triplet_distance=0
    quad_distance=0
    optimal_pair_distance=pair_distance
    optimal_triplet_distance=triplet_distance
    optimal_quad_distance=quad_distance

    cluster_number_list=[[0,0,0]]
    cv_list=[]
    cv_file=file('cv.log','w')
    cv_file.write('#number CV pair triplet quad\n')

    subprocess.check_call('corrdump -clus -2=0.5',shell=True)
    cluster_number=int(subprocess.check_output('getclus | wc -l',shell=True))
    if average==False:
        cv_list.append(float(subprocess.check_output('clusterexpand -e -cv %s | tail -1' %(property_to_expand),shell=True)))
    else:
        cv_list.append(float(subprocess.check_output('clusterexpand -e -pa -cv %s | tail -1' %(property_to_expand),shell=True)))
    #null_point_cluster_number=int(subprocess.check_output('getclus | wc -l',shell=True))

    while 1:
        pair_distance+=scan_step
        subprocess.check_call('corrdump -clus -2=%s' %(pair_distance),shell=True)
        cluster_number=int(subprocess.check_output('getclus | wc -l',shell=True))
        if cluster_number>max_cluster_number:
            break
        for triplet_distance in numpy.arange(0,pair_distance,scan_step):
            for quad_distance in numpy.arange(0,triplet_distance,scan_step):
                subprocess.check_call('corrdump -clus -2=%s -3=%s -4=%s' %(pair_distance,triplet_distance,quad_distance),shell=True)
                #print pair_distance,triplet_distance,quad_distance
                clusters=[]
                for i in [2,3,4]:
                    a,b=read_clusters(i)
                    clusters.append(a)
                cluster_number=int(subprocess.check_output('getclus | wc -l',shell=True))
                #print cluster_number
                if cluster_number>max_cluster_number:
                    break
                if clusters not in cluster_number_list:
                    cluster_number_list.append(clusters)
                    #print cluster_number
                    #print pair_distance,triplet_distance,quad_distance
                    if average==False:
                        cv=float(subprocess.check_output('clusterexpand -e -cv %s | tail -1' %(property_to_expand),shell=True))
                    else:
                        cv=float(subprocess.check_output('clusterexpand -e -pa -cv %s | tail -1' %(property_to_expand),shell=True))
                    cv_file.write('%s %s %s %s %s\n' %(cluster_number,cv,pair_distance,triplet_distance,quad_distance))
                    if cv<min(cv_list):
                        optimal_pair_distance=pair_distance
                        optimal_triplet_distance=triplet_distance
                        optimal_quad_distance=quad_distance
                    cv_list.append(cv)
    
    cv_file.close()
    subprocess.check_call('corrdump -clus -2=%s -3=%s -4=%s' %(optimal_pair_distance,optimal_triplet_distance,optimal_quad_distance),shell=True)
    if average==False:
        optimal_cv=float(subprocess.check_output('clusterexpand -e -cv %s | tail -1' %(property_to_expand),shell=True))
    else:
        optimal_cv=float(subprocess.check_output('clusterexpand -e -pa -cv %s | tail -1' %(property_to_expand),shell=True))
    return optimal_cv

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='optimal cluster expansion construction',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n',dest='max_cluster_number',type=int,help="the maximum number of clusters")
    parser.add_argument('-a',action='store_true',dest='average',help="decide whether the property to be averaged")
    parser.add_argument('-p',default='energy',dest='property',help="the property to expand")
    parser.add_argument('-s',type=float,default=0.5,dest='step',help="the scan step of clusters")

    args=parser.parse_args()
    
    if args.max_cluster_number!=None:
        clusters_optimizer(args.property,args.max_cluster_number,args.step,args.average)
    else:
        structure_number=0
        #scan current directory to collect finished calculation
        for item in os.listdir(os.environ['PWD']):
            fullpath=os.path.join(os.environ['PWD'],item)
            if os.path.isdir(fullpath):
                os.chdir(fullpath)
                if os.path.isfile(args.property) and not os.path.isfile('error'):
                    structure_number+=1
                    os.chdir('../')
        clusters_optimizer(args.property,structure_number,args.step,args.average)
