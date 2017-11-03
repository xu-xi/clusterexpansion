#!/usr/bin/env python
import numpy,argparse,subprocess,os,sys

def Main(Arglist):
    parser=argparse.ArgumentParser(description='Get cluster correlation functions from MC under certain temperature for SQS. Note that the number of clusters for SQS is in accordance with that of MC',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('-f',default='mc.out',dest='mcfile',help="the Monte Carlo data file")
    parser.add_argument('-o',dest='outfile',default='tcorr.out',help="output file")
    parser.add_argument('-T',type=int,dest='temp',required=True,help="temperature")
    args=parser.parse_args()

    subprocess.check_call('getclus > clusters.tmp',shell=True)
    cluster_number=len(file('clusters.tmp').readlines())
    os.remove('clusters.tmp')

    try:
        mcdata=numpy.loadtxt(args.mcfile)
    except IOError,error_msg1:
        try:
            mcdata=numpy.loadtxt('../%s' % (args.mcfile))
        except IOError,error_msg2:
            print error_msg1,'\n',error_msg2
            sys.exit(1)

    mc_temps=map(round,list(mcdata[:,0]))

    if args.temp not in mc_temps:
        print 'ERROR: %s is not available' % (args.temp)
        print mc_temps
        sys.exit(1)
    else:
        ccfs=mcdata[mc_temps.index(args.temp),-cluster_number:]
        numpy.savetxt(args.outfile,ccfs,fmt='%.6f')

if __name__ =="__main__":
    Main(sys.argv)
