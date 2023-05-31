#!/usr/local/bin/python3

from sys import argv
import numpy as np

def confusion_matrix(input_file,threshold):
    f=open(input_file)

    cm=np.zeros([2,2]).astype(int)

    for line in f:
        eval = float(line.split()[1])
        classK = int(line.split()[2])
        
        if eval <= threshold and classK == 1: #TRUE POSITIVE
            cm[0][0]+=1
        elif eval > threshold and classK == 0: #TRUE NEGATIVE
            cm[1][1]+=1
        elif eval <= threshold and classK == 0: #FALSE POSITIVE
            cm[0][1]+=1
        elif eval > threshold and classK == 1: #FALSE NEGATIVE
            cm[1][0]+=1    
    f.close()
    tp,tn,fp,fn = cm[0][0],cm[1][1],cm[0][1],cm[1][0]

    acc=(tp+tn)/(tp+tn+fn+fp)
    mcc=((tp*tn)-(fp*fn))/np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

    return acc,mcc,cm


if __name__ == '__main__':
    conf_file=argv[1]
    threshold=float(argv[2])
    acc,mcc,cm=confusion_matrix(conf_file,threshold)
    print('THR: %.1E\nACC: %f\nMCC: %f'%(threshold,acc,mcc)) 
    print('CM:',cm)    
