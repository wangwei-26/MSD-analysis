#!/bin/python3
import scipy.optimize as optimization

def fourPL(x, A, B, C, D):
    return ((A-D)/(1.0+((x/C)**(B))) + D)
    ###b1=D, b2=A, b3=C, b4=B

def concentration(x,a,b,c,d):
    z=c*(((a-d)/(x-d)-1)**(1/b))
    return z

f=open('RP2_raw.txt') ###
info={}
tag='tag'
temp=[]
row=[]
for line in f:
    x=line.strip('\n').strip().strip('\t').split('\t')
    if len(x)==13:
        info[tag]=temp
        temp=[]
        tag=x[0].strip()
        row.append(tag)
        temp2=[]
        for i in range(1,len(x)):
            temp2.append(x[i])
        temp.append(temp2)
    else:
        temp.append(x)
f.close()
info[tag]=temp
info.pop('tag')

reads={}

for i in range(len(info[row[0]])):
    tag='Prot'+str(i+1)
    reads[tag]={}
    for j in row:
        reads[tag][j]=info[j][i]
layout={}

f=open('new_RP2_layout.txt') ###
pos={}
standard=[]
sample=[]
dilution={}
dilution_fold={}
dilution_fold['score']=[]
for line in f:
    x=line.strip('\n').strip().strip('\t').split('\t')
    if len(x)==13:
        for i in range(1,len(x)):
            tag=x[i].strip().strip('"')
            if tag not in pos.keys():
                pos[tag]=[]
            pos[tag].append(x[0]+'-'+str(i))
            if 'STANDARD' in tag.upper() and tag not in standard:
                standard.append(tag)
            else:
                if tag not in sample:
                    sample.append(tag)
                if ':' in x[i]:
                    b=x[i].strip().split(' ')[1].split(':')[1].strip('"').split(',')
                    c=''
                    for k in b:
                        c+=k
                    dilution[tag]=int(c.strip('"'))
                    if c.strip('"') not in dilution_fold.keys():
                        dilution_fold[c.strip('"')]={}
                        dilution_fold['score'].append(int(c.strip('"')))
                    dilution_fold[c.strip('"')][tag.split(' ')[0].strip()]=tag
dilution_fold['score'].sort()
                    
f.close()
standard_con={}

f=open('Standards_conc_K15659U.txt') ###
for line in f:
    x=line.strip('\n').strip('\t').strip().split('\t')
    temp=[]
    for i in range(1,len(x)):
        temp.append(x[i].strip())
    standard_con[x[0].strip()]=temp
f.close()

#potential bug for more than 2 replicates
reads_avg={}
for i in reads.keys():
    reads_avg[i]={}
    for j in pos.keys():
        x1=pos[j][0].split('-')
        x2=pos[j][0].split('-')
        score1=reads[i][x1[0]][int(x1[1])-1]
        score2=reads[i][x2[0]][int(x2[1])-1]
        reads_avg[i][j]=(int(score1)+int(score2))/2
##for i in reads_avg.keys():
##    print(i, reads_avg[i])

pl={}
prot=[]
for i in range(len(standard_con[standard[0]])):
    tag='Prot'+str(i+1)
    prot.append(tag)
    xdata=[]
    ydata=[]
    for j in standard_con.keys():
        xdata.append(float(standard_con[j][i]))
        ydata.append(reads_avg[tag][j])
    params, params_covariance = optimization.curve_fit(fourPL, xdata, ydata)
    pl[tag]=[list(params)]

score={}
score2={}
for i in pl.keys():
    score[i]={}
    score2[i]={}
    for j in standard:
        x=reads_avg[i][j]
        a=concentration(x,pl[i][0][0],pl[i][0][1],pl[i][0][2],pl[i][0][3])
        score[i][j]=[a]
        score2[i][j]=[a]
    for j in sample:
        x=reads_avg[i][j]
        a=concentration(x,pl[i][0][0],pl[i][0][1],pl[i][0][2],pl[i][0][3])
        score2[i][j]=[a]
        if j in dilution.keys():
            score[i][j]=[a*dilution[j]]
        else:
            score[i][j]=[a]

out=open('RP2_result.txt','w') ###
head='sample\t'

for i in prot:
    head+=i+'-readings\tconc.\t'+i+'-Cal.Con.\t'+i+'-1\t'
out.write(head.strip('\t')+'\n')
used={}
for i in standard:
    result=i+'\t'
    used[i]=1
    for j in prot:
        temp=int(j[4:])-1
        result+=str(reads_avg[j][i])+'\t'+str(standard_con[i][temp])+'\t'+str(score2[j][i][0])+'\t'+str(score[j][i][0])+'\t'
    out.write(result.strip('\t')+'\n')
for i in sample:
    if i not in used.keys():
        result=i+'\t'
        for j in prot:
            result+=str(reads_avg[j][i])+'\t-\t'+str(score2[j][i][0])+'\t'+str(score[j][i][0])+'\t'
        out.write(result.strip('\t')+'\n')
out.close()

out=open('RP2_parameters.txt','w') ###
out.write('prot\tA, B, C, D\n')

for i in pl.keys():
    result=i+'\t'
    for j in pl[i][0]:
        result+=str(j)+','
    result=result.strip(',')+'\n'
    out.write(result)
out.close()

rename={}
f=open('Analytes_K15659U.txt') ###
count=0
for line in f:
    x=line.strip('\n').strip()
    count+=1
    rename['Prot'+str(count)]=x
f.close()
order_list=[]

f=open('samplelist_SP2_RP2.txt') ###
for line in f:
    x=line.strip('\n').strip()
    order_list.append(x)
f.close()

f=open('RP2_result.txt') ###
out=open('output_RP2_result.txt','w') ###
for line in f:
    x=line.strip('\n').strip().split('\t')
    head=''
    for i in x:
        y=i.split('-')
        if y[0] not in rename.keys():
            head+=i+'\t'
        else:
            temp=rename[y[0]]+'-'
            for j in range(1,len(y)):
                temp+=y[j]+'-'
            head+=temp.strip('-')+'\t'
    out.write(head.strip('\t')+'\n')
    break
info2={}
for line in f:
    x=line.strip('\n').strip().split('\t')
    if 'STANDARD' in x[0].upper():
        out.write(line)
    else:
        info2[x[0]]=line
f.close()

for i in dilution_fold['score']:
    for j in order_list:
        if j in dilution_fold[str(i)].keys():
            out.write(info2[dilution_fold[str(i)][j]])
        else:
            print(j)
out.close()


                



    


