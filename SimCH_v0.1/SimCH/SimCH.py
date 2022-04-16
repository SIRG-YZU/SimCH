#Package release: SimCH (version 0.1, April 2022)
#Authors: Gongming Wang, Lei Sun
#Email: sunlei(at)yzu.edu.cn
#Description: SimCH is a Python-based simulator for benchmarking and evaluating scRNA-seq data computational methods.
#License: SimCH is distributed under the GNU GENERAL PUBLIC (GPL) license (Version 3).
#Copyright (C) 2022 Yangzhou University

"""
SimCH completed 10 functions,
-a means NB simulation with parameters
-b means NB-zero simulation wi
-c nocopula means NB simulation without parameters(no copula)
-d means NB-zero simulation without th parameters
-c means NB simulation without p
-e means NB simulation of multiple garameters(copula)parameters(copula)
-d nocopula means NB-zero simulation without parameters(no copula)roups without parameters(copula)
-e nocopula means NB simulation of multiple groups without parameters(no copula)
-f means NB-zero simulation of multiple groups without parameters(copula)
-f nocopula means NB-zero simulation of multiple groups without parameters(no copula))

"""

def Remove0():
    """
    Read the file and remove the 0 gene expressed

    Parameters
    ----------
    path1:str
          File path
    param:list
          All parameters used in simulation
          
    Returns
    -------
    ls2:list
          The data of gene

    """

    f = open(path1,'r')
    if path1[-3:] == 'csv':
        reader = csv.reader(f)
        ls = list(reader)
    else:

        f1 = f.readlines()
        if len(f1) < 250:
            while 1:
                f1.extend(f1[1:])
                if len(f1) > 250:
                    break
        ls = []
        for i in f1:
            k = i.strip().split('\t')
            ls.append(k)
    if len(ls[0]) < len(ls[1]):
        ls[0].insert(0,' ')

    ls1 = []
    for i in range(len(ls)-1):
        if float(max(ls[i+1][1:])) == 0:
            ls1.append(i+1)

    param.append(len(ls)-1)
    global ls2
    ls2 = []
    for i in range(len(ls)):
        if not i in ls1:
            ls2.append(ls[i])
    f.close()

    
def nml():
    """
    Normalize data with size factor

    Parameters
    ----------
    ls2:list
          The data of gene
    
          
    Returns
    -------
    ls3:list
          The data of nomalized gene
    ls7:list
          Size factor

    """
    global ls3
    ls3 = []
    global ls7                              #size factor

    ls7 = []
    ls = []
    for i in range(len(ls2[0])-1):
        sum1 = 0
        for j in range(len(ls2)-1):
            sum1 += float(ls2[j+1][i+1])
        ls.append(sum1)
    ls1 = sorted(ls)
    if len(ls1)%2 == 0:
        k = (ls1[len(ls1)//2-1]+ls1[len(ls1)//2])/2
    else:
        k = ls1[len(ls1)//2]
    for i in ls:
        ls7.append(i/k)
        

    ls3.append(ls2[0])
    for i in range(len(ls2)-1):
        ls8 = []
        ls8.append(ls2[i+1][0])
        for j in range(len(ls2[i+1])-1):
            ls8.append(float(ls2[i+1][j+1])/ls7[j])
        ls3.append(ls8)
    param.append(len(ls3)-1)
    param.append(len(ls3[1])-1)


    

def sizef():
    """
    Estimate distributed of size factor

    Parameters
    ----------
    ls7:list
          Size factor
    
          
    Returns
    -------
    param[3]:list
          GMM3 parameters

    """
    
    import warnings
    warnings.filterwarnings("ignore")
    lsm = sorted(ls7)
    lsn = np.log(lsm)

    x = np.array(lsn)
    n, bins, patches = plt.hist(x, 50,density=1, alpha=0.75)
    binm = []
    nm = []
    for i in range(len(bins)-1):
        binm.append((bins[i]+bins[i+1])/2)
        nm.append(n[i])
    binm.extend(bins[:-1])
    nm.extend(n)
    binm.extend(bins[1:])
    nm.extend(n)

    def fund(x,k1,k2,k3,mu1,sigma1,mu2,sigma2,mu3,sigma3):
        return k1*(1/(((2*math.pi)**(1/2))*sigma1))*np.exp(-((x-mu1)**2)/(2*(sigma1**2)))+k2*(1/(((2*math.pi)**(1/2))*sigma2))*np.exp(-((x-mu2)**2)/(2*(sigma2**2)))+k3*(1/(((2*math.pi)**(1/2))*sigma3))*np.exp(-((x-mu3)**2)/(2*(sigma3**2)))

    
    
        
    

    x = np.array(binm)
    y = np.array(nm)
    
    param_bounds=([0,0,0,-5,0.01,-5,0.01,-5,0.01],[1,1,1,5,np.inf,5,np.inf,5,np.inf])
    poptk, pcovk = curve_fit(fund,x,y,bounds=param_bounds)
    
    k1 = poptk[0]
    k2 = poptk[1]
    k3 = poptk[2]
    mu1 = poptk[3]
    sigma1 = poptk[4]
    mu2 = poptk[5]
    sigma2 = poptk[6]
    mu3 = poptk[7]
    sigma3 = poptk[8]

    
    

    
    x1 = np.linspace(lsn[0],lsn[-1],5000,endpoint=True)

    y2 = k1*(1/(((2*math.pi)**(1/2))*sigma1))*np.exp(-((x1-mu1)**2)/(2*(sigma1**2)))+k2*(1/(((2*math.pi)**(1/2))*sigma2))*np.exp(-((x1-mu2)**2)/(2*(sigma2**2)))+k3*(1/(((2*math.pi)**(1/2))*sigma3))*np.exp(-((x1-mu3)**2)/(2*(sigma3**2)))


    

    #plt.plot(x1,y2)
    #plt.show()

    k4 = k1+k2+k3
    k1 = round(k1/k4,4)
    k2 = round(k2/k4,4)
    k3 = round(1-k1-k2,4)
    

    
    param.append([k1,k2,k3,mu1,sigma1,mu2,sigma2,mu3,sigma3])
    param.append([lsn[0],lsn[-1]])
    


def NBzero_mean_dispersion():
    """
    Estimate the mean, dispersion, p0 andε parameters of -b

    Parameters
    ----------
    ls3:list
          The data of nomalized gene
     
    Returns
    -------
    param[5:7],param[19:24]:list
          Parameters of -b
    """
    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')

    global ls4,ls4p,ls2,ls7,k,k_a
    ls4 = []
    ls4p = []
    for i in range(len(ls3)-1):
        ls4p.append(ls3[i+1][1:])
    for i in range(len(ls4p)):
        sum1 = 0
        k = 0
        for j in range(len(ls4p[i])):
            sum1 += ls4p[i][j]
            k += 1
        b = sum1/k
        ls4.append(b)

    ls = []
    for i in range(len(ls4)):
        ls.append((i,ls4[i]))
    ls.sort(key=lambda x:x[1])



    ls1 = []
    for i in range(len(ls)):
        ls1.append(ls[-(i+1)][0])
        if len(ls1) > len(ls4p)*0.05:
            break


    global lme,ldis
    lme = []
    ldis = []
    lsk1 = []
    lsk2 = []
    k1 = 0
    s = np.array(ls7)
    print('Estimate γ,ε...')
    k3 = len(ls1)//100
    if k3 < 1:
        k3 = len(ls1)/100
    kk = 0
    warnings.filterwarnings("ignore")
    for i in range(len(ls4p)):
        if i in ls1:
            kk += 1
            progress(int(kk/k3),width=50)
            m0 = float(np.mean(ls4p[i]))
            p0 = ls4p[i].count(0)/len(ls4p[i])
            v0 = float(np.var(ls4p[i]))
            d0 = v0/(m0**2)


            l1 = []
            l2 = []
            l3 = []
            l4 = []
            step1 = (1/(1-p0))**(0.2)
            for i1 in range(7):
                for j1 in range(5):
                    for k1 in range(5):
                        for k2 in range(10):
                            if i1 == 0:
                                l1.append([m0])
                                l2.append([0.6*d0+j1*0.1*d0])
                                l3.append([0.1+0.2*k1])
                                l4.append([-5+k2*0.56])
                            elif i1 == 1:
                                if (step1-1)*m0*0.1+m0 > 0:
                                    l1.append([(step1-1)*m0*0.1+m0])
                                    l2.append([0.6*d0+j1*0.1*d0])
                                    l3.append([0.1+0.2*k1])
                                    l4.append([-5+k2*0.56])
                            else:
                                l1.append([m0*(step1**(i1-1))])
                                l2.append([0.6*d0+j1*0.1*d0])
                                l3.append([0.1+0.2*k1])
                                l4.append([-5+k2*0.56])
            l1 = np.array(l1)
            l2 = np.array(l2)
            l3 = np.array(l3)
            l4 = np.array(np.exp(l4))
            d = l2
            m = l1
            p1 = l3
            k = l4
            g = np.exp(gammaln(1/d+1)-gammaln(1/d))
            s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
            if p0 > 0:
                s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
            g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
            m1 = np.transpose(np.array(s2.mean(axis=1)))
            ls = []
            for i1 in m1:
                ls.append([i1])
            m1  = np.array(ls)

            s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
            s4 = s3.mean(axis=1)
            d1 = []
            for i1 in s4:
                d1.append([i1])
            d1 = np.array(d1)
            d1 = d1-m1**2
            if p0 > 0:
                s5 = s1.mean(axis=1)
                ls = []
                for i1 in s5:
                    ls.append([i1])
                p1 = np.array(ls)

            if p0 > 0:

                if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                    l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                else:
                    l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

            else:
                l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            r = np.argmin(l)
            m = l1[r][0]
            d = l2[r][0]
            p1 = l3[r][0]
            k = np.log(l4[r][0])



            l1 = []
            l2 = []
            l3 = []
            l4 = []
            step1 = (1/(1-p0))**(0.2)
            for i1 in range(5):
                for j1 in range(5):
                    for k1 in range(5):
                        for k2 in range(5):
                            if p1+(k1-2)*0.1 >= 0 and p1+(k1-2)*0.1 <=1:
                                if m == m0:
                                    if (step1-1)*m0*0.1*i1*0.25+m0 > 0:
                                        l1.append([(step1-1)*m0*0.1*i1*0.25+m0])
                                        l2.append([d+(j1-2)*0.05*d0])
                                        l3.append([p1+(k1-2)*0.1])
                                        l4.append([k+(k2-2)*0.28])
                                elif m == (step1-1)*m0*0.1+m0:
                                    if (step1-1)*m0*0.1*i1*0.5+m0 > 0:
                                        l1.append([(step1-1)*m0*0.1*i1*0.5+m0])
                                        l2.append([d+(j1-2)*0.05*d0])
                                        l3.append([p1+(k1-2)*0.1])
                                        l4.append([k+(k2-2)*0.28])
                                else:
                                    l1.append([m*(step1**((i1-2)*0.5))])
                                    l2.append([d+(j1-2)*0.05*d0])
                                    l3.append([p1+(k1-2)*0.1])
                                    l4.append([k+(k2-2)*0.28])
            l1 = np.array(l1)
            l2 = np.array(l2)
            l3 = np.array(l3)
            l4 = np.array(np.exp(l4))
            d = l2
            m = l1
            p1 = l3
            k = l4
            g = np.exp(gammaln(1/d+1)-gammaln(1/d))
            s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
            if p0 > 0:
                s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
            g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
            m1 = np.transpose(np.array(s2.mean(axis=1)))
            ls = []
            for i1 in m1:
                ls.append([i1])
            m1  = np.array(ls)

            s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
            s4 = s3.mean(axis=1)
            d1 = []
            for i1 in s4:
                d1.append([i1])
            d1 = np.array(d1)
            d1 = d1-m1**2
            if p0 > 0:
                s5 = s1.mean(axis=1)
                ls = []
                for i1 in s5:
                    ls.append([i1])
                p1 = np.array(ls)

            if p0 > 0:

                if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                    l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                else:
                    l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

            else:
                l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            r = np.argmin(l)
            m = l1[r][0]
            d = l2[r][0]
            p1 = l3[r][0]
            k = np.log(l4[r][0])





            l1 = []
            l2 = []
            l3 = []
            l4 = []
            step1 = (1/(1-p0))**0.2
            for i1 in range(5):
                for j1 in range(5):
                    for k1 in range(5):
                        for k2 in range(5):
                            if p1+(k1-2)*0.1 >= 0 and p1+(k1-2)*0.1 <=1:
                                if m == m0:
                                    if (step1-1)*m0*0.1*i1*0.125+m0 > 0:
                                        l1.append([(step1-1)*m0*0.1*i1*0.125+m0])
                                        l2.append([d+(j1-2)*0.025*d0])
                                        l3.append([p1+(k1-2)*0.05])
                                        l4.append([k+(k2-2)*0.14])
                                elif m <= (step1-1)*m0*0.2+m0:
                                    if (step1-1)*m0*0.05*(i1-2)*0.25+m > 0:
                                        l1.append([(step1-1)*m0*0.05*(i1-2)*0.25+m])
                                        l2.append([d+(j1-2)*0.025*d0])
                                        l3.append([p1+(k1-2)*0.05])
                                        l4.append([k+(k2-2)*0.14])
                                else:
                                    l1.append([m*(step1**((i1-2)*0.25))])
                                    l2.append([d+(j1-2)*0.025*d0])
                                    l3.append([p1+(k1-2)*0.05])
                                    l4.append([k+(k2-2)*0.14])
            l1 = np.array(l1)
            l2 = np.array(l2)
            l3 = np.array(l3)
            l4 = np.array(np.exp(l4))
            d = l2
            m = l1
            p1 = l3
            k = l4
            g = np.exp(gammaln(1/d+1)-gammaln(1/d))
            s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
            if p0 > 0:
                s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
            g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
            m1 = np.transpose(np.array(s2.mean(axis=1)))
            ls = []
            for i1 in m1:
                ls.append([i1])
            m1  = np.array(ls)

            s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
            s4 = s3.mean(axis=1)
            d1 = []
            for i1 in s4:
                d1.append([i1])
            d1 = np.array(d1)
            d1 = d1-m1**2
            if p0 > 0:
                s5 = s1.mean(axis=1)
                ls = []
                for i1 in s5:
                    ls.append([i1])
                p1 = np.array(ls)

            if p0 > 0:

                if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                    l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                else:
                    l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

            else:
                l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            r = np.argmin(l)
            m = l1[r][0]
            d = l2[r][0]
            p1 = l3[r][0]
            k = np.log(l4[r][0])



            l1 = []
            l2 = []
            l3 = []
            l4 = []
            step1 = (1/(1-p0))**0.2
            for i1 in range(5):
                for j1 in range(5):
                    for k1 in range(5):
                        for k2 in range(5):
                            if p1+(k1-2)*0.1 >= 0 and p1+(k1-2)*0.1 <=1:
                                if m == m0:
                                    if (step1-1)*m0*0.1*i1*0.0625+m0 > 0:
                                        l1.append([(step1-1)*m0*0.1*i1*0.0625+m0])
                                        l2.append([d+(j1-2)*0.0125*d0])
                                        l3.append([p1+(k1-2)*0.025])
                                        l4.append([k+(k2-2)*0.07])
                                elif m <= (step1-1)*m0*0.2+m0:
                                    if (step1-1)*m0*0.05*(i1-2)*0.125+m > 0:
                                        l1.append([(step1-1)*m0*0.05*(i1-2)*0.125+m])
                                        l2.append([d+(j1-2)*0.0125*d0])
                                        l3.append([p1+(k1-2)*0.025])
                                        l4.append([k+(k2-2)*0.07])
                                else:
                                    l1.append([m*(step1**((i1-2)*0.125))])
                                    l2.append([d+(j1-2)*0.0125*d0])
                                    l3.append([p1+(k1-2)*0.025])
                                    l4.append([k+(k2-2)*0.07])
            l1 = np.array(l1)
            l2 = np.array(l2)
            l3 = np.array(l3)
            l4 = np.array(np.exp(l4))
            d = l2
            m = l1
            p1 = l3
            k = l4
            g = np.exp(gammaln(1/d+1)-gammaln(1/d))
            s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
            if p0 > 0:
                s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
            g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
            m1 = np.transpose(np.array(s2.mean(axis=1)))
            ls = []
            for i1 in m1:
                ls.append([i1])
            m1  = np.array(ls)

            s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
            s4 = s3.mean(axis=1)
            d1 = []
            for i1 in s4:
                d1.append([i1])
            d1 = np.array(d1)
            d1 = d1-m1**2
            if p0 > 0:
                s5 = s1.mean(axis=1)
                ls = []
                for i1 in s5:
                    ls.append([i1])
                p1 = np.array(ls)

            if p0 > 0:

                if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                    l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                else:
                    l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

            else:
                l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            r = np.argmin(l)












            M = abs(m0-m1)[r]
            D = (abs(v0-d1)/v0)[r]
            if p0 > 0:
                P = (abs(p0-p1)/p0)[r]
            else:
                P = 0.05
            r_1 = r
            m1_1 = m1
            d1_1 = d1
            p1_1 = p1

            m = l1[r][0]
            d = l2[r][0]
            p1 = l3[r][0]
            k = l4[r][0]


            L = min(l)[0]
            L1 = L


            m_result = m
            d_result = d
            p1_result = p1
            k_result = k


            for i2 in range(20):

                if D < 0.05:
                    if M > 0:
                        l1 = np.random.normal(0,M,64)
                        l2 = np.random.normal(0,0.025*d_result,64)
                        if P*p0 > 0:
                            l3 = np.random.normal(0,P*p0,64)
                        else:
                            l3 = np.random.normal(0,p1_result*0.001,64)
                        if (m_result*p1_result) <= 0 or P*p0 <= 0:
                            l4 = np.random.normal(0,k_result*0.1,64)
                        else:
                            if P*p0/(m_result*p1_result) > 0.001 and P*p0/(m_result*p1_result) < 10:
                                l4 = np.random.normal(0,P*p0/(m_result*p1_result),64)
                            else:
                                l4 = np.random.normal(0,k_result*0.1,64)
                    else:
                        l1 = np.random.normal(0,m_result*0.001,64)
                        l2 = np.random.normal(0,0.025*d_result,64)
                        if P*p0 > 0:
                            l3 = np.random.normal(0,P*p0,64)
                        else:
                            l3 = np.random.normal(0,p1_result*0.001,64)
                        if (m_result*p1_result) <= 0 or P*p0 <= 0:
                            l4 = np.random.normal(0,k_result*0.1,64)
                        else:
                            if P*p0/(m_result*p1_result) > 0.001 and P*p0/(m_result*p1_result) < 10:
                                l4 = np.random.normal(0,P*p0/(m_result*p1_result),64)
                            else:
                                l4 = np.random.normal(0,k_result*0.1,64)
                else:
                    l1 = np.random.normal(0,0.01*m_result,64)
                    l2 = np.random.normal(0,D*d_result,64)
                    if P*p0 > 0:
                        l3 = np.random.normal(0,P*p0,64)
                    else:
                        l3 = np.random.normal(0,p1_result*0.001,64)
                    if (m_result*p1_result) <= 0 or P*p0 <= 0:
                        l4 = np.random.normal(0,k_result*0.1,64)
                    else:
                        if P*p0/(m_result*p1_result) > 0.001 and P*p0/(m_result*p1_result) < 10:
                            l4 = np.random.normal(0,P*p0/(m_result*p1_result),64)
                        else:
                            l4 = np.random.normal(0,k_result*0.1,64)


                l1a = []
                l2a = []
                l3a = []
                l4a = []

                for j in range(64):
                    if m_result+l1[j] > 0 and d_result+l2[j] > 0 and p1_result+l3[j] > 0 and p1_result+l3[j] <= 1 and k_result+l4[j] > 0:
                        l1a.append([m_result+l1[j]])
                        l2a.append([d_result+l2[j]])
                        l3a.append([p1_result+l3[j]])
                        l4a.append([k_result+l4[j]])
                l1 = np.array(l1a)
                l2 = np.array(l2a)
                l3 = np.array(l3a)
                l4 = np.array(l4a)

                if len(l1a) > 0:

                    d = l2
                    m = l1
                    p1 = l3
                    k = l4
                    g = np.exp(gammaln(1/d+1)-gammaln(1/d))


                    s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))


                    if p0 > 0:
                        s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
                    g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
                    m1 = np.transpose(np.array(s2.mean(axis=1)))
                    ls = []
                    for i1 in m1:
                        ls.append([i1])
                    m1  = np.array(ls)

                    s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
                    s4 = s3.mean(axis=1)
                    d1 = []
                    for i1 in s4:
                        d1.append([i1])
                    d1 = np.array(d1)
                    d1 = d1-m1**2
                    if p0 > 0:
                        s5 = s1.mean(axis=1)
                        ls = []
                        for i1 in s5:
                            ls.append([i1])
                        p1 = np.array(ls)



                    if p0 > 0:

                        if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                            l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                        else:
                            l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

                    else:
                        l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)



                    if min(l)[0] < L:
                        r = np.argmin(l)
                        da = d1[r]
                        pa = p1[r]
                        ma = m1[r]

                        M = abs(m0-m1)[r]
                        D = (abs(v0-d1)/v0)[r]
                        if p0 > 0:
                            P = (abs(p0-p1)/p0)[r]
                        else:
                            P = 0.05
                        m_result = l1[r][0]
                        d_result = l2[r][0]
                        p1_result = l3[r][0]
                        k_result = l4[r][0]
                        L = min(l)[0]

            if p1_result > 0:
                lsk1.append(p1_result)
                lsk2.append(k_result)
    k = np.median(lsk2)
    p1 = np.median(lsk1)
    print(' ')
    print('γ:{},ε:{}'.format(p1,k))

    param.append([p1,k])

    lme = []
    ldis = []
    k1 = 0
    k3 = len(ls4p)//100
    s = np.array(ls7)
    print('Estimate mean and dispersion...')
    for i in range(len(ls4p)):
        if i%k3 == 0:
            progress(int(i/k3),width=50)

        m0 = float(np.mean(ls4p[i]))
        p0 = ls4p[i].count(0)/len(ls4p[i])
        v0 = float(np.var(ls4p[i]))
        d0 = v0/(m0**2)


        l1 = []
        l2 = []
        step1 = (1/(1-p0))**(0.2)
        for i1 in range(7):
            for j1 in range(5):
                if i1 == 0:
                    l1.append([m0])
                    l2.append([0.6*d0+j1*0.1*d0])
                elif i1 == 1:
                    if (step1-1)*m0*0.1+m0 > 0:
                        l1.append([(step1-1)*m0*0.1+m0])
                        l2.append([0.6*d0+j1*0.1*d0])
                else:
                    l1.append([m0*(step1**(i1-1))])
                    l2.append([0.6*d0+j1*0.1*d0])
        l1 = np.array(l1)
        l2 = np.array(l2)
        d = l2
        m = l1
        g = np.exp(gammaln(1/d+1)-gammaln(1/d))
        s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
        if p0 > 0:
            s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
        g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
        m1 = np.transpose(np.array(s2.mean(axis=1)))
        ls = []
        for i1 in m1:
            ls.append([i1])
        m1  = np.array(ls)

        s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
        s4 = s3.mean(axis=1)
        d1 = []
        for i1 in s4:
            d1.append([i1])
        d1 = np.array(d1)
        d1 = d1-m1**2
        if p0 > 0:
            s5 = s1.mean(axis=1)
            ls = []
            for i1 in s5:
                ls.append([i1])
            p_1 = np.array(ls)

        if p0 > 0:

            if (abs(p0-p_1)/p0)[np.argmin((abs(p0-p_1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                l = (abs(p0-p_1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            else:
                l = (abs(p0-p_1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

        else:
            l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
        r = np.argmin(l)
        m = l1[r][0]
        d = l2[r][0]




        l1 = []
        l2 = []
        step1 = (1/(1-p0))**(0.2)
        for i1 in range(5):
            for j1 in range(5):
                if m == m0:
                    if (step1-1)*m0*0.1*i1*0.25+m0 > 0:
                        l1.append([(step1-1)*m0*0.1*i1*0.25+m0])
                        l2.append([d+(j1-2)*0.05*d0])
                elif m == (step1-1)*m0*0.1+m0:
                    if (step1-1)*m0*0.1*i1*0.5+m0 > 0:
                        l1.append([(step1-1)*m0*0.1*i1*0.5+m0])
                        l2.append([d+(j1-2)*0.05*d0])
                else:
                    l1.append([m*(step1**((i1-2)*0.5))])
                    l2.append([d+(j1-2)*0.05*d0])

        l1 = np.array(l1)
        l2 = np.array(l2)
        d = l2
        m = l1
        g = np.exp(gammaln(1/d+1)-gammaln(1/d))
        s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
        if p0 > 0:
            s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
        g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
        m1 = np.transpose(np.array(s2.mean(axis=1)))
        ls = []
        for i1 in m1:
            ls.append([i1])
        m1  = np.array(ls)

        s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
        s4 = s3.mean(axis=1)
        d1 = []
        for i1 in s4:
            d1.append([i1])
        d1 = np.array(d1)
        d1 = d1-m1**2
        if p0 > 0:
            s5 = s1.mean(axis=1)
            ls = []
            for i1 in s5:
                ls.append([i1])
            p_1 = np.array(ls)

        if p0 > 0:

            if (abs(p0-p_1)/p0)[np.argmin((abs(p0-p_1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                l = (abs(p0-p_1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            else:
                l = (abs(p0-p_1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

        else:
            l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
        r = np.argmin(l)
        m = l1[r][0]
        d = l2[r][0]






        l1 = []
        l2 = []
        l3 = []
        l4 = []
        step1 = (1/(1-p0))**0.2
        for i1 in range(5):
            for j1 in range(5):
                if m == m0:
                    if (step1-1)*m0*0.1*i1*0.125+m0 > 0:
                        l1.append([(step1-1)*m0*0.1*i1*0.125+m0])
                        l2.append([d+(j1-2)*0.025*d0])

                elif m <= (step1-1)*m0*0.2+m0:
                    if (step1-1)*m0*0.05*(i1-2)*0.25+m > 0:
                        l1.append([(step1-1)*m0*0.05*(i1-2)*0.25+m])
                        l2.append([d+(j1-2)*0.025*d0])

                else:
                    l1.append([m*(step1**((i1-2)*0.25))])
                    l2.append([d+(j1-2)*0.025*d0])

        l1 = np.array(l1)
        l2 = np.array(l2)

        d = l2
        m = l1

        g = np.exp(gammaln(1/d+1)-gammaln(1/d))
        s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
        if p0 > 0:
            s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
        g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
        m1 = np.transpose(np.array(s2.mean(axis=1)))
        ls = []
        for i1 in m1:
            ls.append([i1])
        m1  = np.array(ls)

        s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
        s4 = s3.mean(axis=1)
        d1 = []
        for i1 in s4:
            d1.append([i1])
        d1 = np.array(d1)
        d1 = d1-m1**2
        if p0 > 0:
            s5 = s1.mean(axis=1)
            ls = []
            for i1 in s5:
                ls.append([i1])
            p_1 = np.array(ls)

        if p0 > 0:

            if (abs(p0-p_1)/p0)[np.argmin((abs(p0-p_1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                l = (abs(p0-p_1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            else:
                l = (abs(p0-p_1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

        else:
            l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
        r = np.argmin(l)
        m = l1[r][0]
        d = l2[r][0]




        l1 = []
        l2 = []
        l3 = []
        l4 = []
        step1 = (1/(1-p0))**0.2
        for i1 in range(5):
            for j1 in range(5):
                if m == m0:
                    if (step1-1)*m0*0.1*i1*0.0625+m0 > 0:
                        l1.append([(step1-1)*m0*0.1*i1*0.0625+m0])
                        l2.append([d+(j1-2)*0.0125*d0])
                elif m <= (step1-1)*m0*0.2+m0:
                    if (step1-1)*m0*0.05*(i1-2)*0.125+m > 0:
                        l1.append([(step1-1)*m0*0.05*(i1-2)*0.125+m])
                        l2.append([d+(j1-2)*0.0125*d0])
                else:
                    l1.append([m*(step1**((i1-2)*0.125))])
                    l2.append([d+(j1-2)*0.0125*d0])
        l1 = np.array(l1)
        l2 = np.array(l2)
        d = l2
        m = l1
        g = np.exp(gammaln(1/d+1)-gammaln(1/d))
        s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
        if p0 > 0:
            s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
        g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
        m1 = np.transpose(np.array(s2.mean(axis=1)))
        ls = []
        for i1 in m1:
            ls.append([i1])
        m1  = np.array(ls)

        s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
        s4 = s3.mean(axis=1)
        d1 = []
        for i1 in s4:
            d1.append([i1])
        d1 = np.array(d1)
        d1 = d1-m1**2
        if p0 > 0:
            s5 = s1.mean(axis=1)
            ls = []
            for i1 in s5:
                ls.append([i1])
            p_1 = np.array(ls)

        if p0 > 0:

            if (abs(p0-p_1)/p0)[np.argmin((abs(p0-p_1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                l = (abs(p0-p_1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            else:
                l = (abs(p0-p_1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

        else:
            l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
        r = np.argmin(l)
        M = abs(m0-m1)[r]
        D = (abs(v0-d1)/v0)[r]
        m = l1[r][0]
        d = l2[r][0]
        L = min(l)[0]
        m_result = m
        d_result = d


        for i2 in range(20):

            if D < 0.05:
                if M > 0:
                    l1 = np.random.normal(0,M,64)
                    l2 = np.random.normal(0,0.025*d_result,64)

                else:
                    l1 = np.random.normal(0,m_result*0.001,64)
                    l2 = np.random.normal(0,0.025*d_result,64)
            else:
                l1 = np.random.normal(0,0.01*m_result,64)
                l2 = np.random.normal(0,D*d_result,64)
            l1a = []
            l2a = []


            for j in range(64):
                if m_result+l1[j] > 0 and d_result+l2[j] > 0:
                    l1a.append([m_result+l1[j]])
                    l2a.append([d_result+l2[j]])
            l1 = np.array(l1a)
            l2 = np.array(l2a)


            if len(l1a) > 0:

                d = l2
                m = l1
                g = np.exp(gammaln(1/d+1)-gammaln(1/d))


                s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))


                if p0 > 0:
                    s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
                g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
                m1 = np.transpose(np.array(s2.mean(axis=1)))
                ls = []
                for i1 in m1:
                    ls.append([i1])
                m1  = np.array(ls)

                s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
                s4 = s3.mean(axis=1)
                d1 = []
                for i1 in s4:
                    d1.append([i1])
                d1 = np.array(d1)
                d1 = d1-m1**2
                if p0 > 0:
                    s5 = s1.mean(axis=1)
                    ls = []
                    for i1 in s5:
                        ls.append([i1])
                    p_1 = np.array(ls)



                if p0 > 0:

                    if (abs(p0-p_1)/p0)[np.argmin((abs(p0-p_1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                        l = (abs(p0-p_1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                    else:
                        l = (abs(p0-p_1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

                else:
                    l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

                if min(l)[0] < L:
                    r = np.argmin(l)
                    M = abs(m0-m1)[r]
                    D = (abs(v0-d1)/v0)[r]
                    m_result = l1[r][0]
                    d_result = l2[r][0]
                    L = min(l)[0]
        lme.append(m_result)
        ldis.append(d_result)

    progress(100,width=50)
    print(' ')
    lsn = list(np.log(lme))

    ls1a = sorted(lsn)

    x = np.array(ls1a)
    n, bins, patches = plt.hist(x, 100,density=1, alpha=0.75)
    binm = []
    nm = []
    for i in range(len(bins)-1):
        binm.append((bins[i]+bins[i+1])/2)
        nm.append(n[i])
        
    binm.extend(bins[:-1])
    nm.extend(n)
    binm.extend(bins[1:])
    nm.extend(n)

    def fund(x,k1,k2,k3,k4,k5,mu1,sigma1,mu2,sigma2,mu3,sigma3,mu4,sigma4,mu5,sigma5):
        return k1*(1/(((2*math.pi)**(1/2))*sigma1))*np.exp(-((x-mu1)**2)/(2*(sigma1**2)))+k2*(1/(((2*math.pi)**(1/2))*sigma2))*np.exp(-((x-mu2)**2)/(2*(sigma2**2)))+k3*(1/(((2*math.pi)**(1/2))*sigma3))*np.exp(-((x-mu3)**2)/(2*(sigma3**2)))+k4*(1/(((2*math.pi)**(1/2))*sigma4))*np.exp(-((x-mu4)**2)/(2*(sigma4**2)))+k5*(1/(((2*math.pi)**(1/2))*sigma5))*np.exp(-((x-mu5)**2)/(2*(sigma5**2)))
    
    x = np.array(binm+5*binm[:2])
    y = np.array(nm+5*nm[:2])
    
    param_bounds=([0,0,0,0,0,ls1a[0],0.02,ls1a[0],0.02,ls1a[0],0.02,ls1a[0],0.02,ls1a[0],0.02,],[1,1,1,1,1,ls1a[-1],np.inf,ls1a[-1],np.inf,ls1a[-1],np.inf,ls1a[-1],np.inf,ls1a[-1],np.inf])
    try:
        poptk, pcovk = curve_fit(fund,x,y,bounds=param_bounds)
    except:
        param_bounds=([0,0,0,0,0,ls1a[0],0.01,ls1a[0],0.01,ls1a[0],0.01,ls1a[0],0.01,ls1a[0],0.01,],[1,1,1,1,1,ls1a[-1],5,ls1a[-1],5,ls1a[-1],5,ls1a[-1],5,ls1a[-1],5])
        poptk, pcovk = curve_fit(fund,x,y,bounds=param_bounds)

    
    k1 = poptk[0]
    k2 = poptk[1]
    k3 = poptk[2]
    k4 = poptk[3]
    k5 = poptk[4]
    mu1 = poptk[5]
    sigma1 = poptk[6]
    mu2 = poptk[7]
    sigma2 = poptk[8]
    mu3 = poptk[9]
    sigma3 = poptk[10]
    mu4 = poptk[11]
    sigma4 = poptk[12]
    mu5 = poptk[13]
    sigma5 = poptk[14]

    
    
    
    x1 = np.linspace(ls1a[0],ls1a[-1],5000,endpoint=True)
    y1 = k1*(1/(((2*math.pi)**(1/2))*sigma1))*np.exp(-((x1-mu1)**2)/(2*(sigma1**2)))+k2*(1/(((2*math.pi)**(1/2))*sigma2))*np.exp(-((x1-mu2)**2)/(2*(sigma2**2)))+k3*(1/(((2*math.pi)**(1/2))*sigma3))*np.exp(-((x1-mu3)**2)/(2*(sigma3**2)))+k4*(1/(((2*math.pi)**(1/2))*sigma4))*np.exp(-((x1-mu4)**2)/(2*(sigma4**2)))+k5*(1/(((2*math.pi)**(1/2))*sigma5))*np.exp(-((x1-mu5)**2)/(2*(sigma5**2)))
    #plt.plot(x1,y1)
    #plt.title('GMM')
    #plt.show()
    k6 = k1+k2+k3+k4+k5
    k1 = round(k1/k6,4)
    k2 = round(k2/k6,4)
    k3 = round(k3/k6,4)
    k4 = round(k4/k6,4)
    k5 = round(1-k1-k2-k3-k4,4)
    param.append([round(ls1a[0],4),round(ls1a[-1],4),k1,k2,k3,k4,k5,round(mu1,4),round(sigma1,4),round(mu2,4),round(sigma2,4),round(mu3,4),round(sigma3,4),round(mu4,4),round(sigma4,4),round(mu5,4),round(sigma5,4)])
    param.extend([' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '])


    ls = sorted(lsn)
    l = (ls[-1]+0.001-ls[0])/20
    l1 = []
    l2 = []
    for i in range(20):
        l1.append([])
        l2.append([])
        
    for i in range(len(ls4)):
        l1[int((lsn[i]-ls[0])/l)].append(lsn[i])
        l2[int((lsn[i]-ls[0])/l)].append(ldis[i])
    l3 = []
    l4 = []
    for i in range(len(l1)):
        l2[i].sort()
        if  len(l1[i]) > 10:
            l3.append(float(np.mean(l1[i])))
            l4.append(l2[i][int(len(l2[i])/2)])



            
    def funb(x,t1,t2,t3,d0):
        return t1/(t2+np.exp(t3*x))+d0
    x = np.array(l3[3:]+5*l3[8:])
    y = np.array(l4[3:]+5*l4[8:])
    param_bounds=([0,0,0,0.95*min(l4[10:])],[np.inf,np.inf,np.inf,2*min(l4[10:])])
    popt, pcov = curve_fit(funb,x,y,bounds=param_bounds)     

    T1 = popt[0]
    T2 = popt[1]
    T3 = popt[2]
    D0 = popt[3]

    param.extend([round(T1,3),round(T2,3),round(T3,3),round(D0,3)])
    
    l5a = sorted(lsn)
    x1 = np.linspace(l5a[0],l5a[-1],600000,endpoint=True)
    y1 =T1/(T2+np.exp(T3*x1))+D0

    #plt.scatter(lsn,ldis,s = 0.5,c = 'b')
    #plt.plot(x1,y1,c = 'g')
    #plt.xlabel('ln(mean)')
    #plt.ylabel('dispersion')
    #plt.title('Dispersion curve fitting')
    #plt.show()



    def funa(x,n):
        return n*(((n*x)**(n/2-1))*np.exp(-(n*x)/2))/((2**(n/2))*gamma(n/2))
    

    ls = []
    for i in range(len(lsn)):
        if lsn[i] >= l3[4]:
            if ldis[i]/(T1/(T2+np.exp(T3*lsn[i]))+D0) < 50:
                ls.append(ldis[i]/(T1/(T2+np.exp(lsn[i]))+D0))
     
    x = np.array(ls)
    n, bins, patches = plt.hist(x,500,density=1, alpha=0.75)
    
    binm = []
    nm = []
    for i1 in range(len(bins)-1):
        binm.append((bins[i1]+bins[i1+1])/2)
        nm.append(n[i1])
    

    x = np.array(binm)
    y = np.array(nm)
    param_bounds=([0],[200])
    poptk, pcovk = curve_fit(funa,x,y,bounds=param_bounds)
    n = poptk[0]
    


    param.append(round(n,3))
    param.extend([' '])

    
    l5 = sorted(ls)
    x1 = np.linspace(l5[0],l5[-1],600000,endpoint=True)
    y1 =n*(((n*x1)**(n/2-1))*np.exp(-(n*x1)/2))/((2**(n/2))*gamma(n/2))
    
    #plt.title('Estimated degrees of freedom')
    #plt.plot(x1,y1)
    #plt.show()
    
    

    
def mean():
    """
    Estimate mean of each gene

    Parameters
    ----------
    ls3:list
          The data of nomalized gene

    Returns
    -------
    ls4p:list
          Normalized count



    """
    param.append(' ')
    global ls4,ls4p         
    ls4 = []
    ls4p = []
    for i in range(len(ls3)-1):
        ls4p.append(ls3[i+1][1:])
    for i in range(len(ls4p)):
        sum1 = 0
        k = 0
        for j in range(len(ls4p[i])):
            sum1 += ls4p[i][j]
            k += 1
        b = sum1/k
        ls4.append(b)


def fitmean():
    """
    Fit the model of mean expression

    Parameters
    ----------
    ls4:list
          Mean of each gene

    Returns
    -------
    param[6]:list
          Parameters of mean

    """
    ls1 = sorted(ls4)
    ls1a = np.log(ls1)

    x = np.array(ls1a)
    n, bins, patches = plt.hist(x, 100,density=1, alpha=0.75)
    binm = []
    nm = []
    for i in range(len(bins)-1):
        binm.append((bins[i]+bins[i+1])/2)
        nm.append(n[i])
        
    binm.extend(bins[:-1])
    nm.extend(n)
    binm.extend(bins[1:])
    nm.extend(n)




    def fund(x,k1,k2,k3,k4,k5,mu1,sigma1,mu2,sigma2,mu3,sigma3,mu4,sigma4,mu5,sigma5):
        return k1*(1/(((2*math.pi)**(1/2))*sigma1))*np.exp(-((x-mu1)**2)/(2*(sigma1**2)))+k2*(1/(((2*math.pi)**(1/2))*sigma2))*np.exp(-((x-mu2)**2)/(2*(sigma2**2)))+k3*(1/(((2*math.pi)**(1/2))*sigma3))*np.exp(-((x-mu3)**2)/(2*(sigma3**2)))+k4*(1/(((2*math.pi)**(1/2))*sigma4))*np.exp(-((x-mu4)**2)/(2*(sigma4**2)))+k5*(1/(((2*math.pi)**(1/2))*sigma5))*np.exp(-((x-mu5)**2)/(2*(sigma5**2)))
    
    x = np.array(binm)
    y = np.array(nm)

    
    param_bounds=([0,0,0,0,0,ls1a[0],0.01,ls1a[0],0.01,ls1a[0],0.01,ls1a[0],0.01,ls1a[0],0.01,],[1,1,1,1,1,ls1a[-1],np.inf,ls1a[-1],np.inf,ls1a[-1],np.inf,ls1a[-1],np.inf,ls1a[-1],np.inf])
    try:
        poptk, pcovk = curve_fit(fund,x,y,bounds=param_bounds)
    except:
        param_bounds=([0,0,0,0,0,ls1a[0],0.01,ls1a[0],0.01,ls1a[0],0.01,ls1a[0],0.01,ls1a[0],0.01,],[1,1,1,1,1,ls1a[-1],5,ls1a[-1],5,ls1a[-1],5,ls1a[-1],5,ls1a[-1],5])
        poptk, pcovk = curve_fit(fund,x,y,bounds=param_bounds)
    
    k1 = poptk[0]
    k2 = poptk[1]
    k3 = poptk[2]
    k4 = poptk[3]
    k5 = poptk[4]
    mu1 = poptk[5]
    sigma1 = poptk[6]
    mu2 = poptk[7]
    sigma2 = poptk[8]
    mu3 = poptk[9]
    sigma3 = poptk[10]
    mu4 = poptk[11]
    sigma4 = poptk[12]
    mu5 = poptk[13]
    sigma5 = poptk[14]

    
    
    
    

    
    x1 = np.linspace(ls1a[0],ls1a[-1],5000,endpoint=True)
    y1 = k1*(1/(((2*math.pi)**(1/2))*sigma1))*np.exp(-((x1-mu1)**2)/(2*(sigma1**2)))+k2*(1/(((2*math.pi)**(1/2))*sigma2))*np.exp(-((x1-mu2)**2)/(2*(sigma2**2)))+k3*(1/(((2*math.pi)**(1/2))*sigma3))*np.exp(-((x1-mu3)**2)/(2*(sigma3**2)))+k4*(1/(((2*math.pi)**(1/2))*sigma4))*np.exp(-((x1-mu4)**2)/(2*(sigma4**2)))+k5*(1/(((2*math.pi)**(1/2))*sigma5))*np.exp(-((x1-mu5)**2)/(2*(sigma5**2)))
    #plt.plot(x1,y1)
    #plt.title('GMM')
    #plt.show()
    k6 = k1+k2+k3+k4+k5
    k1 = round(k1/k6,4)
    k2 = round(k2/k6,4)
    k3 = round(k3/k6,4)
    k4 = round(k4/k6,4)
    k5 = round(1-k1-k2-k3-k4,4)
    param.append([round(ls1a[0],4),round(ls1a[-1],4),k1,k2,k3,k4,k5,round(mu1,4),round(sigma1,4),round(mu2,4),round(sigma2,4),round(mu3,4),round(sigma3,4),round(mu4,4),round(sigma4,4),round(mu5,4),round(sigma5,4)])
    param.extend([' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '])








def BCV():
    """
    Fit the model of BCV

    Parameters
    ----------
    ls4p:list
          Modified gene expression data
    ls12:list
          BCV²of modified gene expression data

    Returns
    -------
    param[19:24]:list
          Relationship parameter of mean and dispersion and degree of freedom
    """

    from scipy.special import gammaln
    from scipy.special import digamma

    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')


    global ls2,ls3,ls7,ls4p,ls4
    ls12 = []
    for i in range(len(ls4p)):
        ls12.append(float(np.var(ls4p[i])/(ls4[i]**2)-1/ls4[i]))
    ls12a = []
    ls12b = []
    ls12d = []
    i1 = 0
    i2 = 0



    L = len(ls4)//100
    for i in range(len(ls2)-1):
        if i%L == 0:
            if i//L < 100:
                progress(i//L,width=50)
            else:
                progress(100,width=50)

        i2 = 0
        i3 = 0
        i4 = 0
        for j in range(len(ls2[i+1])-1):
            i2 += (ls7[j]*ls4[i])**2
            i3 += ls7[j]
            i4 += (float(ls2[i+1][j+1])-ls7[j]*ls4[i])**2
        da = (i4-ls4[i]*i3)/i2

        ls12b.append(da)

        
        if ls12[i] <= 0:
            d1 = 100
        else:
            d1 = ls12[i]
        mu = ls4[i]


        #l = list(np.random.normal(0,d1/3,15))
        #l1 = list(np.random.normal(0,d1/9,15))
        #l.extend(l1)

        #k = 0
        #for j in range(len(ls2[i+1])-1):
        #    if int(ls2[i+1][j+1]) == 0:
        #        k += 1
        
        #k1 = 0
        #for j in range(len(ls7)):
        #    k1 += ((d1**(-1)/(ls7[j]*mu+d1**(-1)))**(1/d1))
        #k2 = abs(k1-k)
        #d = d1
        #for j in l:
        #    if d+j > 0:
        #        d2 = d+j
        #    k1 = 0
        #    for j in range(len(ls7)):
        #        k1 += ((d2**(-1)/(ls7[j]*mu+d2**(-1)))**(1/d2))
        #    if abs(k1-k) < k2:
        #        d = d2
        #        k2 = abs(k1-k)
        #ls12d.append(d)


        


        
        
        la = []
        lb = []
        lc = []
        ld = []
        for j in range(len(ls7)):
            la.append(mu*ls7[j]*d1)
            lb.append(1/(mu*ls7[j]*d1)+1)
            lc.append(1/d1+ls4p[i][j])
            ld.append(1/d1)
        la1 = list(np.log(la))
        lb1 = list(np.log(lb))
        lc1 = list(digamma(lc))
        ld1 = list(digamma(ld))
        la1 = pd.Series(la1)
        lb1 = pd.Series(lb1)
        lc1 = pd.Series(lc1)
        ld1 = pd.Series(ld1)
        s1 = pd.Series(ls4p[i])
        s2 = pd.Series(ls7)
        #lnF = 0
        d12 = d1**2
        d13 = d1**3
        s = 1/d12*la1-1/d12+1/d12*lb1+(1+s1*d1)/(d12+mu*s2*d13)-1/d12*lc1+1/d12*ld1
        lnF = float(sum(s))

        #for j in range(len(ls7)):
        #    lnF += 1/d12*la1[j]-1/d12+1/d12*lb1[j]+(1+ls4p[i][j]*d1)/(d12+mu*ls7[j]*d13)-1/d12*lc1[j]+1/d12*ld1[j]

        if lnF < 0:
            K1 = 0
            K2 = d1
        else:
            K1 = d1
            K2 = 1000
        for j in range(10):
            d1 = (K1+K2)/2
            la = []
            lb = []
            lc = []
            ld = []
            for i1 in range(len(ls7)):
                la.append(mu*ls7[i1]*d1)
                lb.append(1/(mu*ls7[i1]*d1)+1)
                lc.append(1/d1+ls4p[i][i1])
                ld.append(1/d1)
            la1 = np.log(la)
            lb1 = np.log(lb)
            lc1 = digamma(lc)
            ld1 = digamma(ld)
            la1 = pd.Series(la1)
            lb1 = pd.Series(lb1)
            lc1 = pd.Series(lc1)
            ld1 = pd.Series(ld1)
            s1 = pd.Series(ls4p[i])
            s2 = pd.Series(ls7)
            #lnF = 0
            d12 = d1**2
            d13 = d1**3
            s = 1/d12*la1-1/d12+1/d12*lb1+(1+s1*d1)/(d12+mu*s2*d13)-1/d12*lc1+1/d12*ld1
            lnF = float(sum(s))
            #for i1 in range(len(ls7)):
            #    lnF += 1/d12*la1[i1]-1/d12+1/d12*lb1[i1]+(1+ls4p[i][i1]*d1)/(d12+mu*ls7[i1]*d13)-1/d12*lc1[i1]+1/d12*ld1[i1]
            if lnF > 0:
                K1 = d1
            else:
                K2 = d1
        
        ls12a.append((K1+K2)/2)
    progress(100,width=50)
    print(' ')
    ls4 = np.log(ls4)

    
        
    ls = sorted(ls4)
    l = (ls[-1]+0.001-ls[0])/20
    l1 = []
    l2 = []
    for i in range(20):
        l1.append([])
        l2.append([])
        
    for i in range(len(ls4)):
        l1[int((ls4[i]-ls[0])/l)].append(ls4[i])
        l2[int((ls4[i]-ls[0])/l)].append(ls12a[i])
    l3 = []
    l4 = []
    for i in range(len(l1)):
        l2[i].sort()
        if  len(l1[i]) > 10:
            l3.append(float(np.mean(l1[i])))
            l4.append(l2[i][int(len(l2[i])/2)])

    l5 = []
    l6 = []
    for i in range(20):
        l5.append([])
        l6.append([])
    for i in range(len(ls4)):
        l5[int((ls4[i]-ls[0])/l)].append(ls4[i])
        l6[int((ls4[i]-ls[0])/l)].append(ls12b[i])
    l7 = []
    l8 = []
    for i in range(len(l5)):
        l6[i].sort()
        if  len(l1[i]) > 10:
            l7.append(float(np.mean(l5[i])))
            l8.append(l6[i][int(len(l6[i])/2)])





            
    def funb(x,t1,t2,t3,d0):
        return t1/(t2+np.exp(t3*x))+d0


    ls1 = []
    for i in range(len(l4)):
        if l4[i] != 0:
            ls1.append(l8[i]/l4[i])
        else:
            ls1.append(0)
    x = np.array(l3[8:])
    y = np.array(ls1[8:])
    param_bounds=([0,0,-np.inf,0],[np.inf,np.inf,np.inf,np.inf])
    popt, pcov = curve_fit(funb,x,y,bounds=param_bounds)
    T1 = popt[0]
    T2 = popt[1]
    T3 = popt[2]
    D0 = popt[3]
    l4a = []
    for i in range(len(l3)):
        l4a.append((T1/(T2+np.exp(T3*l3[i]))+D0)*l4[i])
    ls12c = []
    for i in range(len(ls4)):
        ls12c.append((T1/(T2+np.exp(T3*ls4[i]))+D0)*ls12a[i])


    
    x = np.array(l3[8:])
    y = np.array(l4a[8:])
    param_bounds=([0,0,0,0.95*min(l4a[10:])],[np.inf,np.inf,np.inf,1.05*min(l4a[10:])])
    popt, pcov = curve_fit(funb,x,y,bounds=param_bounds)     

    T1 = popt[0]
    T2 = popt[1]
    T3 = popt[2]
    D0 = popt[3]

    param.extend([round(T1,3),round(T2,3),round(T3,3),round(D0,3)])
    
    l5a = sorted(ls4)
    x1 = np.linspace(l5a[0],l5a[-1],600000,endpoint=True)
    y1 =T1/(T2+np.exp(T3*x1))+D0

    #plt.scatter(ls4,ls12c,s = 0.5,c = 'b')      #Maximum likelihood
    #plt.scatter(ls4,ls12b,s = 0.5,c = 'pink')   #Moment
    #plt.scatter(ls4,ls12d,s = 0.5,c = 'y')      #0
    #plt.scatter(ls4,ls12,s = 0.5,c = 'y')
    #plt.scatter(l3,l4a,c = 'r')
    #plt.scatter(l7,l8,c = 'navy')
    #plt.plot(x1,y1,c = 'g')
    #plt.xlabel('ln(mean)')
    #plt.ylabel('dispersion')
    #plt.title('Dispersion curve fitting')
    #plt.show()
    


            
    
    def funa(x,n):
        return n*(((n*x)**(n/2-1))*np.exp(-(n*x)/2))/((2**(n/2))*gamma(n/2))
    def funb(x,n):
        return (1/((math.pi**0.5)*n))*np.exp(-x**2/(2*(n**2)))
    

    ls = []
    for i in range(len(ls4)):
        if ls12c[i]/(T1/(T2+np.exp(T3*ls4[i]))+D0) < 50:
            ls.append(ls12c[i]/(T1/(T2+np.exp(T3*ls4[i]))+D0))
     
    x = np.array(ls)
    n, bins, patches = plt.hist(x,200,density=1, alpha=0.75)
    
    binm = []
    nm = []
    for i1 in range(len(bins)-1):
        binm.append((bins[i1]+bins[i1+1])/2)
        nm.append(n[i1])
    

    x = np.array(binm)
    y = np.array(nm)
    param_bounds=([0],[200])
    poptk, pcovk = curve_fit(funa,x,y,bounds=param_bounds)
    n = poptk[0]
    

    param.append(round(n,3))
    param.extend([' '])

    
    l5 = sorted(ls)
    x1 = np.linspace(l5[0],l5[-1],600000,endpoint=True)
    y1 =n*(((n*x1)**(n/2-1))*np.exp(-(n*x1)/2))/((2**(n/2))*gamma(n/2))
    
    #plt.title('Estimated degrees of freedom')
    #plt.plot(x1,y1)
    #plt.show()
    





    ls2 = []
    ls3 = []
    ls4p = []
    ls7 = []


def createmean():
    """
    Simulate the mean of gene experession

    Parameters
    ----------
    param:list
         All parameters used in simulation
        
    Returns
    -------
    lmean:list
          Simulated gene mean
          
    """




    
    global lmean
    lmean = []

    
    lmean1 = []
    ls = []
    i = 0
    while 1:
        if i == int(param[1]*param[6][2]):
            break
        k = np.random.normal(param[6][7],param[6][8],1)
        if k < param[6][1] and k >= param[6][0]:
            ls.append(np.exp(k))
            i+=1
        
    lmean1 = list(map(float,ls))

    lmean2 = []
    ls = []
    i = 0
    while 1:
        if i == int(param[1]*param[6][3]):
            break
        k = np.random.normal(param[6][9],param[6][10],1)
        if k < param[6][1] and k >= param[6][0]:
            ls.append(np.exp(k))
            i+=1
        
    lmean2 = list(map(float,ls))


    lmean3 = []
    ls = []
    i = 0
    while 1:
        if i == int(param[1]*param[6][4]):
            break
        k = np.random.normal(param[6][11],param[6][12],1)
        if k < param[6][1] and k >= param[6][0]:
            ls.append(np.exp(k))
            i+=1
        
    lmean3 = list(map(float,ls))

    lmean4 = []
    ls = []
    i = 0
    while 1:
        if i == int(param[1]*param[6][5]):
            break
        k = np.random.normal(param[6][13],param[6][14],1)
        if k < param[6][1] and k >= param[6][0]:
            ls.append(np.exp(k))
            i+=1
        
    lmean4 = list(map(float,ls))

    lmean5 = []
    ls = []
    i = 0
    while 1:
        if i == param[1]-int(param[1]*param[6][2])-int(param[1]*param[6][3])-int(param[1]*param[6][4])-int(param[1]*param[6][5]):
            break
        k = np.random.normal(param[6][15],param[6][16],1)
        if k < param[6][1] and k >= param[6][0]:
            ls.append(np.exp(k))
            i+=1
        
    lmean5 = list(map(float,ls))

    lmean.extend(lmean1)
    lmean.extend(lmean2)
    lmean.extend(lmean3)
    lmean.extend(lmean4)
    lmean.extend(lmean5)

    lmean.sort()

    if mod == '-a':
        global ls4
        ls4.sort()
        if abs(max(ls4)/np.exp(param[6][1])-1) > 0.01:
            ls4 = list(np.exp(ls4))
        lmean[-int(len(lmean)*0.01):] = ls4[-int(len(lmean)*0.01):]
    elif mod == '-b':
        global lme
        lme.sort()
        lmean[-int(len(lmean)*0.01):] = lme[-int(len(lmean)*0.01):]
        

    
    
    for i in range(param[0]-param[1]):
        lmean.append(0)
    from random import shuffle
    shuffle(lmean)

    

def creatBCV():
    """
    Simulate the BCV of gene experession

    Parameters
    ----------
    param:list
          All parameters used in simulation
    
    Returns
    -------
    lbcv:list
          Simulated gene BCV
    lgenebcv2:list
          BCV² of each gene
          
    """
    global lbcv,lgenebcv2
    lgenebcv2 = []
    lbcv = []
    lsa = []
    for i in lmean:
        if i == 0:
            lsa.append(0)
        else:
            k1 = np.random.chisquare(param[23],1)/param[23]
            k = param[19]/(param[20]+np.exp(param[21]*np.log(i)))+param[22]
            lsa.append(k*k1)
    

    lgb = []
    lgb2 = []
    lgb3 = []
    if param[25] == [1]:
        for i in range(param[0]):
            ls = []
            for j in range(param[2]):
                ls.append(1)
            lgb3.append(ls)                                        
    else:
        global lgba,lba1
        lgba = []
        lba1 = []
        for i in range(param[2]):
            p1 = random.random()
            ls = [0]
            for j in range(len(param[25])):
                ls.append(sum(param[25][:j+1]))
            for j in range(len(ls)-1):
                if p1 > ls[j] and p1 < ls[j+1]:
                    break
            lgba.append(j+1)                 ####batch_label
        for i in range(param[0]):
            lgb1 = []
            for i in range(len(param[25])):
                k = np.random.normal(param[26][0],param[26][1],1)
                lgb1.append(k)
            lgb2 = list(map(float,lgb1))
            ls = []
            for j in range(param[2]):
                ls.append(np.exp(lgb2[lgba[j]-1]))
            lgb3.append(ls)
        lba1 = []
        lba1 = lgb3
            
    if len(param[27]) == 1:
        lgb.extend(lgb3)
    else:
        a = []
        global lgbb,lgbc
        lgbb = []
        lgbc = []
        for i in range(len(param[27])):
            a.append(i+1)
        lgbb = list(np.random.choice(a, param[2], p=param[27]))  ####group
        lgb4 = []
        b = [0,1]
        lgb4 = list(np.random.choice(b, param[0], p=[1-param[28],param[28]]))
        global lgb5,l_marker
        l_marker = []
        lgb5 = []
        for i in range(len(lgb4)):
            if lgb4[i] == 1 and lmean[i] != 0:
                lgb5.append(i)                              ###DE
        global lgbd
        lgbd = []
        lgb6 = []
        for i in range(len(lgb5)):
            while 1:
                ls = []
                for j in range(len(param[27])):
                    k = np.random.normal(param[29][0],param[29][1],1)
                    ls.append(k)
                if min(ls) < param[29][0] and max(ls) > param[29][0] and (max(ls)-min(ls)>0.5 or max(ls)-min(ls)>1.5*param[29][1]):
                    break
            p1 = random.random()
            if p1 < param[-2]/param[28]:
                l_marker.append(lgb5[i])
                l = max(ls)
                l1 = ls
                for i1 in range(len(l1)):
                    if l1[i1] != l:
                        ls[i1] == -9
            else:
                lgb6.append(ls)
        k = 0
        for i in range(param[0]):
            if i in lgb5:
                ls = []
                ls1 = []
                for j in range(param[2]):
                    ls.append(lgb3[i][j]*np.exp(lgb6[k][lgbb[j]-1]))
                    ls1.append(float(np.exp(lgb6[k][lgbb[j]-1])))
                lgbd.append(ls1)                                         ####DE
                k += 1
                lgb.append(ls)
            else:
                lgb.append(lgb3[i])
        for i in range(len(lgb)):
            if i in lgb5:
                k = 0
                sum1 = 0
                for j in range(len(lgb[i])):
                    k += 1
                    sum1 += lgb[i][j]*lmean[i]
                lgbc.append(float(sum1/k))                                     ###mean group
            else:
                lgbc.append(lmean[i])




    lpa = []
    
    if param[30] == 0:
        lpa.extend(lgb)
    else:
        global lgp,lsp,lde,lnorl,ldea,lnorla
        lgp = []
        lsp = []
        lde = []
        lnorl = []
        ldea = []
        lnorla = []
        for i in range(param[2]):
            k = random.randint(1,param[30])
            lgp.append(k)                        #####group
            k1 = random.random()
            lsp.append(k1)                       #####step
        
        b = [0,1]
        lde1 = []
        lde1 = list(np.random.choice(b, param[0], p=[1-param[28],param[28]]))
        for i in range(param[0]):
            if lde1[i] == 1 and lmean[i] != 0:
                lde.append(i)                   #####DEgene
        lnorl1 = []
        lnorl1 = list(np.random.choice(b, param[0], p=[1-param[31]/(1-param[28]),param[31]/(1-param[28])]))
        for i in range(param[0]):
            if lnorl1[i] == 1 and lmean[i] != 0:
                if i not in lde:
                    lnorl.append(i)              ######non-linear
        for i in range(param[0]):
            ls = []
            if i not in lde and i not in lnorl:
                for j in range(param[2]):
                    ls.append(1)
                lpa.append(ls)
            else:
                ls = []
                if i in lde:
                    while 1:
                        ls1 = []
                        for k in range(param[30]):
                            k1 = np.random.normal(param[29][0],param[29][1],1)
                            ls1.append(k1)
                        if min(ls) < param[29][0] and max(ls) > param[29][0] and (max(ls)-min(ls)>0.5 or max(ls)-min(ls)>1.5*param[29][1]):
                            break
                    for j in range(param[2]):
                        ls.append(1+lsp[j]*(np.exp(ls1[lgp[j]-1])-1))
                    lpa.append(ls)
                    ldea.append(ls)
                else:
                    ls1 = []
                    for j in range(param[30]):
                        l = []
                        t = []
                        w = []
                        t.append(0)
                        w.append(0)
                        for i in range(50):
                            t.append((i+1)*0.02)
                            k = np.random.normal(0,0.1)
                            w.append(w[i]+k)
                        b = []
                        for i in range(len(t)):
                            b.append(w[i]+1-t[i]*(w[-1]))

                        def funa(x,k1,k2,k3):
                            return k1*(x-x**2)+k2*(x**3-x**4)+k3*(x**5-x**6)+1

                        x = np.array(t)
                        y = np.array(b)
                        popt, pcov = curve_fit(funa,x,y)     
                        k1 = popt[0]
                        k2 = popt[1]
                        k3 = popt[2]
                        l.append(k1)
                        l.append(k2)
                        l.append(k3)
                        ls1.append(l)

                    for j in range(param[2]):
                        if ls1[(lgp[j]-1)][0]*(lsp[j]-lsp[j]**2)+ls1[(lgp[j]-1)][1]*(lsp[j]**3-lsp[j]**4)+ls1[(lgp[j]-1)][2]*(lsp[j]**5-lsp[j]**6)+1 >= 0:
                            ls.append(ls1[(lgp[j]-1)][0]*(lsp[j]-lsp[j]**2)+ls1[(lgp[j]-1)][1]*(lsp[j]**3-lsp[j]**4)+ls1[(lgp[j]-1)][2]*(lsp[j]**5-lsp[j]**6)+1)
                        else:
                            ls.append(0)
                    lpa.append(ls)
                    lnorla.append(ls)
                    

    if param[25] == [1] and param[27][0] == 1 and param[30] == 0:
        for i in range(len(lmean)):
            ls1 = []
            ls11 = []
            if lmean[i] == 0:
                for j in range(param[2]):
                    ls11.append(0)
                lgenebcv2.append(0)
            else:
                if lsa[i] == 0:
                    for i in range(param[2]):
                        ls11.append(lmean[i])
                    lgenebcv2.append(0)
                else:
                    ls11 = np.random.gamma(1/lsa[i],lmean[i]*lsa[i],param[2])
                    lgenebcv2.append(float(lsa[i]))
            ls1 = list(map(float,ls11))
            lbcv.append(ls1)
    else:
        if len(param[25]) > 1:
            for i in range(len(lmean)):
                ls1 = []
                ls11 = []
                if lmean[i] == 0:
                    for j in range(param[2]):
                        ls11.append(0)
                    lgenebcv2.append(0)
                else:
                    if lsa[i] == 0:
                        for j in range(param[2]):
                            ls11.append(lmean[i])
                    else:
                        for j in range(param[2]):
                            ls11.append(np.random.gamma(1/lsa[i],lpa[i][j]*lmean[i]*lsa[i],1))
                    lgenebcv2.append(float(lsa[i]))
                ls1 = list(map(float,ls11))
                lbcv.append(ls1)

            
        else:
            for i in range(len(lmean)):
                ls1 = []
                ls11 = []
                if lmean[i] == 0:
                    for j in range(param[2]):
                        ls11.append(0)
                    lgenebcv2.append(0)
                else:
                    if i not in lgb5:
                        if lsa[i] == 0:
                            for i in range(param[2]):
                                ls11.append(lmean[i])
                        else:
                            ls11 = np.random.gamma(1/lsa[i],lmean[i]*lsa[i],param[2])
                    else:
                        if lsa[i] == 0:
                            for i in range(param[2]):
                                ls11.append(lpa[i][j]*lmean[i])
                        else:
                            for i in range(param[2]):
                                ls11.append(np.random.gamma(1/lsa[i],lpa[i][j]*lmean[i]*lsa[i],1))
                    lgenebcv2.append(float(lsa[i]))
                ls1 = list(map(float,ls11))
                lbcv.append(ls1)
            
    
    
    
    

def creatsizef():
    """
    Simulate the size factor

    Parameters
    ----------
    param:list
          All parameters used in simulation
    
    Returns
    -------
    lsizef:list
          Simulated size factor gene
    lsf:list
          Simulated size factor of each cell
          
    """
    global lsizef,lsf
    lsizef = []
    ls = []
    lsf = []
    lsf1 = []
    while 1:
        k = float(np.random.normal(param[3][3],param[3][4],1))
        if k > param[4][0] and k < param[4][1]:
            lsf1.append(k)
        if len(lsf1) == int(param[2]*param[3][0]):
            break
    lsf2 = []
    while 1:
        k = float(np.random.normal(param[3][5],param[3][6],1))
        if k > param[4][0] and k < param[4][1]:
            lsf2.append(k)
        if len(lsf2) == int(param[2]*param[3][1]):
            break
    lsf3 = []
    while 1:
        k = float(np.random.normal(param[3][7],param[3][8],1))
        if k > param[4][0] and k < param[4][1]:
            lsf3.append(k)
        if len(lsf3) == param[2]-int(param[2]*param[3][0])-int(param[2]*param[3][1]):
            break
    lsf.extend(lsf1)
    lsf.extend(lsf2)
    lsf.extend(lsf3)

    shuffle(lsf)
    
    ls1 = []
    for i in lsf:
        ls1.append(float(np.exp(i)))
    if param[-1] == 1:
        for i in range(len(lbcv)):
            ls2 = []
            for j in range(len(lbcv[i])):
                ls2.append(lbcv[i][j]*ls1[j])
            lsizef.append(ls2)
    else:
        for i in range(len(lbcv)):
            ls2 = []
            for j in range(len(lbcv[i])):
                ls2.append(lbcv[i][j]*ls1[j]*param[-1])
            lsizef.append(ls2)

def creatpoi():
    """
    Simulate poisson sampling

    Parameters
    ----------
    param:list
          All parameters used in simulation
    
    Returns
    -------
    lpoi:list
          Simulated poisson sampling gene
          
    """
    global lpoi,lmean
    lpoi = []
    i = 0
    kkk = 0
    while 1:
        ls = []
        ls1 = []
        if lmean[i] == 0:
            for  j in range(len(lsizef[i])):
                ls1.append(0)
            lpoi.append(ls1)
            i += 1
        else:
            kkk += 1
            x = np.random.poisson(tuple(lsizef[i]),(1,param[2]))
            ls1 = list(map(int,x[0]))
            if max(ls1) > 0 or kkk >= 20:
                i += 1
                lpoi.append(ls1)
                kkk = 0
        if len(lpoi) == param[0]:
            break
    if mod == '-a':
        global llisize
        llisize = []
        for i in range(len(lpoi[0])):
            k = 0
            for j in range(len(lpoi)):
                k += lpoi[j][i]
            llisize.append(k)
        

def creatdropout():
    """
    Simulate dropout

    Parameters
    ----------
    param:list
          All parameters used in simulation
    
    Returns
    -------
    ldp:list
          Simulated dropouted gene
    ldropout:list
          Dropout information
    llisize:list
          Library size
          
    """
    global ldp,ldropout,llisize
    ldropout = []
    ldp = []
    llisize = []
    def randomm(m): 
        k = random.random()
        if k < m:
            return 0
        else:
            return 1

    kkk = 0
    i = 0
    while 1:
        ls1 = []
        ls2 = []
        if lmean[i] == 0:
            for j in range(len(lpoi[i])):
                ls1.append(0)
                ls2.append('uncertain')
            ldp.append(ls1)
            ldropout.append(ls2)
            i += 1
        else:
            kkk += 1
            for j in range(len(lpoi[i])):
                if lsizef[i][j] == 0:
                    ls1.append(0)
                    ls2.append('true')
                else:
                    k = randomm(param[5][0]*np.exp(-param[5][1]*lsizef[i][j]))
                    ls1.append(lpoi[i][j]*k)
                    if k == 0:
                        ls2.append('true')
                    else:
                        ls2.append('false')
            if max(ls1) > 0 or kkk >=20:
                i += 1
                ldp.append(ls1)
                ldropout.append(ls2)
                kkk = 0
            else:
                continue

        if len(ldp) == param[0]:
            break

    for i in range(len(ldp[0])):
        k = 0
        for j in range(len(ldp)):
            k += ldp[j][i]
        llisize.append(k)

    
def write():
    """
    Write file

    Parameters
    ----------
    lmean:list
          Simulated gene mean
    lbcv:list
          Simulated gene BCV
    lgenebcv2:list
          BCV² of each gene
    lsf:list
          Simulated size factor of each cell
    lpoi:list
          Simulated poisson sampling gene
    ldp:list
          Simulated dropouted gene
    ldropout:list
          Dropout information
    llisize:list
          Library size
    
    Returns
    -------
    All simulated file

    """
    global lmean,lpoi,ldp,lsizef,lbcv,ldropout,lgenebcv2,llisize,lsf
    f = open(path2+'/Simulation information.txt','w')
    if mod == '-a':
        f.write('Method:SimCH-flex-NB')
        f.write('\n')
        f.write('Data:'+path1)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('-1     '+'Number of genes:{}'.format(param[0]))
        f.write('\n')
        f.write('-2     '+'Number of non-zero genes:{}'.format(param[1]))
        f.write('\n')
        f.write('-3     '+'Number of cells:{}'.format(param[2]))
        f.write('\n')
        f.write('-4     '+'Ln(Sj) GMM3:{}'.format(param[3]))
        f.write('\n')
        f.write('-5     '+'Ln(mean) GMM5:{}'.format(param[6][2:]))
        f.write('\n')
        f.write('-6     '+'Relationship parameter of mean and dispersion:{} '.format(param[19:23]))
        f.write('\n')
        f.write('-7     '+'Degree of freedom:{}'.format(param[23]))
        f.write('\n')
        f.write('-8     '+'Group ratio:{}'.format(param[27]))
        f.write('\n')
        f.write('-9     '+'Batch ratio:{}'.format(param[25]))
        f.write('\n')
        f.write('-10    '+'Batch variation:{}'.format(param[26]))
        f.write('\n')
        f.write('-11    '+'Number of path:{}'.format(param[30]))
        f.write('\n')
        f.write('-12    '+'DEgene ratio:{}'.format(param[28]))
        f.write('\n')
        f.write('-13    '+'DEgene variation:{}'.format(param[29]))
        f.write('\n')
        f.write('-14    '+'Non-linear gene ratio:{}'.format(param[31]))
        f.write('\n')
        f.write('-15    '+'Marker gene ratio:{}'.format(param[32]))
        f.write('\n')
        f.write('-16    '+'Library magnification:{}'.format(param[33]))
        f.write('\n')
        f.close()
    else:
        f.write('Method:SimCH-flex-NBZI')
        f.write('\n')
        f.write('Data:'+path1)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('-1     '+'Number of genes:{}'.format(param[0]))
        f.write('\n')
        f.write('-2     '+'Number of non-zero genes:{}'.format(param[1]))
        f.write('\n')
        f.write('-3     '+'Number of cells:{}'.format(param[2]))
        f.write('\n')
        f.write('-4     '+'Ln(Sj) GMM3:{}'.format(param[3]))
        f.write('\n')
        f.write('-5     '+'γ and ε:{} '.format(param[5]))
        f.write('\n')
        f.write('-6     '+'Ln(mean) GMM5:{}'.format(param[6][2:]))
        f.write('\n')
        f.write('-7     '+'Relationship parameter of mean and dispersion:{} '.format(param[19:23]))
        f.write('\n')
        f.write('-8     '+'Degree of freedom:{}'.format(param[23]))
        f.write('\n')
        f.write('-9     '+'Group ratio:{}'.format(param[27]))
        f.write('\n')
        f.write('-10     '+'Batch ratio:{}'.format(param[25]))
        f.write('\n')
        f.write('-11    '+'Batch variation:{}'.format(param[26]))
        f.write('\n')
        f.write('-12    '+'Number of path:{}'.format(param[30]))
        f.write('\n')
        f.write('-13    '+'DEgene ratio:{}'.format(param[28]))
        f.write('\n')
        f.write('-14    '+'DEgene variation:{}'.format(param[29]))
        f.write('\n')
        f.write('-15    '+'Non-linear gene ratio:{}'.format(param[31]))
        f.write('\n')
        f.write('-16    '+'Marker gene ratio:{}'.format(param[32]))
        f.write('\n')
        f.write('-17    '+'Library magnification:{}'.format(param[33]))
        f.write('\n')
        f.close()




    if mod == '-b':
        f = open(path2+'/TrueCounts.txt','w')
    else:
        f = open(path2+'/Counts.txt','w')
    for i in range(len(lpoi[0])):
        f.write('\t')
        f.write('cell'+str(i+1))
    f.write('\n')
    for i in range(len(lpoi)):
        f.write('gene'+str(i+1))
        f.write('\t')
        for j in range(len(lpoi[i])):
            f.write(str(lpoi[i][j]))
            if j < len(lpoi[i])-1:
                f.write('\t')
        f.write('\n')
    f.close()
    if mod == '-b':
        f1 = open(path2+'/Counts.txt','w')
        for i in range(len(lpoi[0])):
            f1.write('\t')
            f1.write('cell'+str(i+1))
        f1.write('\n')
        for i in range(len(ldp)):
            f1.write('gene'+str(i+1))
            f1.write('\t')
            for j in range(len(ldp[i])):
                f1.write(str(ldp[i][j]))
                if j < len(ldp[i])-1:
                    f1.write('\t')
            f1.write('\n')
        f1.close()

    f2 = open(path2+'/Cell gene mean(size factor).txt','w')
    for i in range(len(lsizef[0])):
        f2.write('\t')
        f2.write('cell'+str(i+1))
    f2.write('\n')
    for i in range(len(lsizef)):
        f2.write('gene'+str(i+1))
        f2.write('\t')
        for j in range(len(lsizef[i])):
            f2.write(str(lsizef[i][j]))
            if j < len(lsizef[i])-1:
                f2.write('\t')
        f2.write('\n')
    f2.close()

    f2 = open(path2+'/Cell gene mean(bcv2).txt','w')
    for i in range(len(lbcv[0])):
        f2.write('\t')
        f2.write('cell'+str(i+1))
    f2.write('\n')
    for i in range(len(lbcv)):
        f2.write('gene'+str(i+1))
        f2.write('\t')
        for j in range(len(lbcv[i])):
            f2.write(str(lbcv[i][j]))
            if j < len(lbcv[i])-1:
                f2.write('\t')
        f2.write('\n')
    f2.close()

    if mod == '-b':
        f3 = open(path2+'/Dropout.txt','w')
        for i in range(len(ldropout[0])):
            f3.write('\t')
            f3.write('cell'+str(i+1))
        f3.write('\n')
        for i in range(len(ldropout)):
            f3.write('gene'+str(i+1))
            f3.write('\t')
            for j in range(len(ldropout[i])):
                f3.write(str(ldropout[i][j]))
                if j < len(ldropout[i])-1:
                    f3.write('\t')
            f3.write('\n')
        f3.close()

    f4 = open(path2+'/Gene mean.txt','w')
    for i in range(len(lmean)):
        f4.write('gene'+str(i+1))
        f4.write('\t')
        f4.write(str(lmean[i]))
        f4.write('\n')
    f4.close()

    f5 = open(path2+'/Gene dispersion.txt','w')
    for i in range(len(lgenebcv2)):
        f5.write('gene'+str(i+1))
        f5.write('\t')
        f5.write(str(lgenebcv2[i]))
        f5.write('\n')
    f5.close()
        
    f6 = open(path2+'/Library size.txt','w')
    for i in range(len(llisize)):
        if i < len(llisize)-1:
            f6.write('cell'+str(i+1))
            f6.write('\t')
        else:
            f6.write('cell'+str(i+1))
            f6.write('\n')
    for i in range(len(llisize)):
        if i < len(llisize)-1:
            f6.write(str(llisize[i]))
            f6.write('\t')
        else:
            f6.write(str(llisize[i]))
            f6.write('\n')
    f6.close()

    f7 = open(path2+'/Size factor.txt','w')
    for i in range(len(lsf)):
        if i < len(lsf)-1:
            f7.write('cell'+str(i+1))
            f7.write('\t')
        else:
            f7.write('cell'+str(i+1))
            f7.write('\n')
    for i in range(len(lsf)):
        if i < len(lsf)-1:
            f7.write(str(np.exp(lsf[i])))
            f7.write('\t')
        else:
            f7.write(str(np.exp(lsf[i])))
            f7.write('\n')
    f7.close()
    if param[25] != [1]:
        global lgba,lba1
        f8 = open(path2+'/Batch.txt','w')
        for i in range(len(lgba)):
            if i < len(lgba)-1:
                f8.write('cell'+str(i+1))
                f8.write('\t')
            else:
                f8.write('cell'+str(i+1))
                f8.write('\n')
        for i in range(len(lgba)):
            if i < len(lgba)-1:
                f8.write(str(lgba[i]))
                f8.write('\t')
            else:
                f8.write(str(lgba[i]))
                f8.write('\n')
        f8.close()


        f8 = open(path2+'/Batch factor.txt','w')

        for i in range(len(lgba)):
            f8.write('\t')
            f8.write('cell'+str(i+1))
        f8.write('\n')
        for i in range(len(lba1)):
            f8.write('gene'+str(i+1))
            for j in range(len(lba1[i])):
                f8.write('\t')
                f8.write(str(lba1[i][j]))
            f8.write('\n')
        f8.close()

        lgba = []
        lba1 = []
        
    if len(param[27]) != 1:
        global lgbb,lgb5,lgbd,lgbc
        f8 = open(path2+'/Group.txt','w')
        for i in range(len(lgbb)):
            if i < len(lgbb)-1:
                f8.write('cell'+str(i+1))
                f8.write('\t')
            else:
                f8.write('cell'+str(i+1))
                f8.write('\n')
        for i in range(len(lgbb)):
            if i < len(lgbb)-1:
                f8.write(str(lgbb[i]))
                f8.write('\t')
            else:
                f8.write(str(lgbb[i]))
                f8.write('\n')
        f8.close()


        f8 = open(path2+'/Marker gene.txt','w')
        for i in range(param[0]):
            if i in l_marker:
                f8.write('gene'+str(i+1))
                f.write('\n')
        f8.close()


        f8 = open(path2+'/DEgene factor.txt','w')
        for i in range(param[2]):
            f8.write('\t')
            f8.write('cell'+str(i+1))
        f8.write('\n')
        for i in range(len(lgb5)):
            f8.write('gene'+str(lgb5[i]+1))
            for j in range(param[2]):
                f8.write('\t')
                f8.write(str(lgbd[i][j]))
            f8.write('\n')
        f8.close()


        f8 = open(path2+'/Gene mean(group).txt','w')
        for i in range(len(lgbc)):
            f8.write('gene'+str(i+1)+'\t')
            f8.write(str(lgbc[i])+'\n')
        f8.close()
        lgbc = []
        lgbb = []
        lgb5 = []
        lgbd = []


    if param[30] != 0:
        global lgp,lsp,lde,ldea
        f8 = open(path2+'/Path.txt','w')
        for i in range(param[2]):
            if i < len(lgp)-1:
                f8.write('cell'+str(i+1))
                f8.write('\t')
            else:
                f8.write('cell'+str(i+1))
                f8.write('\n')
        for i in range(len(lgp)):
            if i < len(lgp)-1:
                f8.write(str(lgp[i]))
                f8.write('\t')
            else:
                f8.write(str(lgp[i]))
                f8.write('\n')
        for i in range(len(lsp)):
            if i < len(lsp)-1:
                f8.write(str(lsp[i]))
                f8.write('\t')
            else:
                f8.write(str(lsp[i]))
                f8.write('\n')
        f8.close()

        
        f8 = open(path2+'/DEgene factor.txt','w')
        for i in range(param[2]):
            f8.write('\t')
            f8.write('cell'+str(i+1))
        f8.write('\n')
        for i in range(len(lde)):
            f8.write('gene'+str(lde[i]+1))
            for j in range(param[2]):
                f8.write('\t')
                f8.write(str(list(ldea[i][j])[0]))
            f8.write('\n')
        f8.close()
        lgp = []
        lsp = []
        lde = []
        ldea = []

        global lnorl,lnorla
        if len(lnorl) > 0:
            f8 = open(path2+'/noline gene factor','w')
            for i in range(param[2]):
                f8.write('\t')
                f8.write('cell'+str(i+1))
            f8.write('\n')
            for i in range(len(lnorl)):
                f8.write('gene'+str(lnorl[i]+1))
                for j in range(param[2]):
                    f8.write('\t')
                    f8.write(str(lnorla[i][j]))
                f8.write('\n')
            f8.close()
            lnorl = []
            lnorla = []
    lmean = []
    lbcv = []
    lgenebcv2 = []
    lsizef = []
    lsf = []
    lpoi = []
    ldp = []
    ldropout = []
    llisize = []
    

def Remove0_noparam():
    """
    Read the file and remove the 0 gene expressed

    Parameters
    ----------
    path1:str
          File path
    param:list
          All parameters used in simulation
          
    Returns
    -------
    ls2:list
          The data of gene

    """
    f = open(path1,'r')
    if path1[-3:] == 'csv':
        reader = csv.reader(f)
        ls = list(reader)
    else:

        f1 = f.readlines()

        ls = []
        for i in f1:
            k = i.strip().split('\t')
            ls.append(k)

    if len(ls[0]) < len(ls[1]):
        ls[0].insert(0,' ')

    ls1 = []
    for i in range(len(ls)-1):
        if float(max(ls[i+1][1:])) == 0:
            ls1.append(i+1)

    print('Number of genes:{}'.format(len(ls)-1))
    param.append(len(ls)-1)
    global ls2,ls2_0
    ls2 = []
    ls2_0 = []
    for i in range(len(ls)):
        if not i in ls1:
            ls2.append(ls[i])
    ls2_0 = [i-1 for i in ls1]
    f.close()

    
def nml_noparam():
    """
    Normalize data with size factor
    If there are too few housekeeping genes,we will reduce the selection conditions about it

    Parameters
    ----------
    ls2:list
          The data of gene
    
          
    Returns
    -------
    ls3:list
          The data of nomalized gene
    ls7:list
          Size factor

    """
    global ls3
    ls3 = []
    global ls7                              #size factor

    ls7 = []
    ls = []
    for i in range(len(ls2[0])-1):
        sum1 = 0
        for j in range(len(ls2)-1):
            sum1 += float(ls2[j+1][i+1])
        ls.append(sum1)
    ls1 = sorted(ls)
    if len(ls1)%2 == 0:
        k = (ls1[len(ls1)//2-1]+ls1[len(ls1)//2])/2
    else:
        k = ls1[len(ls1)//2]
    for i in ls:
        ls7.append(i/k)
        

    ls3.append(ls2[0])
    for i in range(len(ls2)-1):
        ls8 = []
        ls8.append(ls2[i+1][0])
        for j in range(len(ls2[i+1])-1):
            ls8.append(float(ls2[i+1][j+1])/ls7[j])
        ls3.append(ls8)
    print('Number of non-zero genes:{}'.format(len(ls3)-1))
    param.append(len(ls3)-1)
    print('Number of cells:{}'.format(len(ls3[1])-1))
    param.append(len(ls3[1])-1)


    


def NBzero_mean_dispersion_noparam():
    """
    Estimate mean and dispersion of -d or -d nocopula

    Parameters
    ----------
    ls3:list
          The data of nomalized gene
     
    Returns
    -------
    ls4:list
          Mean of gene
    ls12c:list
          Dispersion of gene

    """
    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')


    global ls4,ls4p,ls2,ls7,l_c,dictx,ls12c,k,k_a,correlation,lsk_a,lsk_b,lp0


    lsk_a = []
    lsk_b = []
    ls4 = []
    ls4p = []
    lp0 = []
    for i in range(len(ls3)-1):
        ls4p.append(ls3[i+1][1:])
    for i in range(len(ls4p)):
        ls4.append(float(np.mean(ls4p[i])))
    ls = []
    for i in range(len(ls4)):
        ls.append((i,ls4[i]))
    ls.sort(key=lambda x:x[1])


    ls1 = []
    ls1a = []
    l = []
    for i in range(len(ls)):

        if ls4p[ls[-(i+1)][0]].count(0)/len(ls4p[0]) > 0:
            ls1.append(ls[-(i+1)][0])
            ls1a.append(len(l))
            l = []
            l.append(1)
        else:
            l.append(1)
        if len(ls1) > len(ls4p)*0.02:
            ls1a.append(len(l))
            break



    ls = []
    ls_2 = []
    for i in range(len(ls1)):
        ls.append(ls4p[ls1[i]])
        ls_2.append(ls1a[i+1])
    ls1 = []
    for i in ls:
        ls1.append(i)
    ls = []
    ls_1 = []
    for i in range(len(ls1)):
        p = random.random()
        if p < 0.1:
            ls.append(ls1[i])
            ls_1.append(ls_2[i])


    lsd = []
    lsm = []
    lsd2 = []
    for i in range(len(ls)):
        lsd.append(float(np.var(ls[i])/(np.mean(ls[i])**2)-1/np.mean(ls[i])))
        lsm.append(np.mean(ls[i]))
        lsd2.append(ls[i])



    x = Symbol('x')
    global lme,ldis
    lme = []
    ldis = []
    k1 = 0
    s = np.array(ls7)
    print('Estimate mean,dispersion,γ,ε...')
    k3 = len(ls4p)//100
    kk = 0
    warnings.filterwarnings("ignore")
    for i in range(len(ls4p)):
        if i%k3 == 0:
            progress(int(i/k3),width=50)
            #print('{}%'.format(int(i/k1)))
        m0 = float(np.mean(ls4p[i]))
        p0 = ls4p[i].count(0)/len(ls4p[i])
        lp0.append(p0)
        v0 = float(np.var(ls4p[i]))
        d0 = v0/(m0**2)
        

        ##########################################################################
            
            
        l1 = []
        l2 = []
        l3 = []
        l4 = []
        step1 = (1/(1-p0))**(0.2)
        for i1 in range(7):
            for j1 in range(5):
                for k1 in range(5):
                    for k2 in range(10):
                        if i1 == 0:
                            l1.append([m0])
                            l2.append([0.6*d0+j1*0.1*d0])
                            l3.append([0.1+0.2*k1])
                            l4.append([-5+k2*0.56])
                        elif i1 == 1:
                            if (step1-1)*m0*0.1+m0 > 0:
                                l1.append([(step1-1)*m0*0.1+m0])
                                l2.append([0.6*d0+j1*0.1*d0])
                                l3.append([0.1+0.2*k1])
                                l4.append([-5+k2*0.56])
                        else:
                            l1.append([m0*(step1**(i1-1))])
                            l2.append([0.6*d0+j1*0.1*d0])
                            l3.append([0.1+0.2*k1])
                            l4.append([-5+k2*0.56])
        l1 = np.array(l1)
        l2 = np.array(l2)
        l3 = np.array(l3)
        l4 = np.array(np.exp(l4))
        d = l2
        m = l1
        p1 = l3
        k = l4
        g = np.exp(gammaln(1/d+1)-gammaln(1/d))
        s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
        if p0 > 0:
            s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
        g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
        m1 = np.transpose(np.array(s2.mean(axis=1)))
        ls = []
        for i1 in m1:
            ls.append([i1])
        m1  = np.array(ls)

        s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
        s4 = s3.mean(axis=1)
        d1 = []
        for i1 in s4:
            d1.append([i1])
        d1 = np.array(d1)
        d1 = d1-m1**2
        if p0 > 0:
            s5 = s1.mean(axis=1)
            ls = []
            for i1 in s5:
                ls.append([i1])
            p1 = np.array(ls)

        if p0 > 0:

            if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            else:
                l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

        else:
            l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
        r = np.argmin(l)
        m = l1[r][0]
        d = l2[r][0]
        p1 = l3[r][0]
        k = np.log(l4[r][0])



        l1 = []
        l2 = []
        l3 = []
        l4 = []
        step1 = (1/(1-p0))**(0.2)
        for i1 in range(5):
            for j1 in range(5):
                for k1 in range(5):
                    for k2 in range(5):
                        if p1+(k1-2)*0.1 >= 0 and p1+(k1-2)*0.1 <=1:
                            if m == m0:
                                if (step1-1)*m0*0.1*i1*0.25+m0 > 0:
                                    l1.append([(step1-1)*m0*0.1*i1*0.25+m0])
                                    l2.append([d+(j1-2)*0.05*d0])
                                    l3.append([p1+(k1-2)*0.1])
                                    l4.append([k+(k2-2)*0.28])
                            elif m == (step1-1)*m0*0.1+m0:
                                if (step1-1)*m0*0.1*i1*0.5+m0 > 0:
                                    l1.append([(step1-1)*m0*0.1*i1*0.5+m0])
                                    l2.append([d+(j1-2)*0.05*d0])
                                    l3.append([p1+(k1-2)*0.1])
                                    l4.append([k+(k2-2)*0.28])
                            else:
                                l1.append([m*(step1**((i1-2)*0.5))])
                                l2.append([d+(j1-2)*0.05*d0])
                                l3.append([p1+(k1-2)*0.1])
                                l4.append([k+(k2-2)*0.28])
        l1 = np.array(l1)
        l2 = np.array(l2)
        l3 = np.array(l3)
        l4 = np.array(np.exp(l4))
        d = l2
        m = l1
        p1 = l3
        k = l4
        g = np.exp(gammaln(1/d+1)-gammaln(1/d))
        s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
        if p0 > 0:
            s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
        g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
        m1 = np.transpose(np.array(s2.mean(axis=1)))
        ls = []
        for i1 in m1:
            ls.append([i1])
        m1  = np.array(ls)

        s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
        s4 = s3.mean(axis=1)
        d1 = []
        for i1 in s4:
            d1.append([i1])
        d1 = np.array(d1)
        d1 = d1-m1**2
        if p0 > 0:
            s5 = s1.mean(axis=1)
            ls = []
            for i1 in s5:
                ls.append([i1])
            p1 = np.array(ls)

        if p0 > 0:

            if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            else:
                l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

        else:
            l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
        r = np.argmin(l)
        m = l1[r][0]
        d = l2[r][0]
        p1 = l3[r][0]
        k = np.log(l4[r][0])





        l1 = []
        l2 = []
        l3 = []
        l4 = []
        step1 = (1/(1-p0))**0.2
        for i1 in range(5):
            for j1 in range(5):
                for k1 in range(5):
                    for k2 in range(5):
                        if p1+(k1-2)*0.1 >= 0 and p1+(k1-2)*0.1 <=1:
                            if m == m0:
                                if (step1-1)*m0*0.1*i1*0.125+m0 > 0:
                                    l1.append([(step1-1)*m0*0.1*i1*0.125+m0])
                                    l2.append([d+(j1-2)*0.025*d0])
                                    l3.append([p1+(k1-2)*0.05])
                                    l4.append([k+(k2-2)*0.14])
                            elif m <= (step1-1)*m0*0.2+m0:
                                if (step1-1)*m0*0.05*(i1-2)*0.25+m > 0:
                                    l1.append([(step1-1)*m0*0.05*(i1-2)*0.25+m])
                                    l2.append([d+(j1-2)*0.025*d0])
                                    l3.append([p1+(k1-2)*0.05])
                                    l4.append([k+(k2-2)*0.14])
                            else:
                                l1.append([m*(step1**((i1-2)*0.25))])
                                l2.append([d+(j1-2)*0.025*d0])
                                l3.append([p1+(k1-2)*0.05])
                                l4.append([k+(k2-2)*0.14])
        l1 = np.array(l1)
        l2 = np.array(l2)
        l3 = np.array(l3)
        l4 = np.array(np.exp(l4))
        d = l2
        m = l1
        p1 = l3
        k = l4
        g = np.exp(gammaln(1/d+1)-gammaln(1/d))
        s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
        if p0 > 0:
            s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
        g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
        m1 = np.transpose(np.array(s2.mean(axis=1)))
        ls = []
        for i1 in m1:
            ls.append([i1])
        m1  = np.array(ls)

        s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
        s4 = s3.mean(axis=1)
        d1 = []
        for i1 in s4:
            d1.append([i1])
        d1 = np.array(d1)
        d1 = d1-m1**2
        if p0 > 0:
            s5 = s1.mean(axis=1)
            ls = []
            for i1 in s5:
                ls.append([i1])
            p1 = np.array(ls)

        if p0 > 0:

            if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            else:
                l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

        else:
            l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
        r = np.argmin(l)
        m = l1[r][0]
        d = l2[r][0]
        p1 = l3[r][0]
        k = np.log(l4[r][0])



        l1 = []
        l2 = []
        l3 = []
        l4 = []
        step1 = (1/(1-p0))**0.2
        for i1 in range(5):
            for j1 in range(5):
                for k1 in range(5):
                    for k2 in range(5):
                        if p1+(k1-2)*0.1 >= 0 and p1+(k1-2)*0.1 <=1:
                            if m == m0:
                                if (step1-1)*m0*0.1*i1*0.0625+m0 > 0:
                                    l1.append([(step1-1)*m0*0.1*i1*0.0625+m0])
                                    l2.append([d+(j1-2)*0.0125*d0])
                                    l3.append([p1+(k1-2)*0.025])
                                    l4.append([k+(k2-2)*0.07])
                            elif m <= (step1-1)*m0*0.2+m0:
                                if (step1-1)*m0*0.05*(i1-2)*0.125+m > 0:
                                    l1.append([(step1-1)*m0*0.05*(i1-2)*0.125+m])
                                    l2.append([d+(j1-2)*0.0125*d0])
                                    l3.append([p1+(k1-2)*0.025])
                                    l4.append([k+(k2-2)*0.07])
                            else:
                                l1.append([m*(step1**((i1-2)*0.125))])
                                l2.append([d+(j1-2)*0.0125*d0])
                                l3.append([p1+(k1-2)*0.025])
                                l4.append([k+(k2-2)*0.07])
        l1 = np.array(l1)
        l2 = np.array(l2)
        l3 = np.array(l3)
        l4 = np.array(np.exp(l4))
        d = l2
        m = l1
        p1 = l3
        k = l4
        g = np.exp(gammaln(1/d+1)-gammaln(1/d))
        s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
        if p0 > 0:
            s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
        g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
        m1 = np.transpose(np.array(s2.mean(axis=1)))
        ls = []
        for i1 in m1:
            ls.append([i1])
        m1  = np.array(ls)

        s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
        s4 = s3.mean(axis=1)
        d1 = []
        for i1 in s4:
            d1.append([i1])
        d1 = np.array(d1)
        d1 = d1-m1**2
        if p0 > 0:
            s5 = s1.mean(axis=1)
            ls = []
            for i1 in s5:
                ls.append([i1])
            p1 = np.array(ls)

        if p0 > 0:

            if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            else:
                l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

        else:
            l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
        r = np.argmin(l)












        M = abs(m0-m1)[r]
        D = (abs(v0-d1)/v0)[r]
        if p0 > 0:
            P = (abs(p0-p1)/p0)[r]
        else:
            P = 0.05
        r_1 = r
        m1_1 = m1
        d1_1 = d1
        p1_1 = p1

        m = l1[r][0]
        d = l2[r][0]
        p1 = l3[r][0]
        k = l4[r][0]


        L = min(l)[0]
        L1 = L


        m_result = m
        d_result = d
        p1_result = p1
        k_result = k


        for i2 in range(20):

            if D < 0.05:
                if M > 0:
                    l1 = np.random.normal(0,M,64)
                    l2 = np.random.normal(0,0.025*d_result,64)
                    if P*p0 > 0:
                        l3 = np.random.normal(0,P*p0,64)
                    else:
                        l3 = np.random.normal(0,p1_result*0.001,64)
                    if (m_result*p1_result) <= 0 or P*p0 <= 0:
                        l4 = np.random.normal(0,k_result*0.1,64)
                    else:
                        if P*p0/(m_result*p1_result) > 0.001 and P*p0/(m_result*p1_result) < 10:
                            l4 = np.random.normal(0,P*p0/(m_result*p1_result),64)
                        else:
                            l4 = np.random.normal(0,k_result*0.1,64)
                else:
                    l1 = np.random.normal(0,m_result*0.001,64)
                    l2 = np.random.normal(0,0.025*d_result,64)
                    if P*p0 > 0:
                        l3 = np.random.normal(0,P*p0,64)
                    else:
                        l3 = np.random.normal(0,p1_result*0.001,64)
                    if (m_result*p1_result) <= 0 or P*p0 <= 0:
                        l4 = np.random.normal(0,k_result*0.1,64)
                    else:
                        if P*p0/(m_result*p1_result) > 0.001 and P*p0/(m_result*p1_result) < 10:
                            l4 = np.random.normal(0,P*p0/(m_result*p1_result),64)
                        else:
                            l4 = np.random.normal(0,k_result*0.1,64)
            else:
                l1 = np.random.normal(0,0.01*m_result,64)
                l2 = np.random.normal(0,D*d_result,64)
                if P*p0 > 0:
                    l3 = np.random.normal(0,P*p0,64)
                else:
                    l3 = np.random.normal(0,p1_result*0.001,64)
                if (m_result*p1_result) <= 0 or P*p0 <= 0:
                    l4 = np.random.normal(0,k_result*0.1,64)
                else:
                    if P*p0/(m_result*p1_result) > 0.001 and P*p0/(m_result*p1_result) < 10:
                        l4 = np.random.normal(0,P*p0/(m_result*p1_result),64)
                    else:
                        l4 = np.random.normal(0,k_result*0.1,64)


            l1a = []
            l2a = []
            l3a = []
            l4a = []

            for j in range(64):
                if m_result+l1[j] > 0 and d_result+l2[j] > 0 and p1_result+l3[j] > 0 and p1_result+l3[j] <= 1 and k_result+l4[j] > 0:
                    l1a.append([m_result+l1[j]])
                    l2a.append([d_result+l2[j]])
                    l3a.append([p1_result+l3[j]])
                    l4a.append([k_result+l4[j]])
            l1 = np.array(l1a)
            l2 = np.array(l2a)
            l3 = np.array(l3a)
            l4 = np.array(l4a)

            if len(l1a) > 0:

                d = l2
                m = l1
                p1 = l3
                k = l4
                g = np.exp(gammaln(1/d+1)-gammaln(1/d))


                s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))


                if p0 > 0:
                    s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
                g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
                m1 = np.transpose(np.array(s2.mean(axis=1)))
                ls = []
                for i1 in m1:
                    ls.append([i1])
                m1  = np.array(ls)

                s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
                s4 = s3.mean(axis=1)
                d1 = []
                for i1 in s4:
                    d1.append([i1])
                d1 = np.array(d1)
                d1 = d1-m1**2
                if p0 > 0:
                    s5 = s1.mean(axis=1)
                    ls = []
                    for i1 in s5:
                        ls.append([i1])
                    p1 = np.array(ls)



                if p0 > 0:

                    if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                        l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                    else:
                        l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

                else:
                    l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)



                if min(l)[0] < L:
                    r = np.argmin(l)
                    da = d1[r]
                    pa = p1[r]
                    ma = m1[r]

                    M = abs(m0-m1)[r]
                    D = (abs(v0-d1)/v0)[r]
                    if p0 > 0:
                        P = (abs(p0-p1)/p0)[r]
                    else:
                        P = 0.05
                    m_result = l1[r][0]
                    d_result = l2[r][0]
                    p1_result = l3[r][0]
                    k_result = l4[r][0]
                    L = min(l)[0]
        # print(D,M,P)
        # if L1 != L:
        #     print('Yes')
        #     if p0 > 0:
        #         print(L,m0,ma,v0,da,p0,pa)
        #         print(p1_result,k_result,m_result,d_result)
        #     else:
        #         print(L,m0,ma,v0,da)
        #         print(p1_result,k_result,m_result,d_result)
        # else:
        #     print('No')
        #     if p0 > 0:
        #         print(L1,m0,m1_1[r_1],v0,d1_1[r_1],p0,p1_1[r_1])
        #         print(p1_result,k_result,m_result,d_result)
        #     else:
        #         print(L1,m0,m1_1[r_1],v0,d1_1[r_1])
        #         print(p1_result,k_result,m_result,d_result)




        l1 = []
        l2 = []

        step1 = (1/(1-p0))**0.2
        for i1 in range(11):
            for j1 in range(11):
                if m_result == m0:
                    if (step1-1)*m0*0.1*i1*0.03125+m0 > 0 and d_result+(j1-5)*0.00625*d0 > 0:
                        l1.append([(step1-1)*m0*0.1*i1*0.03125+m0])
                        l2.append([d_result+(j1-5)*0.00625*d0])

                elif m_result <= (step1-1)*m0*0.2+m0:
                    if (step1-1)*m0*0.025*(i1-5)*0.125+m_result > 0 and d_result+(j1-5)*0.00625*d0 > 0:
                        l1.append([(step1-1)*m0*0.025*(i1-5)*0.125+m_result])
                        l2.append([d_result+(j1-5)*0.00625*d0])

                else:
                    if m_result*(step1**((i1-5)*0.0625)) > 0 and d_result+(j1-5)*0.00625*d0 > 0:
                        l1.append([m_result*(step1**((i1-5)*0.0625))])
                        l2.append([d_result+(j1-5)*0.00625*d0])
        if len(l1) > 0:
            l1 = np.array(l1)
            l2 = np.array(l2)

            d = l2
            m = l1

            g = np.exp(gammaln(1/d+1)-gammaln(1/d))
            s2 = d*m*g*(1-p1_result*(1+k_result*m*s*d)**(-1/d-1))
            if p0 > 0:
                s1 = (1+d*m*s)**(-1/d)*(1-p1_result*(k_result/(1/(d*m*s)+1)+1)**(-1/d))+(1+k_result*d*m*s)**(-1/d)*p1_result
            g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
            m1 = np.transpose(np.array(s2.mean(axis=1)))
            ls = []
            for i1 in m1:
                ls.append([i1])
            m1  = np.array(ls)

            s3 = ((d*m)**2*g1*(1-(1+k_result*d*m*s)**(-1/d-2)*p1_result)+(d*m/s)*g*(1-(1+k_result*d*m*s)**(-1/d-1)*p1_result))
            s4 = s3.mean(axis=1)
            d1 = []
            for i1 in s4:
                d1.append([i1])
            d1 = np.array(d1)
            d1 = d1-m1**2
            if p0 > 0:
                s5 = s1.mean(axis=1)
                ls = []
                for i1 in s5:
                    ls.append([i1])
                p1 = np.array(ls)

            if p0 > 0:

                if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                    l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5
                else:
                    l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5

            else:
                l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
            r = np.argmin(l)
            m_result = l1[r][0]
            d_result = l2[r][0]



        lme.append(m_result)
        ldis.append(d_result)
        lsk_a.append(k_result)
        lsk_b.append(p1_result)



    ls4 = lme
    ls12c = ldis
    
    dictx = {}

    l2 = []
    l3 = []

    l2 = list(reversed(list(np.argsort(lme))))
    
    l_c = []
    if mod != '-d nocopula':
        if len(ls4) > copula_number:
            for i in range(copula_number):
                l3.append(ls4p[l2[i]])
                l_c.append(l2[i])
            for i in range(len(l_c)):
                dictx[l_c[i]] = i
        else:
            for i in range(len(ls4)):
                l3.append(ls4p[l2[i]])
                l_c.append(l2[i])
            for i in range(len(l_c)):
                dictx[l_c[i]] = i

        Mat = np.array(l3)
        correlation = np.corrcoef(Mat)
        correlation1 = np.corrcoef(Mat)
        #sns.set()
        #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        #plt.show()
        
        Mat = np.array(l3[:50])
        correlation1 = np.corrcoef(Mat)
        #sns.set()
        #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        #plt.show()


def creatcount_NBzero_noparam():
    """
    Simulate count of -d or -d nocopula

    Parameters
    ----------
    ls4:list
          Mean of gene
    ls12c:list
          Dispersion of gene
    ls7:list
          Size factor

    Returns
    -------
    lcount:list
          Count

    """
    from scipy.stats import nbinom
    from scipy.stats import norm
    global lcount,ls4,correlation,ls4p,ls12c,ls7,dictx,ltruecount,ldp,k,lm,ls7_1
    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')

    lcount = []
    ltruecount = []
    ldp = []
    ls7_1 = []
    if cell_number == len(ls7):
        ls7_1 = ls7
    else:
        if cell_number < len(ls7):
            ls7_1 = ls7[:cell_number]
        else:
            ls7_2 = np.log(ls7)
            x_mean, x_std = norm.fit(ls7_2)
            ls7_1 = [i for i in ls7]
            while 1:
                k = float(np.exp(np.random.normal(x_mean, x_std,1)))
                ls7_1.append(k)
                if len(ls7_1) == cell_number:
                    break


    print(' ')
    print('Start simulation...')
    def randomm(m): 
        k = random.random()
        if k < m:
            return 0
        else:
            return 1
    def integral(m,d,p,k,k_a):
        s = np.array(p)
        return list(np.exp((-1/d)*np.log(m*d)+(-1/d-s)*np.log(1/(d*m)+1)+gammaln(1/d+s)-gammaln(1/d)-gammaln(1+s))*(1-np.exp(-k_a)*(1+k/(1/(d*m)+1))**(-1/d-s)))
        
    if mod != '-d nocopula':
        l5 = []
        for i in range(len(correlation)):
            l5.append(0)
        p = np.random.multivariate_normal(mean=l5,cov=correlation,size=len(ls7_1))

        k_2  = list(np.transpose(p))
        # k1 = np.corrcoef(np.array(k2))
        # sns.set()
        # ax = sns.clustermap(k1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        # plt.show()
        #
        # k1 = np.corrcoef(np.array(k2[:50]))
        # sns.set()
        # ax = sns.clustermap(k1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        # plt.show()


    global batch_factor,batch_label
    batch_factor = []
    batch_label = []
    if param_batch == [1] and param_library == 1:
        for i in range(len(ls7_1)):
            batch_factor.append(1)
    else:
        if param_library != 1 and param_batch == [1]:
            for i in range(len(ls7_1)):
                batch_factor.append(param_library)
        else:
            for i in range(len(ls7_1)):
                p1 = random.random()
                ls = [0]
                for j in range(len(param_batch)):
                    ls.append(sum(param_batch[:j+1]))
                for j in range(len(ls)-1):
                    if p1 > ls[j] and p1 < ls[j+1]:
                        break
                batch_label.append(j+1)                 ####batch_label
            for i in range(len(ls4p)):
                ls = []
                for j in range(len(param_batch)):
                    k1 = np.random.normal(param_batch_var[0],param_batch_var[1],1)
                    ls.append(k1)
                ls1 = []
                for j in range(len(ls7_1)):
                    ls1.append(param_library*float(np.exp(ls[batch_label[j]-1])))
                batch_factor.append(ls1)                             #######batch_factor
    global DEgene_factor,Group_label,DEgene_label,Path,Noline_gene,l_marker
    l_marker = []
    DEgene_label = []
    DEgene_factor = []
    Group_label = []
    Path = []
    Noline_gene = []
    if param_group == [1]:
        for i in range(len(ls4)):
            l = [1 for j in range(len(ls7_1))]
            DEgene_factor.append(l)
    else:
        for i in range(len(ls4)):
            p1 = random.random()
            if p1 < param_DEgene*(len(ls4)+len(ls2_0))/len(ls4):
                DEgene_label.append(1)
            else:
                DEgene_label.append(0)              ####DEgene_label

        param_group1 = []
        tree = 0
        for i in range(len(param_group)):
            if type(param_group[i]) == type([]):
                param_group1.append(sum(param_group[i]))
                tree = 1
            else:
                param_group1.append(param_group[i])


        for i in range(len(ls7_1)):
            p1 = random.random()
            ls = [0]
            for j in range(len(param_group1)):
                ls.append(sum(param_group1[:j+1]))
            for j in range(len(ls)-1):
                if p1 > ls[j] and p1 < ls[j+1]:
                    break
            Group_label.append(j+1)                 ####Group_label
        if param_path == 'No':
            for i in range(len(ls4)):
                if DEgene_label[i] == 0:
                    l = [1 for i1 in range(len(ls7_1))]
                else:
                    while 1:
                        ls = []
                        for i1 in range(len(param_group1)):
                            k1 = np.random.normal(param_DEgene_var[0],param_DEgene_var[1],1)
                            ls.append(k1)
                        if min(ls)<param_DEgene_var[0] and max(ls)>param_DEgene_var[0] and (max(ls)-min(ls)>0.2 or max(ls)-min(ls)>param_DEgene_var[1]):
                            break
                    l = []
                    for i1 in range(len(ls7_1)):
                        l.append(float(np.exp(ls[Group_label[i1]-1])))
                    p1 = random.random()
                    if p1 < param_marker/param_DEgene:
                        l_marker.append(i)
                        l2 = l
                        l1 = max(l)
                        for i1 in range(len(l2)):
                            if l2[i1] != l1:
                                l[i1] = 0.0001
                DEgene_factor.append(l)
            if tree == 1:

                for i1 in range(len(ls7_1)):
                    if type(param_group[Group_label[i1]-1]) == type([]):
                        p = random.random()
                        if p < 1.1:
                            ls = []
                            for i2 in range(len(param_group[Group_label[i1]-1])):
                                ls.append(param_group[Group_label[i1]-1][i2]/sum(param_group[Group_label[i1]-1]))

                            ls1 = [0]
                            for i2 in range(len(ls)):
                                ls1.append(sum(ls[:i2+1]))

                            p1 = random.random()
                            for i2 in range(len(ls1)-1):
                                if p1 > ls1[i2] and p1 <= ls1[i2+1]:
                                    break
                            Group_label[i1] = [Group_label[i1],i2+1]



                for i in range(len(ls4)):
                    if DEgene_label[i] != 0:
                        ls1 = []
                        dict1 = {}
                        for i1 in range(len(param_group)):
                            if type(param_group[i1]) == type([]):
                                ls = []
                                for i2 in range(len(param_group[i1])):
                                    k1 = np.random.normal(param_DEgene_var[0],0.8*param_DEgene_var[1],1)
                                    ls.append(float(np.exp(k1)))
                                ls1.append(ls)
                                dict1[i1] = ls

                        for i1 in range(len(ls7_1)):
                            if type(Group_label[i1]) == type([]):
                                DEgene_factor[i][i1] = DEgene_factor[i][i1]*dict1[Group_label[i1][0]-1][Group_label[i1][1]-1]



        else:
            l1a = []
            l2a = []
            for i in range(len(ls7_1)):
                p1 = random.random()
                Path.append(float(p1))
            ls = []
            while 1:
                p1 = random.randint(0,len(ls4)-1)
                if DEgene_label[p1] == 0 and p1 not in ls:
                    ls.append(p1)
                if len(ls) >= param_noline*(len(ls4)+len(ls2_0)):
                    break
            for i in range(len(ls4)):
                if i in ls:
                    Noline_gene.append(1)
                else:
                    Noline_gene.append(0)
            for i in range(len(ls4)):
                if DEgene_label[i] == 0 and Noline_gene[i] == 0:
                    l = [1 for i1 in range(len(ls7_1))]
                    l1a.append(0)
                    l2a.append(0)
                else:
                    l = []
                    if DEgene_label[i] == 1:
                        while 1:
                            ls = []
                            for i1 in range(len(param_group)):
                                k1 = np.random.normal(param_DEgene_var[0],param_DEgene_var[1],1)
                                ls.append(k1)
                            if min(ls)<param_DEgene_var[0] and max(ls)>param_DEgene_var[0] and (max(ls)-min(ls)>0.2 or max(ls)-min(ls)>param_DEgene_var[1]):
                                break
                        for i1 in range(len(ls7_1)):
                            l.append((float(np.exp(ls[Group_label[i1]-1]))-1)*Path[i1]+1)
                        l1a.append(list(np.exp(ls)))
                        l2a.append(0)
                    else:
                        ls = []
                        for i in range(len(param_group)):
                            k1 = np.random.normal(param_DEgene_var[0],0.5*param_DEgene_var[1],1)
                            ls.append(float(np.exp(k1)))
                        K1 = []
                        K2 = []
                        K3 = []
                        for i1 in range(len(ls)):
                            t = []
                            w = []
                            t.append(0)
                            w.append(0)
                            for i in range(50):
                                t.append((i+1)*0.02)
                                k1 = np.random.normal(0,0.01)
                                w.append(w[i]+k1)
                            b = []
                            for i in range(len(t)):
                                b.append(w[i]+1-t[i]*(w[-1])+t[i]*(ls[i1]-1))
                            p2 = ls[i1]
                            def funa(x,k1,k2,k3):
                                return k1*(x-x**2)+k2*(x**3-x**4)+k3*(x**5-x**6)+1+(p2-1)*x

                            x = np.array(t)
                            y = np.array(b)
                            popt, pcov = curve_fit(funa,x,y)
                            k1 = popt[0]
                            k2 = popt[1]
                            k3 = popt[2]
                            K1.append(k1)
                            K2.append(k2)
                            K3.append(k3)
                        for i1 in range(len(ls7_1)):
                            x1 = Path[i1]
                            if K1[Group_label[i1]-1]*(x1-x1**2)+K2[Group_label[i1]-1]*(x1**3-x1**4)+K3[Group_label[i1]-1]*(x1**5-x1**6)+1+(ls[Group_label[i1]-1]-1)*x1 > 0:
                                l.append(K1[Group_label[i1]-1]*(x1-x1**2)+K2[Group_label[i1]-1]*(x1**3-x1**4)+K3[Group_label[i1]-1]*(x1**5-x1**6)+1+(ls[Group_label[i1]-1]-1)*x1)
                            else:
                                l.append(0)
                        l1a.append(0)
                        l2a.append(list(ls))
                DEgene_factor.append(l)

            if tree == 1:
                for i1 in range(len(ls7_1)):
                    if type(param_group[Group_label[i1]-1]) == type([]):
                        p = random.random()
                        if p < 0.5:
                            ls = []
                            for i2 in range(len(param_group[Group_label[i1]-1])):
                                ls.append(param_group[Group_label[i1]-1][i2]/sum(param_group[Group_label[i1]-1]))

                            ls1 = [0]
                            for i2 in range(len(ls)):
                                ls1.append(sum(ls[:i2+1]))

                            p1 = random.random()
                            for i2 in range(len(ls1)-1):
                                if p1 > ls1[i2] and p1 <= ls1[i2+1]:
                                    break
                            Group_label[i1] = [Group_label[i1],i2+1]
                dictx1 = {}
                i1 = 0
                for i in range(len(param_group)):
                    if type(param_group[i]) == type([]):
                        dictx1[i+1] = i1
                        i1 += 1

                ls1b = []
                for i in range(len(param_group)):
                    if type(param_group[i]) == type([]):
                        ls1b.append(i)

                for i in range(len(ls4)):
                    if DEgene_label[i] != 0:
                        ls1 = []
                        dict1 = {}
                        for i1 in range(len(param_group)):
                            if type(param_group[i1]) == type([]):
                                ls = []
                                for i2 in range(len(param_group[i1])):
                                    k1 = np.random.normal(param_DEgene_var[0],0.8*param_DEgene_var[1],1)
                                    ls.append(float(np.exp(k1)))
                                ls1.append(ls)
                                dict1[i1] = ls
                        for i1 in range(len(ls7_1)):
                            if type(Group_label[i1]) == type([]):
                                start = l1a[i][Group_label[i1][0]-1]
                                end = l1a[i][Group_label[i1][0]-1]*dict1[Group_label[i1][0]-1][Group_label[i1][1]-1]
                                DEgene_factor[i][i1] = float(start+(end-start)*Path[i1])
                    if Noline_gene[i] != 0:
                        KA = []
                        KB = []
                        KC = []
                        for i1 in range(len(param_group)):
                            if type(param_group[i1]) == type([]):
                                ls = []
                                for i2 in range(len(param_group[i1])):
                                    k1 = np.random.normal(param_DEgene_var[0],0.8*param_DEgene_var[1],1)
                                    ls.append(float(np.exp(k1)))
                                K1 = []
                                K2 = []
                                K3 = []
                                for i2 in range(len(ls)):
                                    t = []
                                    w = []
                                    t.append(0)
                                    w.append(0)
                                    for i3 in range(50):
                                        t.append((i3+1)*0.02)
                                        k1 = np.random.normal(0,0.02)
                                        w.append(w[i3]+k1)
                                    b = []
                                    for i3 in range(len(t)):
                                        b.append(w[i3]+l2a[i][i1]-t[i3]*(w[-1])+t[i3]*(ls[i2]*l2a[i][i1]-l2a[i][i1]))
                                    p2 = ls[i2]*l2a[i][i1]
                                    def funa(x,k1,k2,k3):
                                        return k1*(x-x**2)+k2*(x**3-x**4)+k3*(x**5-x**6)+l2a[i][i1]+(p2-l2a[i][i1])*x

                                    x = np.array(t)
                                    y = np.array(b)
                                    popt, pcov = curve_fit(funa,x,y)
                                    k1 = popt[0]
                                    k2 = popt[1]
                                    k3 = popt[2]
                                    K1.append(k1)
                                    K2.append(k2)
                                    K3.append(k3)
                                KA.append(K1)
                                KB.append(K2)
                                KC.append(K3)
                        for i1 in range(len(ls7_1)):
                            if type(Group_label[i1]) == type([]):
                                x1 = Path[i1]
                                if KA[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1-x1**2)+KB[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**3-x1**4)+KC[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**5-x1**6)+l2a[i][Group_label[i1][0]-1]+(ls[Group_label[i1][1]-1]*l2a[i][Group_label[i1][0]-1]-l2a[i][Group_label[i1][0]-1])*x1 > 0:
                                    DEgene_factor[i][i1] = KA[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1-x1**2)+KB[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**3-x1**4)+KC[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**5-x1**6)+l2a[i][Group_label[i1][0]-1]+(ls[Group_label[i1][1]-1]*l2a[i][Group_label[i1][0]-1]-l2a[i][Group_label[i1][0]-1])*x1
                                else:
                                    DEgene_factor[i][i1] = 0


                for i1 in range(len(ls7_1)):
                    if type(Group_label[i1]) == type([]):
                        Path[i1] = 0.5+Path[i1]/2
                    else:
                        if Group_label[i1]-1 in ls1b:
                            ls = []
                            for i2 in range(len(param_group[Group_label[i1]-1])):
                                ls.append(param_group[Group_label[i1]-1][i2]/sum(param_group[Group_label[i1]-1]))

                            ls1 = [0]
                            for i2 in range(len(ls)):
                                ls1.append(sum(ls[:i2+1]))

                            p1 = random.random()
                            for i2 in range(len(ls1)-1):
                                if p1 > ls1[i2] and p1 <= ls1[i2+1]:
                                    break
                            Group_label[i1] = [Group_label[i1],i2+1]
                            Path[i1] = Path[i1]/2







    ls7a = []
    ls7b = []
    ls1 = np.array(ls4)
    ls1a = np.array(ls12c)
    ls1b = np.array(lsk_a)
    ls1c = np.array(lsk_b)
    for i in range(len(ls7)):
        g = np.exp(gammaln(1/ls1a+1)-gammaln(1/ls1a))
        ls7a.append(sum(ls1a*ls1*g*(1-ls1c*(1+ls1b*ls1*ls7[i]*ls1a)**(-1/ls1a-1))))
    ls7b = [i/np.median(ls7a) for i in ls7a]
    Lj = np.median(ls7a)




    Lj1 = sum(ls4)


    lm = []
    ls4a = []
    ls4b = []
    s = np.array(ls7)
    for i in range(len(ls4p)):
        ls4a = np.array(ls4p[i])
        ls4b.append(np.mean(np.log2(ls4a/s*(1000000/Lj1)+1)))
    ls4c = []
    for i in range(len(ls4p)):
        ls4a = np.array(ls4p[i])
        ls4c.append(np.var(np.log2(ls4a/s*(1000000/Lj1)+1)))

    lp = []
    lm1 = []
    s = np.array(ls7b)

    k3 = (len(ls4)+len(ls2_0))//100
    i1 = 0
    for i in range(len(ls4)+len(ls2_0)):
        if i%k3 == 0:
            progress(int(i/k3),width=50)
        l = []
        if i in ls2_0:
            for j in range(len(ls7_1)):
                l.append(0)
            if mod == '-d nocopula':
                l2 = []
                l3 = []
                for i2 in range(len(l)):
                    l2.append(int(0))
                    l3.append('uncertain')
                ltruecount.append(l2)
                ldp.append(l3)
            lp.append(0)
            lm1.append(0)
        else:
            if i1 in l_c:
                l4 = []
                l5 = []
                l4 = k_2[dictx[i1]]
                l5 = list(norm.cdf(l4))
                for j in range(len(ls7_1)):
                    j1 = 0
                    if param_batch == [1]:
                        p1 = float((1+ls12c[i1]*ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[j])**(-1/ls12c[i1])*(1-lsk_b[i1]*(lsk_a[i1]/(1/(ls12c[i1]*ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[j])+1)+1)**(-1/ls12c[i1]))+(1+lsk_a[i1]*ls12c[i1]*ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[j])**(-1/ls12c[i1])*lsk_b[i1])
                    else:
                        p1 = float((1+ls12c[i1]*ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[i1][j])**(-1/ls12c[i1])*(1-lsk_b[i1]*(lsk_a[i1]/(1/(ls12c[i1]*ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[i1][j])+1)+1)**(-1/ls12c[i1]))+(1+lsk_a[i1]*ls12c[i1]*ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[i1][j])**(-1/ls12c[i1])*lsk_b[i1])

                    j2 = p1
                    if l5[j] <= p1:
                        i4 = 0
                    else:
                        while 1:
                            ls = []
                            for j3 in range(500):
                                ls.append(int(j3+1+500*j1))
                            j1 += 1
                            if param_batch == [1]:
                                ls1 = integral(ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[j],ls12c[i1],ls,lsk_a[i1],float(-np.log(lsk_b[i1])))
                            else:
                                ls1 = integral(ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[i1][j],ls12c[i1],ls,lsk_a[i1],float(-np.log(lsk_b[i1])))

                            for j3 in range(len(ls1)):
                                j5 = j2
                                j2 += ls1[j3]
                                if j2 > l5[j] or (j2 > 0.99 and j5 == j2):
                                    i4 = ls[j3]
                                    break
                            else:
                                continue
                            break
                    l.append(i4)
                if cell_number == len(ls7):
                    lp.append(abs(l.count(0)-ls4p[i1].count(0))/len(l))
                    m1 = np.mean(np.mean(np.log2(l/s*(1000000/Lj)+1)))
                    m2 = ls4b[i1]
                    lm1.append(abs(m1-m2))


                l6 = []
                for i3 in range(len(l)):
                    l6.append(l[i3]/ls7_1[i3])
                if max(l6) > 0:
                    lm.append(l6)



            else:
                L = []
                L2 = []
                L3 = []
                if cell_number == len(ls7):
                    for j1 in range(20):
                        ls1 = np.random.gamma(1/ls12c[i1],ls4[i1]*ls12c[i1],len(ls4p[0]))
                        l1 = []
                        for j in range(len(ls1)):
                            if param_batch == [1]:
                                l1.append(float(ls1[j]*ls7[j]*batch_factor[j]*DEgene_factor[i1][j]))
                            else:
                                l1.append(float(ls1[j]*ls7[j]*batch_factor[i1][j]*DEgene_factor[i1][j]))
                        c = 0
                        while 1:
                            c += 1
                            x = np.random.poisson(tuple(l1),(1,len(ls4p[0])))
                            ls1 = list(map(int,x[0]))
                            if max(ls1) > 0 or c >10:
                                break
                        l2 = ls1


                        if mod == '-d nocopula':
                            c = 0
                            while 1:
                                c += 1
                                l = []
                                l3 = []
                                for i2 in range(len(l1)):
                                    m = randomm(np.exp(-lsk_a[i1]*l1[i2])*lsk_b[i1])
                                    l.append(int(l2[i2]*m))
                                    if m == 0:
                                        l3.append('true')
                                    else:
                                        l3.append('false')
                                if max(l) > 0 or c >10:
                                    break
                            L2.append(l2)
                            L3.append(l3)
                            L.append(l)
                            #ltruecount.append(l2)
                            #ldp.append(l3)
                        else:
                            c = 0
                            while 1:
                                c += 1
                                l = []
                                for i2 in range(len(l1)):
                                    m = randomm(np.exp(-lsk_a[i1]*l1[i2])*lsk_b[i1])
                                    l.append(int(l2[i2]*m))
                                if max(l) > 0 or c >10:
                                    break
                            L.append(l)
                    La = [i3.count(0)/len(i3) for i3 in L]
                    Lb = []
                    Lb1 = []
                    s = np.array(ls7b)

                    for i3 in range(len(L)):
                        l = np.array(L[i3])
                        Lb.append(np.mean(np.log2(l/s*(1000000/Lj)+1)))
                        Lb1.append(np.var(np.log2(l/s*(1000000/Lj)+1)))
                    Lc = []


                    for i3 in range(len(L)):
                        p0 = ls4p[i1].count(0)/len(ls4p[i1])
                        if p0 > 0:
                            Lc.append((abs(La[i3]-p0)/p0)+abs(Lb[i3]-ls4b[i1])/ls4b[i1])
                        else:
                            if ls4c[i1] > 0:
                                Lc.append(abs(Lb[i3]-ls4b[i1])/ls4b[i1]+(abs(Lb1[i3]-ls4c[i1])/ls4c[i1])**2)
                            else:
                                Lc.append(abs(Lb[i3]-ls4b[i1])/ls4b[i1])
                    i3 = np.argmin(Lc)
                    if mod == '-d nocopula':
                        ltruecount.append(L2[i3])
                        ldp.append(L3[i3])
                        l = L[i3]
                    else:
                        l = L[i3]
                    lp.append(abs(l.count(0)-ls4p[i1].count(0))/len(l))
                    m1 = np.mean(np.mean(np.log2(l/s*(1000000/Lj)+1)))
                    m2 = ls4b[i1]
                    lm1.append(abs(m1-m2))
                else:
                    ls1 = np.random.gamma(1/ls12c[i1],ls4[i1]*ls12c[i1],len(ls7_1))
                    l1 = []
                    for j in range(len(ls1)):
                        if param_batch == [1]:
                            l1.append(float(ls1[j]*ls7_1[j]*batch_factor[j]*DEgene_factor[i1][j]))
                        else:
                            l1.append(float(ls1[j]*ls7_1[j]*batch_factor[i1][j]*DEgene_factor[i1][j]))
                    c = 0
                    while 1:
                        c += 1
                        x = np.random.poisson(tuple(l1),(1,len(ls7_1)))
                        ls1 = list(map(int,x[0]))
                        if max(ls1) > 0 or c >10:
                            break
                    l2 = ls1



                    if mod == '-d nocopula':
                        c = 0
                        while 1:
                            c += 1
                            l = []
                            l3 = []
                            for i2 in range(len(l1)):
                                m = randomm(np.exp(-lsk_a[i1]*l1[i2])*lsk_b[i1])
                                l.append(int(l2[i2]*m))
                                if m == 0:
                                    l3.append('true')
                                else:
                                    l3.append('false')
                            if max(l) > 0 or c >10:
                                break
                        ltruecount.append(l2)
                        ldp.append(l3)

            i1 += 1
        lcount.append(l)

    if param_group == [1] and param_batch == [1] and cell_number == len(ls7):
        s1 = np.array(ls7b)
        la = np.argsort(lp)
        lb = np.sort(lm1)
        i1 = 0
        for i in range(len(ls4)+len(ls2_0)):
            if i not in ls2_0:
                if i in la[-int(len(ls4)*0.05):] or i in lb[-int(len(ls4)*0.06):]:
                    if i in lb[-int(len(ls4)*0.06):]:
                        l = []
                        l1 = np.random.normal(0,ls4[i1]*0.05,10)
                        l.extend(list(l1))
                        l1 = np.random.normal(0,ls4[i1]*0.01,20)
                        l.extend(list(l1))
                        L = lm1[i]
                        for i2 in range(len(l)):
                            if ls4[i1]+l[i2] > 0:
                                lc = []
                                if i1 in l_c:
                                    l4 = []
                                    l5 = []
                                    l4 = k_2[dictx[i1]]
                                    l5 = list(norm.cdf(l4))
                                    for j in range(len(ls4p[0])):
                                        j1 = 0
                                        if param_batch == [1]:
                                            p1 = float((1+ls12c[i1]*(ls4[i1]+l[i2])*ls7[j]*DEgene_factor[i1][j]*batch_factor[j])**(-1/ls12c[i1])*(1-lsk_b[i1]*(lsk_a[i1]/(1/(ls12c[i1]*(ls4[i1]+l[i2])*ls7[j]*DEgene_factor[i1][j]*batch_factor[j])+1)+1)**(-1/ls12c[i1]))+(1+lsk_a[i1]*ls12c[i1]*(ls4[i1]+l[i2])*ls7[j]*DEgene_factor[i1][j]*batch_factor[j])**(-1/ls12c[i1])*lsk_b[i1])
                                        else:
                                            p1 = float((1+ls12c[i1]*(ls4[i1]+l[i2])*ls7[j]*DEgene_factor[i1][j]*batch_factor[i1][j])**(-1/ls12c[i1])*(1-lsk_b[i1]*(lsk_a[i1]/(1/(ls12c[i1]*(ls4[i1]+l[i2])*ls7[j]*DEgene_factor[i1][j]*batch_factor[i1][j])+1)+1)**(-1/ls12c[i1]))+(1+lsk_a[i1]*ls12c[i1]*(ls4[i1]+l[i2])*ls7[j]*DEgene_factor[i1][j]*batch_factor[i1][j])**(-1/ls12c[i1])*lsk_b[i1])

                                        j2 = p1
                                        if l5[j] <= p1:
                                            i4 = 0
                                        else:
                                            while 1:
                                                ls = []
                                                for j3 in range(500):
                                                    ls.append(int(j3+1+500*j1))
                                                j1 += 1
                                                if param_batch == [1]:
                                                    ls1 = integral((ls4[i1]+l[i2])*ls7[j]*DEgene_factor[i1][j]*batch_factor[j],ls12c[i1],ls,lsk_a[i1],float(-np.log(lsk_b[i1])))
                                                else:
                                                    ls1 = integral((ls4[i1]+l[i2])*ls7[j]*DEgene_factor[i1][j]*batch_factor[i1][j],ls12c[i1],ls,lsk_a[i1],float(-np.log(lsk_b[i1])))

                                                for j3 in range(len(ls1)):
                                                    j5 = j2
                                                    j2 += ls1[j3]
                                                    if j2 > l5[j] or (j2 > 0.99 and j5 == j2):
                                                        i4 = ls[j3]
                                                        break
                                                else:
                                                    continue
                                                break
                                        lc.append(i4)
                                    if abs(np.mean(np.log2(np.array(lc)/s1*(1000000/Lj)+1))-ls4b[i1]) < L:
                                        L = abs(np.mean(np.log2(np.array(lc)/s1*(1000000/Lj)+1))-ls4b[i1])
                                        lcount[i] = lc
                                        ls4[i1] = ls4[i1]+l[i2]
                                else:
                                    ls1 = np.random.gamma(1/ls12c[i1],(ls4[i1]+l[i2])*ls12c[i1],len(ls4p[0]))
                                    l1 = []
                                    for j in range(len(ls1)):
                                        if param_batch == [1]:
                                            l1.append(float(ls1[j]*ls7[j]*batch_factor[j]*DEgene_factor[i1][j]))
                                        else:
                                            l1.append(float(ls1[j]*ls7[j]*batch_factor[i1][j]*DEgene_factor[i1][j]))
                                    c = 0
                                    while 1:
                                        c += 1
                                        x = np.random.poisson(tuple(l1),(1,len(ls4p[0])))
                                        ls1 = list(map(int,x[0]))
                                        if max(ls1) > 0 or c >10:
                                            break
                                    l2 = ls1


                                    if mod == '-d nocopula':
                                        c = 0
                                        while 1:
                                            c += 1
                                            lc = []
                                            l3 = []
                                            for i2 in range(len(l1)):
                                                m = randomm(np.exp(-lsk_a[i1]*l1[i2])*lsk_b[i1])
                                                lc.append(int(l2[i2]*m))
                                                if m == 0:
                                                    l3.append('true')
                                                else:
                                                    l3.append('false')
                                            if max(l) > 0 or c >10:
                                                break
                                        if abs(np.mean(np.log2(np.array(lc)/s1*(1000000/Lj)+1))-ls4b[i1]) < L:
                                            L = abs(np.mean(np.log2(np.array(lc)/s1*(1000000/Lj)+1))-ls4b[i1])
                                            lcount[i] = lc
                                            ls4[i1] = ls4[i1]+l[i2]
                                            ltruecount[i] = l2
                                            ldp[i] = l3
                                        #ltruecount.append(l2)
                                        #ldp.append(l3)
                                    else:
                                        c = 0
                                        while 1:
                                            c += 1
                                            lc = []
                                            for i2 in range(len(l1)):
                                                m = randomm(np.exp(-lsk_a[i1]*l1[i2])*lsk_b[i1])
                                                lc.append(int(l2[i2]*m))
                                            if max(l) > 0 or c >10:
                                                break
                                        if abs(np.mean(np.log2(np.array(lc)/s1*(1000000/Lj)+1))-ls4b[i1]) < L:
                                            L = abs(np.mean(np.log2(np.array(lc)/s1*(1000000/Lj)+1))-ls4b[i1])
                                            lcount[i] = lc
                                            ls4[i1] = ls4[i1]+l[i2]
                    else:
                        l = []
                        l1 = np.random.normal(0,ls12c[i1]*0.05,10)
                        l.extend(list(l1))
                        l1 = np.random.normal(0,ls12c[i1]*0.01,20)
                        l.extend(list(l1))
                        L = lp[i]
                        for i2 in range(len(l)):
                            if ls12c[i1]+l[i2] > 0:
                                lc = []
                                if i1 in l_c:
                                    l4 = []
                                    l5 = []
                                    l4 = k_2[dictx[i1]]
                                    l5 = list(norm.cdf(l4))
                                    for j in range(len(ls4p[0])):
                                        j1 = 0
                                        if param_batch == [1]:
                                            p1 = float((1+(ls12c[i1]+l[i2])*ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j])**(-1/(ls12c[i1]+l[i2]))*(1-lsk_b[i1]*(lsk_a[i1]/(1/((ls12c[i1]+l[i2])*ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j])+1)+1)**(-1/(ls12c[i1]+l[i2])))+(1+lsk_a[i1]*(ls12c[i1]+l[i2])*ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j])**(-1/(ls12c[i1]+l[i2]))*lsk_b[i1])
                                        else:
                                            p1 = float((1+(ls12c[i1]+l[i2])*ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[i1][j])**(-1/(ls12c[i1]+l[i2]))*(1-lsk_b[i1]*(lsk_a[i1]/(1/((ls12c[i1]+l[i2])*ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[i1][j])+1)+1)**(-1/(ls12c[i1]+l[i2])))+(1+lsk_a[i1]*(ls12c[i1]+l[i2])*ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[i1][j])**(-1/(ls12c[i1]+l[i2]))*lsk_b[i1])

                                        j2 = p1
                                        if l5[j] <= p1:
                                            i4 = 0
                                        else:
                                            while 1:
                                                ls = []
                                                for j3 in range(500):
                                                    ls.append(int(j3+1+500*j1))
                                                j1 += 1
                                                if param_batch == [1]:
                                                    ls1 = integral(ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j],(ls12c[i1]+l[i2]),ls,lsk_a[i1],float(-np.log(lsk_b[i1])))
                                                else:
                                                    ls1 = integral(ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[i1][j],(ls12c[i1]+l[i2]),ls,lsk_a[i1],float(-np.log(lsk_b[i1])))

                                                for j3 in range(len(ls1)):
                                                    j5 = j2
                                                    j2 += ls1[j3]
                                                    if j2 > l5[j] or (j2 > 0.99 and j5 == j2):
                                                        i4 = ls[j3]
                                                        break
                                                else:
                                                    continue
                                                break
                                        lc.append(i4)
                                    if abs(lc.count(0)/len(lc)-ls4p[i1].count(0)/len(lc)) < L:
                                        L = abs(lc.count(0)/len(lc)-ls4p[i1].count(0)/len(lc))
                                        lcount[i] = lc
                                        ls12c[i1] = ls12c[i1]+l[i2]
                                else:
                                    ls1 = np.random.gamma(1/(ls12c[i1]+l[i2]),ls4[i1]*(ls12c[i1]+l[i2]),len(ls4p[0]))
                                    l1 = []
                                    for j in range(len(ls1)):
                                        if param_batch == [1]:
                                            l1.append(float(ls1[j]*ls7[j]*batch_factor[j]*DEgene_factor[i1][j]))
                                        else:
                                            l1.append(float(ls1[j]*ls7[j]*batch_factor[i1][j]*DEgene_factor[i1][j]))
                                    c = 0
                                    while 1:
                                        c += 1
                                        x = np.random.poisson(tuple(l1),(1,len(ls4p[0])))
                                        ls1 = list(map(int,x[0]))
                                        if max(ls1) > 0 or c >10:
                                            break
                                    l2 = ls1


                                    if mod == '-d nocopula':
                                        c = 0
                                        while 1:
                                            c += 1
                                            lc = []
                                            l3 = []
                                            for i3 in range(len(l1)):
                                                m = randomm(np.exp(-lsk_a[i1]*l1[i3])*lsk_b[i1])
                                                lc.append(int(l2[i3]*m))
                                                if m == 0:
                                                    l3.append('true')
                                                else:
                                                    l3.append('false')
                                            if max(l) > 0 or c >10:
                                                break
                                        if abs(lc.count(0)/len(lc)-ls4p[i1].count(0)/len(lc)) < L:
                                            L = abs(lc.count(0)/len(lc)-ls4p[i1].count(0)/len(lc))
                                            lcount[i] = lc
                                            ls12c[i1] = ls12c[i1]+l[i2]
                                            ltruecount[i] = l2
                                            ldp[i] = l3
                                        #ltruecount.append(l2)
                                        #ldp.append(l3)
                                    else:
                                        c = 0
                                        while 1:
                                            c += 1
                                            lc = []
                                            for i3 in range(len(l1)):
                                                m = randomm(np.exp(-lsk_a[i1]*l1[i3])*lsk_b[i1])
                                                lc.append(int(l2[i3]*m))
                                            if max(l) > 0 or c >10:
                                                break
                                        if abs(lc.count(0)/len(lc)-ls4p[i1].count(0)/len(lc)) < L:
                                            L = abs(lc.count(0)/len(lc)-ls4p[i1].count(0)/len(lc))
                                            lcount[i] = lc
                                            ls12c[i1] = ls12c[i1]+l[i2]




                i1+=1




    progress(100,width=50)
    print(' ')
    print('Simulation completed')
    if mod != '-d nocopula':
        Mat = np.array(lm)
        correlation1 = np.corrcoef(Mat)
        
        #sns.set()
        #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        #plt.show()

        ln = []
        for i in range(len(lm)):
            ln.append(float(np.mean(lm[i])))
        ln1 = []
        ln1 = list(reversed(list(np.argsort(ln))))
        ln2 = []
        for i in range(len(lm)):
            if i in ln1[:50]:
                ln2.append(lm[i])
        Mat = np.array(ln2)
        correlation1 = np.corrcoef(Mat)
        
        #sns.set()
        #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        #plt.show()


    



    

    
def mean_noparam():
    """
    Estimate mean and dispersion of -c or -c nocopula

    Parameters
    ----------
    ls4p:list
          Normalized count

    Returns
    -------
    ls4:list
          Mean of gene
    ls12c:list
          Dispersion of gene

    """
    print('Estimate mean and dispersion...')
    global ls4,ls4p         
    ls4 = []
    ls4p = []
    for i in range(len(ls3)-1):
        ls4p.append(ls3[i+1][1:])
    for i in range(len(ls4p)):
        sum1 = 0
        k = 0
        for j in range(len(ls4p[i])):
            sum1 += ls4p[i][j]
            k += 1
        b = sum1/k
        ls4.append(b)

def BCV_noparam():
    """
    Estimate mean and dispersion of -c or -c nocopula

    Parameters
    ----------
    ls4p:list
          Normalized count

    Returns
    -------
    ls4:list
          Mean of gene
    ls12c:list
          Dispersion of gene
    ls7:list
          Size factor

    """
    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')


    from scipy.special import gammaln
    from scipy.special import digamma
    global ls2,ls3,ls7,ls4p,ls4,ls12c,l_c,correlation,dictx
    l_c = []
    ls12 = []
    dictx  ={}
    for i in range(len(ls4p)):
        ls12.append(float(np.var(ls4p[i])/(ls4[i]**2)-1/ls4[i]))
    ls12a = []
    ls12b = []
    ls12d = []
    i1 = 0
    i2 = 0

    L = len(ls4)//100
    for i in range(len(ls2)-1):
        if i%L == 0:
            progress(i//L,width=50)
            #print('{}%'.format(i//L))

        i2 = 0
        i3 = np.mean(ls4p[i])
        i4 = np.var(ls4p[i])
        for j in range(len(ls4p[i])):
            i2 += 1/(ls7[j]*i3)

        da = i4/(i3**2)-i2/len(ls7)
        if da > 0:
            ls12b.append(da)
        else:
            ls12b.append(0.001)

        
        if ls12[i] <= 0:
            d1 = 100
        else:
            d1 = ls12[i]
        mu = ls4[i]


        l = list(np.random.normal(0,d1/3,15))
        l1 = list(np.random.normal(0,d1/9,15))
        l.extend(l1)

        k = 0
        for j in range(len(ls2[i+1])-1):
            if float(ls2[i+1][j+1]) == 0:
                k += 1
        
        k1 = 0
        for j in range(len(ls7)):
            k1 += ((d1**(-1)/(ls7[j]*mu+d1**(-1)))**(1/d1))
        k2 = abs(k1-k)
        d = d1
        for j in l:
            if d+j > 0:
                d2 = d+j
            k1 = 0
            for j in range(len(ls7)):
                k1 += ((d2**(-1)/(ls7[j]*mu+d2**(-1)))**(1/d2))
            if abs(k1-k) < k2:
                d = d2
                k2 = abs(k1-k)
        ls12d.append(d)


        
        
        la = []
        lb = []
        lc = []
        ld = []
        for j in range(len(ls7)):
            la.append(mu*ls7[j]*d1)
            lb.append(1/(mu*ls7[j]*d1)+1)
            lc.append(1/d1+ls4p[i][j])
            ld.append(1/d1)
        la1 = np.log(la)
        lb1 = np.log(lb)
        lc1 = digamma(lc)
        ld1 = digamma(ld)
        la1 = pd.Series(la1)
        lb1 = pd.Series(lb1)
        lc1 = pd.Series(lc1)
        ld1 = pd.Series(ld1)
        s1 = pd.Series(ls4p[i])
        s2 = pd.Series(ls7)
        #lnF = 0
        d12 = d1**2
        d13 = d1**3
        s = 1/d12*la1-1/d12+1/d12*lb1+(1+s1*d1)/(d12+mu*s2*d13)-1/d12*lc1+1/d12*ld1
        lnF = float(sum(s))

        
        #for j in range(len(ls7)):
        #    lnF += 1/d12*la1[j]-1/d12+1/d12*lb1[j]+(1+ls4p[i][j]*d1)/(d12+mu*ls7[j]*d13)-1/d12*lc1[j]+1/d12*ld1[j]

        if lnF < 0:
            K1 = 0
            K2 = d1
        else:
            K1 = d1
            K2 = 1000
        for j in range(20):
            d1 = (K1+K2)/2
            la = []
            lb = []
            lc = []
            ld = []
            for i1 in range(len(ls7)):
                la.append(mu*ls7[i1]*d1)
                lb.append(1/(mu*ls7[i1]*d1)+1)
                lc.append(1/d1+ls4p[i][i1])
                ld.append(1/d1)
            la1 = np.log(la)
            lb1 = np.log(lb)
            lc1 = digamma(lc)
            ld1 = digamma(ld)
            la1 = pd.Series(la1)
            lb1 = pd.Series(lb1)
            lc1 = pd.Series(lc1)
            ld1 = pd.Series(ld1)
            s1 = pd.Series(ls4p[i])
            s2 = pd.Series(ls7)
            #lnF = 0
            d12 = d1**2
            d13 = d1**3
            s = 1/d12*la1-1/d12+1/d12*lb1+(1+s1*d1)/(d12+mu*s2*d13)-1/d12*lc1+1/d12*ld1
            lnF = float(sum(s))
            
            #for i1 in range(len(ls7)):
            #    lnF += 1/d12*la1[i1]-1/d12+1/d12*lb1[i1]+(1+ls4p[i][i1]*d1)/(d12+mu*ls7[i1]*d13)-1/d12*lc1[i1]+1/d12*ld1[i1]
            if lnF > 0:
                K1 = d1
            else:
                K2 = d1
        
        ls12a.append((K1+K2)/2)


    ls4 = np.log(ls4)

    # plt.scatter(ls4,ls12a,s = 1,c = 'b')
    # plt.scatter(ls4,ls12b,s = 1,c = 'r')
    # plt.show()
        
    ls = sorted(ls4)
    l = (ls[-1]+0.001-ls[0])/20
    l1 = []
    l2 = []
    for i in range(20):
        l1.append([])
        l2.append([])
        
    for i in range(len(ls4)):
        l1[int((ls4[i]-ls[0])/l)].append(ls4[i])
        l2[int((ls4[i]-ls[0])/l)].append(ls12a[i])
    l3 = []
    l4 = []
    for i in range(len(l1)):
        l2[i].sort()
        if  len(l1[i]) > 10:
            l3.append(float(np.mean(l1[i])))
            l4.append(l2[i][int(len(l2[i])/2)])

    l5 = []
    l6 = []
    for i in range(20):
        l5.append([])
        l6.append([])
    for i in range(len(ls4)):
        l5[int((ls4[i]-ls[0])/l)].append(ls4[i])
        l6[int((ls4[i]-ls[0])/l)].append(ls12b[i])
    l7 = []
    l8 = []
    for i in range(len(l5)):
        l6[i].sort()
        if  len(l1[i]) > 10:
            l7.append(float(np.mean(l5[i])))
            l8.append(l6[i][int(len(l6[i])/2)])





            
    def funb(x,t1,t2,t3,d0):
        return t1/(t2+np.exp(t3*x))+d0

    ls1 = []
    for i in range(len(l4)):
        if l4[i] != 0:
            ls1.append(l8[i]/l4[i])
        else:
            ls1.append(0)
    x = np.array(l3[8:])
    y = np.array(ls1[8:])
    param_bounds=([0,0,-np.inf,0],[np.inf,np.inf,np.inf,np.inf])
    try:
        popt, pcov = curve_fit(funb,x,y,bounds=param_bounds)
    except:
        param_bounds=([0,0,-100,0],[500,500,100,200])
        popt, pcov = curve_fit(funb,x,y,bounds=param_bounds)
    T1 = popt[0]
    T2 = popt[1]
    T3 = popt[2]
    D0 = popt[3]

    #x1 = np.linspace(l3[0],l3[-1],600000,endpoint=True)
    #y1 = T1/(T2+np.exp(T3*x1))+D0
    #plt.plot(x1,y1)
    #plt.scatter(l3,ls1)
    #plt.show()



    ls12c = []
    for i in range(len(ls4)):
        ls12c.append((T1/(T2+np.exp(T3*ls4[i]))+D0)*ls12a[i])








    l2 = []
    l3 = []
    
    l2 = list(reversed(list(np.argsort(ls4))))

    if mod != '-c nocopula':
        if len(ls4) > copula_number:
            for i in range(copula_number):
                l3.append(ls4p[l2[i]])
                l_c.append(l2[i])
            for i in range(len(l_c)):
                dictx[l_c[i]] = i
        else:
            for i in range(len(ls4)):
                l3.append(ls4p[l2[i]])
                l_c.append(l2[i])
            for i in range(len(l_c)):
                dictx[l_c[i]] = i
        warnings.filterwarnings("ignore")
        Mat = np.array(l3)
        correlation = np.corrcoef(Mat)
        correlation1 = np.corrcoef(Mat)
        #sns.set()
        #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        #plt.show()
        
        Mat = np.array(l3[:50])
        correlation1 = np.corrcoef(Mat)
        #sns.set()
        #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        #plt.show()

    ls4 = list(np.exp(ls4))

def creatcount_noparam():
    """

    Simulate count

    Parameters
    ----------
    ls4:list
          Mean of gene
    ls12:list
          Dispersion of gene
    ls7:list
          Size factor


    Returns
    -------
    ldp:list
          Simulated dropouted gene
    ldropout:list
          Dropout information
    llisize:list
          Library size
    lcount:list
          Count

    """
    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')

    from scipy.stats import nbinom
    from scipy.stats import norm
    global lcount,ls4,correlation,ls4p,ls12c,ls7,dictx,ls7_1
    print(' ')
    print('Start simulation...')

    ls7_1 = []
    if cell_number == len(ls7):
        ls7_1 = ls7
    else:
        if cell_number < len(ls7):
            ls7_1 = ls7[:cell_number]
        else:
            ls7_2 = np.log(ls7)
            x_mean, x_std = norm.fit(ls7_2)
            ls7_1 = [i for i in ls7]
            while 1:
                k = float(np.exp(np.random.normal(x_mean, x_std,1)))
                ls7_1.append(k)
                if len(ls7_1) == cell_number:
                    break

    lcount = []
    if mod != '-c nocopula':
        l5 = []
        for i in range(len(correlation)):
            l5.append(0)
        p = np.random.multivariate_normal(mean=l5,cov=correlation,size=len(ls7_1))

        k  = list(np.transpose(p))
        #k1 = np.corrcoef(np.array(k))
        # sns.set()
        # ax = sns.clustermap(k1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        # plt.show()
        #
        # k1 = np.corrcoef(np.array(k[:50]))
        # sns.set()
        # ax = sns.clustermap(k1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        # plt.show()

    
    global batch_factor,batch_label
    batch_factor = []
    batch_label = []
    if param_batch == [1] and param_library == 1:
        for i in range(len(ls7_1)):
            batch_factor.append(1)
    else:
        if param_library != 1 and param_batch == [1]:
            for i in range(len(ls7_1)):
                batch_factor.append(param_library)
        else:
            for i in range(len(ls7_1)):
                p1 = random.random()
                ls = [0]
                for j in range(len(param_batch)):
                    ls.append(sum(param_batch[:j+1]))
                for j in range(len(ls)-1):
                    if p1 > ls[j] and p1 < ls[j+1]:
                        break
                batch_label.append(j+1)                 ####batch_label
            for i in range(len(ls4p)):
                ls = []
                for j in range(len(param_batch)):
                    k1 = np.random.normal(param_batch_var[0],param_batch_var[1],1)
                    ls.append(k1)
                ls1 = []
                for j in range(len(ls7_1)):
                    ls1.append(param_library*float(np.exp(ls[batch_label[j]-1])))
                batch_factor.append(ls1)                             #######batch_factor
    global DEgene_factor,Group_label,DEgene_label,Path,Noline_gene,l_marker
    l_marker = []
    DEgene_label = []
    DEgene_factor = []
    Group_label = []
    Path = []
    Noline_gene = []
    if param_group == [1]:
        for i in range(len(ls4)):
            l = [1 for j in range(len(ls7_1))]
            DEgene_factor.append(l)
    else:
        for i in range(len(ls4)):
            p1 = random.random()
            if p1 < param_DEgene*(len(ls4)+len(ls2_0))/len(ls4):
                DEgene_label.append(1)
            else:
                DEgene_label.append(0)              ####DEgene_label

        param_group1 = []
        tree = 0
        for i in range(len(param_group)):
            if type(param_group[i]) == type([]):
                param_group1.append(sum(param_group[i]))
                tree = 1
            else:
                param_group1.append(param_group[i])


        for i in range(len(ls7_1)):
            p1 = random.random()
            ls = [0]
            for j in range(len(param_group1)):
                ls.append(sum(param_group1[:j+1]))
            for j in range(len(ls)-1):
                if p1 > ls[j] and p1 < ls[j+1]:
                    break
            Group_label.append(j+1)                 ####Group_label
        if param_path == 'No':
            for i in range(len(ls4)):
                if DEgene_label[i] == 0:
                    l = [1 for i1 in range(len(ls7_1))]
                else:
                    while 1:
                        ls = []
                        for i1 in range(len(param_group1)):
                            k1 = np.random.normal(param_DEgene_var[0],param_DEgene_var[1],1)
                            ls.append(k1)
                        if min(ls)<param_DEgene_var[0] and max(ls)>param_DEgene_var[0] and (max(ls)-min(ls)>0.2 or max(ls)-min(ls)>param_DEgene_var[1]):
                            break
                    l = []
                    for i1 in range(len(ls7_1)):
                        l.append(float(np.exp(ls[Group_label[i1]-1])))
                    p1 = random.random()
                    if p1 < param_marker/param_DEgene:
                        l_marker.append(i)
                        l2 = l
                        l1 = max(l)
                        for i1 in range(len(l2)):
                            if l2[i1] != l1:
                                l[i1] = 0.0001
                DEgene_factor.append(l)
            if tree == 1:

                for i1 in range(len(ls7_1)):
                    if type(param_group[Group_label[i1]-1]) == type([]):
                        p = random.random()
                        if p < 1.1:
                            ls = []
                            for i2 in range(len(param_group[Group_label[i1]-1])):
                                ls.append(param_group[Group_label[i1]-1][i2]/sum(param_group[Group_label[i1]-1]))

                            ls1 = [0]
                            for i2 in range(len(ls)):
                                ls1.append(sum(ls[:i2+1]))

                            p1 = random.random()
                            for i2 in range(len(ls1)-1):
                                if p1 > ls1[i2] and p1 <= ls1[i2+1]:
                                    break
                            Group_label[i1] = [Group_label[i1],i2+1]



                for i in range(len(ls4)):
                    if DEgene_label[i] != 0:
                        ls1 = []
                        dict1 = {}
                        for i1 in range(len(param_group)):
                            if type(param_group[i1]) == type([]):
                                ls = []
                                for i2 in range(len(param_group[i1])):
                                    k1 = np.random.normal(param_DEgene_var[0],0.8*param_DEgene_var[1],1)
                                    ls.append(float(np.exp(k1)))
                                ls1.append(ls)
                                dict1[i1] = ls

                        for i1 in range(len(ls7_1)):
                            if type(Group_label[i1]) == type([]):
                                DEgene_factor[i][i1] = DEgene_factor[i][i1]*dict1[Group_label[i1][0]-1][Group_label[i1][1]-1]


        else:
            l1a = []
            l2a = []
            for i in range(len(ls7_1)):
                p1 = random.random()
                Path.append(float(p1))
            ls = []
            while 1:
                p1 = random.randint(0,len(ls4)-1)
                if DEgene_label[p1] == 0 and p1 not in ls:
                    ls.append(p1)
                if len(ls) >= param_noline*(len(ls4)+len(ls2_0)):
                    break
            for i in range(len(ls4)):
                if i in ls:
                    Noline_gene.append(1)
                else:
                    Noline_gene.append(0)
            for i in range(len(ls4)):
                if DEgene_label[i] == 0 and Noline_gene[i] == 0:
                    l = [1 for i1 in range(len(ls7_1))]
                    l1a.append(0)
                    l2a.append(0)
                else:
                    l = []
                    if DEgene_label[i] == 1:
                        while 1:
                            ls = []
                            for i1 in range(len(param_group)):
                                k1 = np.random.normal(param_DEgene_var[0],param_DEgene_var[1],1)
                                ls.append(k1)
                            if min(ls)<param_DEgene_var[0] and max(ls)>param_DEgene_var[0] and (max(ls)-min(ls)>0.2 or max(ls)-min(ls)>param_DEgene_var[1]):
                                break
                        for i1 in range(len(ls7_1)):
                            l.append((float(np.exp(ls[Group_label[i1]-1]))-1)*Path[i1]+1)
                        l1a.append(list(np.exp(ls)))
                        l2a.append(0)
                    else:
                        ls = []
                        for i in range(len(param_group)):
                            k1 = np.random.normal(param_DEgene_var[0],param_DEgene_var[1],1)
                            ls.append(float(np.exp(k1)))
                        K1 = []
                        K2 = []
                        K3 = []
                        for i1 in range(len(ls)):
                            t = []
                            w = []
                            t.append(0)
                            w.append(0)
                            for i in range(50):
                                t.append((i+1)*0.02)
                                k1 = np.random.normal(0,0.02)
                                w.append(w[i]+k1)
                            b = []
                            for i in range(len(t)):
                                b.append(w[i]+1-t[i]*(w[-1])+t[i]*(ls[i1]-1))
                            p2 = ls[i1]
                            def funa(x,k1,k2,k3):
                                return k1*(x-x**2)+k2*(x**3-x**4)+k3*(x**5-x**6)+1+(p2-1)*x

                            x = np.array(t)
                            y = np.array(b)
                            popt, pcov = curve_fit(funa,x,y)
                            k1 = popt[0]
                            k2 = popt[1]
                            k3 = popt[2]
                            K1.append(k1)
                            K2.append(k2)
                            K3.append(k3)
                        for i1 in range(len(ls7_1)):
                            x1 = Path[i1]
                            if K1[Group_label[i1]-1]*(x1-x1**2)+K2[Group_label[i1]-1]*(x1**3-x1**4)+K3[Group_label[i1]-1]*(x1**5-x1**6)+1+(ls[Group_label[i1]-1]-1)*x1 > 0:
                                l.append(K1[Group_label[i1]-1]*(x1-x1**2)+K2[Group_label[i1]-1]*(x1**3-x1**4)+K3[Group_label[i1]-1]*(x1**5-x1**6)+1+(ls[Group_label[i1]-1]-1)*x1)
                            else:
                                l.append(0)
                        l1a.append(0)
                        l2a.append(list(ls))
                DEgene_factor.append(l)


            if tree == 1:
                for i1 in range(len(ls7_1)):
                    if type(param_group[Group_label[i1]-1]) == type([]):
                        p = random.random()
                        if p < 0.5:
                            ls = []
                            for i2 in range(len(param_group[Group_label[i1]-1])):
                                ls.append(param_group[Group_label[i1]-1][i2]/sum(param_group[Group_label[i1]-1]))

                            ls1 = [0]
                            for i2 in range(len(ls)):
                                ls1.append(sum(ls[:i2+1]))

                            p1 = random.random()
                            for i2 in range(len(ls1)-1):
                                if p1 > ls1[i2] and p1 <= ls1[i2+1]:
                                    break
                            Group_label[i1] = [Group_label[i1],i2+1]
                dictx1 = {}
                i1 = 0
                for i in range(len(param_group)):
                    if type(param_group[i]) == type([]):
                        dictx1[i+1] = i1
                        i1 += 1
                ls1b = []
                for i in range(len(param_group)):
                    if type(param_group[i]) == type([]):
                        ls1b.append(i)

                for i in range(len(ls4)):
                    if DEgene_label[i] != 0:
                        ls1 = []
                        dict1 = {}
                        for i1 in range(len(param_group)):
                            if type(param_group[i1]) == type([]):
                                ls = []
                                for i2 in range(len(param_group[i1])):
                                    k1 = np.random.normal(param_DEgene_var[0],0.8*param_DEgene_var[1],1)
                                    ls.append(float(np.exp(k1)))
                                ls1.append(ls)
                                dict1[i1] = ls
                        for i1 in range(len(ls7_1)):
                            if type(Group_label[i1]) == type([]):
                                start = l1a[i][Group_label[i1][0]-1]
                                end = l1a[i][Group_label[i1][0]-1]*dict1[Group_label[i1][0]-1][Group_label[i1][1]-1]
                                DEgene_factor[i][i1] = float(start+(end-start)*Path[i1])
                    if Noline_gene[i] != 0:
                        KA = []
                        KB = []
                        KC = []
                        for i1 in range(len(param_group)):
                            if type(param_group[i1]) == type([]):
                                ls = []
                                for i2 in range(len(param_group[i1])):
                                    k1 = np.random.normal(param_DEgene_var[0],0.8*param_DEgene_var[1],1)
                                    ls.append(float(np.exp(k1)))
                                K1 = []
                                K2 = []
                                K3 = []
                                for i2 in range(len(ls)):
                                    t = []
                                    w = []
                                    t.append(0)
                                    w.append(0)
                                    for i3 in range(50):
                                        t.append((i3+1)*0.02)
                                        k1 = np.random.normal(0,0.02)
                                        w.append(w[i3]+k1)
                                    b = []
                                    for i3 in range(len(t)):
                                        b.append(w[i3]+l2a[i][i1]-t[i3]*(w[-1])+t[i3]*(ls[i2]*l2a[i][i1]-l2a[i][i1]))
                                    p2 = ls[i2]*l2a[i][i1]
                                    def funa(x,k1,k2,k3):
                                        return k1*(x-x**2)+k2*(x**3-x**4)+k3*(x**5-x**6)+l2a[i][i1]+(p2-l2a[i][i1])*x

                                    x = np.array(t)
                                    y = np.array(b)
                                    popt, pcov = curve_fit(funa,x,y)
                                    k1 = popt[0]
                                    k2 = popt[1]
                                    k3 = popt[2]
                                    K1.append(k1)
                                    K2.append(k2)
                                    K3.append(k3)
                                KA.append(K1)
                                KB.append(K2)
                                KC.append(K3)
                        for i1 in range(len(ls7_1)):
                            if type(Group_label[i1]) == type([]):
                                x1 = Path[i1]
                                if KA[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1-x1**2)+KB[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**3-x1**4)+KC[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**5-x1**6)+l2a[i][Group_label[i1][0]-1]+(ls[Group_label[i1][1]-1]*l2a[i][Group_label[i1][0]-1]-l2a[i][Group_label[i1][0]-1])*x1 > 0:
                                    DEgene_factor[i][i1] = KA[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1-x1**2)+KB[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**3-x1**4)+KC[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**5-x1**6)+l2a[i][Group_label[i1][0]-1]+(ls[Group_label[i1][1]-1]*l2a[i][Group_label[i1][0]-1]-l2a[i][Group_label[i1][0]-1])*x1
                                else:
                                    DEgene_factor[i][i1] = 0


                for i1 in range(len(ls7_1)):
                    if type(Group_label[i1]) == type([]):
                        Path[i1] = 0.5+Path[i1]/2
                    else:
                        if Group_label[i1]-1 in ls1b:
                            ls = []
                            for i2 in range(len(param_group[Group_label[i1]-1])):
                                ls.append(param_group[Group_label[i1]-1][i2]/sum(param_group[Group_label[i1]-1]))

                            ls1 = [0]
                            for i2 in range(len(ls)):
                                ls1.append(sum(ls[:i2+1]))

                            p1 = random.random()
                            for i2 in range(len(ls1)-1):
                                if p1 > ls1[i2] and p1 <= ls1[i2+1]:
                                    break
                            Group_label[i1] = [Group_label[i1],i2+1]
                            Path[i1] = Path[i1]/2





    lm = []
    LJ = 0
    for i in range(len(ls4p)):
        LJ += ls4p[i][0]

    lp = []
    lv = []
    s = np.array(ls7)
    i1 = 0
    i2a = (len(ls4)+len(ls2_0))//100
    for i in range(len(ls4)+len(ls2_0)):
        if i%i2a == 0:
            progress(int(i/i2a),width=50)
            #print('{}%'.format(int(i/i2)))
        l = []
        if i in ls2_0:
            for j in range(len(ls7_1)):
                l.append(0)
            lp.append(0)
            lv.append(0)
        else:
            if param_batch == [1] and param_group == [1] and param_path == 'No' and param_library == 1 and cell_number == len(ls7):

                if i1 in l_c:
                    l1 = []
                    l2 = []
                    l4 = []
                    l5 = []
                    for j in range(len(ls4p[0])):
                        l1.append(1/ls12c[i1])
                        l2.append(1/(ls12c[i1]*ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]+1))
                    l4 = k[dictx[i1]]
                    l5 = list(norm.cdf(l4))
                    l = list(nbinom.ppf(l5, l1, l2))





                else:
                    ls1 = np.random.gamma(1/ls12c[i1],ls4[i1]*ls12c[i1],len(ls4p[0]))


                    l1 = []
                    for j in range(len(ls1)):
                        l1.append(float(ls1[j]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]))
                    kkk = 0
                    while 1:
                        kkk+=1
                        x = np.random.poisson(tuple(l1),(1,len(ls4p[0])))
                        ls1 = list(map(int,x[0]))
                        if max(ls1) > 0 or kkk > 20:
                            break
                    l = ls1







                mean1 = np.mean(np.log2(np.array(ls4p[i1])*1000000/LJ+1))
                p_1 = ls4p[i1].count(0)/len(ls4p[i1])
                v_1 = np.var(np.log2(np.array(ls4p[i1])*1000000/LJ+1))
                mean2 = np.mean(np.log2(np.array(l)/s*1000000/LJ+1))
                L = abs(mean2-mean1)

                l7 = []
                l8 = list(np.random.normal(0,ls12c[i1]/5,10))
                l7.extend(l8)
                l8 = list(np.random.normal(0,ls12c[i1]/30,20))
                l7.extend(l8)
                for i2 in range(len(l7)):
                    if ls12c[i1] + l7[i2] > 0:
                        if i1 in l_c:
                            l1 = []
                            l2 = []
                            l4 = []
                            l5 = []
                            for j in range(len(ls4p[0])):
                                l1.append(1/(ls12c[i1]+l7[i2]))
                                l2.append(1/((ls12c[i1]+l7[i2])*ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]+1))
                            l4 = k[dictx[i1]]
                            l5 = list(norm.cdf(l4))
                            la = list(nbinom.ppf(l5, l1, l2))





                        else:
                            ls1 = np.random.gamma(1/(ls12c[i1]+l7[i2]),ls4[i1]*(ls12c[i1]+l7[i2]),len(ls4p[0]))


                            l1 = []
                            for j in range(len(ls1)):
                                l1.append(float(ls1[j]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]))
                            kkk = 0
                            while 1:
                                kkk+=1
                                x = np.random.poisson(tuple(l1),(1,len(ls4p[0])))
                                ls1 = list(map(int,x[0]))
                                if max(ls1) > 0 or kkk > 20:
                                    break
                            la = ls1
                        mean2 = np.mean(np.log2(np.array(la)/s*1000000/LJ+1))
                        if abs(mean2-mean1) < L:
                            L = abs(mean2-mean1)
                            l = la
                            ls12c[i1] = ls12c[i1] + l7[i2]

                if i1 not in l_c:
                    for i2 in range(10):
                        ls1 = np.random.gamma(1/(ls12c[i1]),ls4[i1]*(ls12c[i1]),len(ls4p[0]))


                        l1 = []
                        for j in range(len(ls1)):
                            l1.append(float(ls1[j]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]))
                        kkk = 0
                        while 1:
                            kkk+=1
                            x = np.random.poisson(tuple(l1),(1,len(ls4p[0])))
                            ls1 = list(map(int,x[0]))
                            if max(ls1) > 0 or kkk > 20:
                                break
                        la = ls1
                        mean2 = np.mean(np.log2(np.array(la)/s*1000000/LJ+1))
                        if abs(mean2-mean1) < L:
                            L = abs(mean2-mean1)
                            l = la
                lp.append(abs(l.count(0)/len(l)-p_1))
                lv.append(abs(np.var(np.log2(np.array(l)/s*1000000/LJ+1))-v_1))

                if i1 in l_c:
                    l6 = []
                    for i3 in range(len(l)):
                        l6.append(l[i3]/(ls7[i3]*batch_factor[i3]))
                    if max(l6) > 0:
                        lm.append(l6)
                i1 += 1

            else:
                if i1 in l_c:
                    l1 = []
                    l2 = []
                    l4 = []
                    l5 = []
                    for j in range(len(ls7_1)):
                        l1.append(1/ls12c[i1])
                        if param_batch == [1]:
                            l2.append(1/(ls12c[i1]*ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[j]+1))
                        else:
                            l2.append(1/(ls12c[i1]*ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[i1][j]+1))
                    l4 = k[dictx[i1]]
                    l5 = list(norm.cdf(l4))
                    try:
                        l = list(nbinom.ppf(l5, l1, l2))
                    except:
                        print('erro')

                    l6 = []
                    for i3 in range(len(l)):
                        l6.append(l[i3]/(ls7_1[i3]))
                    if max(l6) > 0:
                        lm.append(l6)



                else:
                    ls1 = np.random.gamma(1/ls12c[i1],ls4[i1]*ls12c[i1],len(ls7_1))


                    l1 = []
                    for j in range(len(ls1)):
                        if param_batch == [1]:
                            l1.append(float(ls1[j]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[j]))
                        else:
                            l1.append(float(ls1[j]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[i1][j]))
                    kkk = 0
                    while 1:
                        kkk+=1
                        x = np.random.poisson(tuple(l1),(1,len(ls7_1)))
                        ls1 = list(map(int,x[0]))
                        if max(ls1) > 0 or kkk > 20:
                            break
                    l = ls1
                i1 += 1
        lcount.append(l)



    if len(lp) > 0:
        lp1 = np.argsort(lp)
        lv1 = np.argsort(lv)

    i1 = 0
    for i in range(len(ls4)+len(ls2_0)):
        l = []
        if i not in ls2_0:
            if param_batch == [1] and param_group == [1] and param_path == 'No' and param_library == 1 and cell_number == len(ls7):
                if i1 in l_c:
                    if i in lp1[-int(len(lp1)*0.08):] or i in lv1[-int(len(lv1)*0.04):]:
                        l7 = []
                        l8 = list(np.random.normal(0,ls12c[i1]/2,5))
                        l7.extend(l8)
                        l8 = list(np.random.normal(0,ls12c[i1]/5,10))
                        l7.extend(l8)
                        l8 = list(np.random.normal(0,ls12c[i1]/30,20))
                        l7.extend(l8)
                        if i in lv1[-int(len(lv1)*0.04):]:
                            v_1 = np.var(np.log2(np.array(ls4p[i1])*1000000/LJ+1))
                            v_2 = np.var(np.log2(np.array(lcount[i])/s*1000000/LJ+1))
                            L = abs(v_2-v_1)
                            for i2 in range(len(l7)):
                                if ls12c[i1] + l7[i2] > 0:
                                    l1 = []
                                    l2 = []
                                    l4 = []
                                    l5 = []
                                    for j in range(len(ls4p[0])):
                                        l1.append(1/(ls12c[i1]+l7[i2]))
                                        l2.append(1/((ls12c[i1]+l7[i2])*ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]+1))
                                    l4 = k[dictx[i1]]
                                    l5 = list(norm.cdf(l4))
                                    la = list(nbinom.ppf(l5, l1, l2))
                                    if abs(np.var(np.log2(np.array(la)/s*1000000/LJ+1))-v_1) < L:
                                        L = abs(np.var(np.log2(np.array(la)/s*1000000/LJ+1))-v_1)
                                        ls12c[i1] = ls12c[i1] + l7[i2]
                                        l = la
                                        lcount[i] = l


                        else:
                            p_1 = ls4p[i1].count(0)/len(ls4p[i1])
                            p_2 = lcount[i].count(0)/len(lcount[i])
                            L = abs(p_2-p_1)
                            for i2 in range(len(l7)):
                                if ls12c[i1] + l7[i2] > 0:
                                    l1 = []
                                    l2 = []
                                    l4 = []
                                    l5 = []
                                    for j in range(len(ls4p[0])):
                                        l1.append(1/(ls12c[i1]+l7[i2]))
                                        l2.append(1/((ls12c[i1]+l7[i2])*ls4[i1]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]+1))
                                    l4 = k[dictx[i1]]
                                    l5 = list(norm.cdf(l4))
                                    la = list(nbinom.ppf(l5, l1, l2))
                                    if abs(la.count(0)/len(la)-p_1) < L:
                                        L = abs(la.count(0)/len(la)-p_1)
                                        ls12c[i1] = ls12c[i1] + l7[i2]
                                        l = la
                                        lcount[i] = l

                else:
                    if i in lp1[-int(len(lp1)*0.08):] or i in lv1[-int(len(lv1)*0.04):]:
                        l7 = []
                        l8 = list(np.random.normal(0,ls12c[i1]/2,5))
                        l7.extend(l8)
                        l8 = list(np.random.normal(0,ls12c[i1]/5,10))
                        l7.extend(l8)
                        l8 = list(np.random.normal(0,ls12c[i1]/30,20))
                        l7.extend(l8)
                        if i in lv1[-int(len(lv1)*0.04):]:
                            v_1 = np.var(np.log2(np.array(ls4p[i1])*1000000/LJ+1))
                            v_2 = np.var(np.log2(np.array(lcount[i])/s*1000000/LJ+1))
                            L = abs(v_2-v_1)
                            for i2 in range(len(l7)):
                                if ls12c[i1] + l7[i2] > 0:
                                    ls1 = np.random.gamma(1/(ls12c[i1]+l7[i2]),ls4[i1]*(ls12c[i1]+l7[i2]),len(ls4p[0]))
                                    l1 = []
                                    for j in range(len(ls1)):
                                        l1.append(float(ls1[j]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]))
                                    kkk = 0
                                    while 1:
                                        kkk+=1
                                        x = np.random.poisson(tuple(l1),(1,len(ls4p[0])))
                                        ls1 = list(map(int,x[0]))
                                        if max(ls1) > 0 or kkk > 20:
                                            break
                                    la = ls1
                                    if abs(np.var(np.log2(np.array(la)/s*1000000/LJ+1))-v_1) < L:
                                        L = abs(np.var(np.log2(np.array(la)/s*1000000/LJ+1))-v_1)
                                        ls12c[i1] = ls12c[i1] + l7[i2]
                                        l = la
                                        lcount[i] = l
                            for i2 in range(10):
                                ls1 = np.random.gamma(1/ls12c[i1],ls4[i1]*ls12c[i1],len(ls4p[0]))
                                l1 = []
                                for j in range(len(ls1)):
                                    l1.append(float(ls1[j]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]))
                                kkk = 0
                                while 1:
                                    kkk+=1
                                    x = np.random.poisson(tuple(l1),(1,len(ls4p[0])))
                                    ls1 = list(map(int,x[0]))
                                    if max(ls1) > 0 or kkk > 20:
                                        break
                                la = ls1
                                if abs(np.var(np.log2(np.array(la)/s*1000000/LJ+1))-v_1) < L:
                                    L = abs(np.var(np.log2(np.array(la)/s*1000000/LJ+1))-v_1)
                                    l = la
                                    lcount[i] = l





                        else:
                            p_1 = ls4p[i1].count(0)/len(ls4p[i1])
                            p_2 = lcount[i].count(0)/len(lcount[i])
                            L = abs(p_2-p_1)
                            for i2 in range(len(l7)):
                                if ls12c[i1] + l7[i2] > 0:
                                    ls1 = np.random.gamma(1/(ls12c[i1]+l7[i2]),ls4[i1]*(ls12c[i1]+l7[i2]),len(ls4p[0]))
                                    l1 = []
                                    for j in range(len(ls1)):
                                        l1.append(float(ls1[j]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]))
                                    kkk = 0
                                    while 1:
                                        kkk+=1
                                        x = np.random.poisson(tuple(l1),(1,len(ls4p[0])))
                                        ls1 = list(map(int,x[0]))
                                        if max(ls1) > 0 or kkk > 20:
                                            break
                                    la = ls1
                                    if abs(la.count(0)/len(la)-p_1) < L:
                                        L = abs(la.count(0)/len(la)-p_1)
                                        ls12c[i1] = ls12c[i1] + l7[i2]
                                        l = la
                                        lcount[i] = l
                            for i2 in range(10):
                                ls1 = np.random.gamma(1/ls12c[i1],ls4[i1]*ls12c[i1],len(ls4p[0]))
                                l1 = []
                                for j in range(len(ls1)):
                                    l1.append(float(ls1[j]*ls7[j]*DEgene_factor[i1][j]*batch_factor[j]))
                                kkk = 0
                                while 1:
                                    kkk+=1
                                    x = np.random.poisson(tuple(l1),(1,len(ls4p[0])))
                                    ls1 = list(map(int,x[0]))
                                    if max(ls1) > 0 or kkk > 20:
                                        break
                                la = ls1
                                if abs(la.count(0)/len(la)-p_1) < L:
                                    L = abs(la.count(0)/len(la)-p_1)
                                    l = la
                                    lcount[i] = l
            i1 += 1




    progress(100,width=50)
    print(' ')
    print('Simulation completed')
    #if mod != '-c nocopula':

        #Mat = np.array(lm)
        #correlation1 = np.corrcoef(Mat)


        #sns.set()
        #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        #plt.show()

        #ln = []
        #for i in range(len(lm)):
        #    ln.append(float(np.mean(lm[i])))
        #ln1 = []
        #ln1 = list(reversed(list(np.argsort(ln))))
        #ln2 = []
        #for i in range(len(lm)):
        #    if i in ln1[:50]:
        #        ln2.append(lm[i])
        #Mat = np.array(ln2)
        #correlation1 = np.corrcoef(Mat)

        #sns.set()
        #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        #plt.show()






def creatcount_noparam_test_imputation_cor():
    """

    Simulate count

    Parameters
    ----------
    ls4:list
          Mean of gene
    ls12:list
          Dispersion of gene
    ls7:list
          Size factor


    Returns
    -------
    ldp:list
          Simulated dropouted gene
    ldropout:list
          Dropout information
    llisize:list
          Library size
    lcount:list
          Count

    """
    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')

    from scipy.stats import nbinom
    from scipy.stats import norm
    from scipy.special import gdtr, gdtrix
    global lcount,ls4,correlation,ls4p,ls12c,ls7,dictx,lgamma,lgamma_sj,ls7_1
    lgamma = []
    lgamma_sj = []
    print(' ')
    print('Start simulation...')
    ls7_1 = []
    if cell_number == len(ls7):
        ls7_1 = ls7
    else:
        if cell_number < len(ls7):
            ls7_1 = ls7[:cell_number]
        else:
            ls7_2 = np.log(ls7)
            x_mean, x_std = norm.fit(ls7_2)
            ls7_1 = [i for i in ls7]
            while 1:
                k = float(np.exp(np.random.normal(x_mean, x_std,1)))
                ls7_1.append(k)
                if len(ls7_1) == cell_number:
                    break
    lcount = []
    if mod != '-g nocopula':
        l5 = []
        for i in range(len(correlation)):
            l5.append(0)
        p = np.random.multivariate_normal(mean=l5,cov=correlation,size=len(ls7_1))

        k  = list(np.transpose(p))
        #k1 = np.corrcoef(np.array(k))
        # sns.set()
        # ax = sns.clustermap(k1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        # plt.show()
        #
        # k1 = np.corrcoef(np.array(k[:50]))
        # sns.set()
        # ax = sns.clustermap(k1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
        # plt.show()


    global batch_factor,batch_label
    batch_factor = []
    batch_label = []
    if param_batch == [1] and param_library == 1:
        for i in range(len(ls7_1)):
            batch_factor.append(1)
    else:
        if param_library != 1 and param_batch == [1]:
            for i in range(len(ls7_1)):
                batch_factor.append(param_library)
        else:
            for i in range(len(ls7_1)):
                p1 = random.random()
                ls = [0]
                for j in range(len(param_batch)):
                    ls.append(sum(param_batch[:j+1]))
                for j in range(len(ls)-1):
                    if p1 > ls[j] and p1 < ls[j+1]:
                        break
                batch_label.append(j+1)                 ####batch_label
            for i in range(len(ls4p)):
                ls = []
                for j in range(len(param_batch)):
                    k1 = np.random.normal(param_batch_var[0],param_batch_var[1],1)
                    ls.append(k1)
                ls1 = []
                for j in range(len(ls7_1)):
                    ls1.append(param_library*float(np.exp(ls[batch_label[j]-1])))
                batch_factor.append(ls1)                             #######batch_factor
    global DEgene_factor,Group_label,DEgene_label,Path,Noline_gene,l_marker
    l_marker = []
    DEgene_label = []
    DEgene_factor = []
    Group_label = []
    Path = []
    Noline_gene = []
    if param_group == [1]:
        for i in range(len(ls4)):
            l = [1 for j in range(len(ls7_1))]
            DEgene_factor.append(l)
    else:
        for i in range(len(ls4)):
            p1 = random.random()
            if p1 < param_DEgene*(len(ls4)+len(ls2_0))/len(ls4):
                DEgene_label.append(1)
            else:
                DEgene_label.append(0)              ####DEgene_label

        param_group1 = []
        tree = 0
        for i in range(len(param_group)):
            if type(param_group[i]) == type([]):
                param_group1.append(sum(param_group[i]))
                tree = 1
            else:
                param_group1.append(param_group[i])


        for i in range(len(ls7_1)):
            p1 = random.random()
            ls = [0]
            for j in range(len(param_group1)):
                ls.append(sum(param_group1[:j+1]))
            for j in range(len(ls)-1):
                if p1 > ls[j] and p1 < ls[j+1]:
                    break
            Group_label.append(j+1)                 ####Group_label
        if param_path == 'No':
            for i in range(len(ls4)):
                if DEgene_label[i] == 0:
                    l = [1 for i1 in range(len(ls7_1))]
                else:
                    while 1:
                        ls = []
                        for i1 in range(len(param_group1)):
                            k1 = np.random.normal(param_DEgene_var[0],param_DEgene_var[1],1)
                            ls.append(k1)
                        if min(ls)<param_DEgene_var[0] and max(ls)>param_DEgene_var[0] and (max(ls)-min(ls)>0.2 or max(ls)-min(ls)>param_DEgene_var[1]):
                            break
                    l = []
                    for i1 in range(len(ls7_1)):
                        l.append(float(np.exp(ls[Group_label[i1]-1])))
                    p1 = random.random()
                    if p1 < param_marker/param_DEgene:
                        l_marker.append(i)
                        l2 = l
                        l1 = max(l)
                        for i1 in range(len(l2)):
                            if l2[i1] != l1:
                                l[i1] = 0.0001
                DEgene_factor.append(l)
            if tree == 1:

                for i1 in range(len(ls7_1)):
                    if type(param_group[Group_label[i1]-1]) == type([]):
                        p = random.random()
                        if p < 1.1:
                            ls = []
                            for i2 in range(len(param_group[Group_label[i1]-1])):
                                ls.append(param_group[Group_label[i1]-1][i2]/sum(param_group[Group_label[i1]-1]))

                            ls1 = [0]
                            for i2 in range(len(ls)):
                                ls1.append(sum(ls[:i2+1]))

                            p1 = random.random()
                            for i2 in range(len(ls1)-1):
                                if p1 > ls1[i2] and p1 <= ls1[i2+1]:
                                    break
                            Group_label[i1] = [Group_label[i1],i2+1]



                for i in range(len(ls4)):
                    if DEgene_label[i] != 0:
                        ls1 = []
                        dict1 = {}
                        for i1 in range(len(param_group)):
                            if type(param_group[i1]) == type([]):
                                ls = []
                                for i2 in range(len(param_group[i1])):
                                    k1 = np.random.normal(param_DEgene_var[0],0.8*param_DEgene_var[1],1)
                                    ls.append(float(np.exp(k1)))
                                ls1.append(ls)
                                dict1[i1] = ls

                        for i1 in range(len(ls7_1)):
                            if type(Group_label[i1]) == type([]):
                                DEgene_factor[i][i1] = DEgene_factor[i][i1]*dict1[Group_label[i1][0]-1][Group_label[i1][1]-1]


        else:
            l1a = []
            l2a = []
            for i in range(len(ls7_1)):
                p1 = random.random()
                Path.append(float(p1))
            ls = []
            while 1:
                p1 = random.randint(0,len(ls4)-1)
                if DEgene_label[p1] == 0 and p1 not in ls:
                    ls.append(p1)
                if len(ls) >= param_noline*(len(ls4)+len(ls2_0)):
                    break
            for i in range(len(ls4)):
                if i in ls:
                    Noline_gene.append(1)
                else:
                    Noline_gene.append(0)
            for i in range(len(ls4)):
                if DEgene_label[i] == 0 and Noline_gene[i] == 0:
                    l = [1 for i1 in range(len(ls7_1))]
                    l1a.append(0)
                    l2a.append(0)
                else:
                    l = []
                    if DEgene_label[i] == 1:
                        while 1:
                            ls = []
                            for i1 in range(len(param_group)):
                                k1 = np.random.normal(param_DEgene_var[0],param_DEgene_var[1],1)
                                ls.append(k1)
                            if min(ls)<param_DEgene_var[0] and max(ls)>param_DEgene_var[0] and (max(ls)-min(ls)>0.2 or max(ls)-min(ls)>param_DEgene_var[1]):
                                break
                        for i1 in range(len(ls7_1)):
                            l.append((float(np.exp(ls[Group_label[i1]-1]))-1)*Path[i1]+1)
                        l1a.append(list(np.exp(ls)))
                        l2a.append(0)
                    else:
                        ls = []
                        for i in range(len(param_group)):
                            k1 = np.random.normal(param_DEgene_var[0],0.5*param_DEgene_var[1],1)
                            ls.append(float(np.exp(k1)))
                        K1 = []
                        K2 = []
                        K3 = []
                        for i1 in range(len(ls)):
                            t = []
                            w = []
                            t.append(0)
                            w.append(0)
                            for i in range(50):
                                t.append((i+1)*0.02)
                                k1 = np.random.normal(0,0.01)
                                w.append(w[i]+k1)
                            b = []
                            for i in range(len(t)):
                                b.append(w[i]+1-t[i]*(w[-1])+t[i]*(ls[i1]-1))
                            p2 = ls[i1]
                            def funa(x,k1,k2,k3):
                                return k1*(x-x**2)+k2*(x**3-x**4)+k3*(x**5-x**6)+1+(p2-1)*x

                            x = np.array(t)
                            y = np.array(b)
                            popt, pcov = curve_fit(funa,x,y)
                            k1 = popt[0]
                            k2 = popt[1]
                            k3 = popt[2]
                            K1.append(k1)
                            K2.append(k2)
                            K3.append(k3)
                        for i1 in range(len(ls7_1)):
                            x1 = Path[i1]
                            if K1[Group_label[i1]-1]*(x1-x1**2)+K2[Group_label[i1]-1]*(x1**3-x1**4)+K3[Group_label[i1]-1]*(x1**5-x1**6)+1+(ls[Group_label[i1]-1]-1)*x1 > 0:
                                l.append(K1[Group_label[i1]-1]*(x1-x1**2)+K2[Group_label[i1]-1]*(x1**3-x1**4)+K3[Group_label[i1]-1]*(x1**5-x1**6)+1+(ls[Group_label[i1]-1]-1)*x1)
                            else:
                                l.append(0)
                        l1a.append(0)
                        l2a.append(list(ls))
                DEgene_factor.append(l)


            if tree == 1:
                for i1 in range(len(ls7_1)):
                    if type(param_group[Group_label[i1]-1]) == type([]):
                        p = random.random()
                        if p < 0.5:
                            ls = []
                            for i2 in range(len(param_group[Group_label[i1]-1])):
                                ls.append(param_group[Group_label[i1]-1][i2]/sum(param_group[Group_label[i1]-1]))

                            ls1 = [0]
                            for i2 in range(len(ls)):
                                ls1.append(sum(ls[:i2+1]))

                            p1 = random.random()
                            for i2 in range(len(ls1)-1):
                                if p1 > ls1[i2] and p1 <= ls1[i2+1]:
                                    break
                            Group_label[i1] = [Group_label[i1],i2+1]

                dictx1 = {}
                i1 = 0
                for i in range(len(param_group)):
                    if type(param_group[i]) == type([]):
                        dictx1[i+1] = i1
                        i1 += 1

                ls1b = []
                for i in range(len(param_group)):
                    if type(param_group[i]) == type([]):
                        ls1b.append(i)

                for i in range(len(ls4)):
                    if DEgene_label[i] != 0:
                        ls1 = []
                        dict1 = {}
                        for i1 in range(len(param_group)):
                            if type(param_group[i1]) == type([]):
                                ls = []
                                for i2 in range(len(param_group[i1])):
                                    k1 = np.random.normal(param_DEgene_var[0],0.8*param_DEgene_var[1],1)
                                    ls.append(float(np.exp(k1)))
                                ls1.append(ls)
                                dict1[i1] = ls
                        for i1 in range(len(ls7_1)):
                            if type(Group_label[i1]) == type([]):
                                start = l1a[i][Group_label[i1][0]-1]
                                end = l1a[i][Group_label[i1][0]-1]*dict1[Group_label[i1][0]-1][Group_label[i1][1]-1]
                                DEgene_factor[i][i1] = float(start+(end-start)*Path[i1])
                    if Noline_gene[i] != 0:
                        KA = []
                        KB = []
                        KC = []
                        for i1 in range(len(param_group)):
                            if type(param_group[i1]) == type([]):
                                ls = []
                                for i2 in range(len(param_group[i1])):
                                    k1 = np.random.normal(param_DEgene_var[0],0.8*param_DEgene_var[1],1)
                                    ls.append(float(np.exp(k1)))
                                K1 = []
                                K2 = []
                                K3 = []
                                for i2 in range(len(ls)):
                                    t = []
                                    w = []
                                    t.append(0)
                                    w.append(0)
                                    for i3 in range(50):
                                        t.append((i3+1)*0.02)
                                        k1 = np.random.normal(0,0.02)
                                        w.append(w[i3]+k1)
                                    b = []
                                    for i3 in range(len(t)):
                                        b.append(w[i3]+l2a[i][i1]-t[i3]*(w[-1])+t[i3]*(ls[i2]*l2a[i][i1]-l2a[i][i1]))
                                    p2 = ls[i2]*l2a[i][i1]
                                    def funa(x,k1,k2,k3):
                                        return k1*(x-x**2)+k2*(x**3-x**4)+k3*(x**5-x**6)+l2a[i][i1]+(p2-l2a[i][i1])*x

                                    x = np.array(t)
                                    y = np.array(b)
                                    popt, pcov = curve_fit(funa,x,y)
                                    k1 = popt[0]
                                    k2 = popt[1]
                                    k3 = popt[2]
                                    K1.append(k1)
                                    K2.append(k2)
                                    K3.append(k3)
                                KA.append(K1)
                                KB.append(K2)
                                KC.append(K3)
                        for i1 in range(len(ls7_1)):
                            if type(Group_label[i1]) == type([]):
                                x1 = Path[i1]
                                if KA[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1-x1**2)+KB[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**3-x1**4)+KC[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**5-x1**6)+l2a[i][Group_label[i1][0]-1]+(ls[Group_label[i1][1]-1]*l2a[i][Group_label[i1][0]-1]-l2a[i][Group_label[i1][0]-1])*x1 > 0:
                                    DEgene_factor[i][i1] = KA[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1-x1**2)+KB[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**3-x1**4)+KC[dictx1[Group_label[i1][0]]][Group_label[i1][1]-1]*(x1**5-x1**6)+l2a[i][Group_label[i1][0]-1]+(ls[Group_label[i1][1]-1]*l2a[i][Group_label[i1][0]-1]-l2a[i][Group_label[i1][0]-1])*x1
                                else:
                                    DEgene_factor[i][i1] = 0



                for i1 in range(len(ls7_1)):
                    if type(Group_label[i1]) == type([]):
                        Path[i1] = 0.5+Path[i1]/2
                    else:
                        if Group_label[i1]-1 in ls1b:
                            ls = []
                            for i2 in range(len(param_group[Group_label[i1]-1])):
                                ls.append(param_group[Group_label[i1]-1][i2]/sum(param_group[Group_label[i1]-1]))

                            ls1 = [0]
                            for i2 in range(len(ls)):
                                ls1.append(sum(ls[:i2+1]))

                            p1 = random.random()
                            for i2 in range(len(ls1)-1):
                                if p1 > ls1[i2] and p1 <= ls1[i2+1]:
                                    break
                            Group_label[i1] = [Group_label[i1],i2+1]
                            Path[i1] = Path[i1]/2






    lm = []

    i1 = 0
    i2 = (len(ls4)+len(ls2_0))//100
    for i in range(len(ls4)+len(ls2_0)):
        if i%i2 == 0:
            progress(int(i/i2),width=50)
            #print('{}%'.format(int(i/i2)))
        l = []
        l_1 = []
        l_2 = []
        if i in ls2_0:
            for j in range(len(ls7_1)):
                l_1.append(0)
                l_2.append(0)
                l.append(0)

        else:
            if i1 in l_c:
                l1 = []
                l2 = []
                l3 = []
                l4 = []
                l5 = []
                for j in range(len(ls7_1)):
                    l1.append(1/ls12c[i1])
                    if param_batch == [1]:
                        l2.append(1/(ls12c[i1]*ls4[i1]*DEgene_factor[i1][j]*batch_factor[j]))
                        l3.append(1/(ls12c[i1]*ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[j]))
                    else:
                        l2.append(1/(ls12c[i1]*ls4[i1]*DEgene_factor[i1][j]*batch_factor[i1][j]))
                        l3.append(1/(ls12c[i1]*ls4[i1]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[i1][j]))
                l4 = k[dictx[i1]]
                l5 = list(norm.cdf(l4))
                l_1 = list(gdtrix(l2, l1, l5))
                l_2 = list(gdtrix(l3, l1, l5))
                l_1 = np.nan_to_num(l_1)
                l_2 = np.nan_to_num(l_2)




                kkk = 0
                while 1:
                    kkk+=1
                    x = np.random.poisson(tuple(l_2),(1,len(ls7_1)))
                    ls1 = list(map(int,x[0]))
                    if max(ls1) > 0 or kkk > 20:
                        break
                l = ls1

                l6 = []
                for i3 in range(len(l)):
                    l6.append(l[i3]/(ls7_1[i3]))
                if max(l6) > 0:
                    lm.append(l6)

            else:
                la = list(np.random.gamma(1/ls12c[i1],ls4[i1]*ls12c[i1],len(ls7_1)))


                for j in range(len(la)):
                    if param_batch == [1]:
                        l_1.append(float(la[j]*DEgene_factor[i1][j]*batch_factor[j]))
                        l_2.append(float(la[j]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[j]))
                    else:
                        l_1.append(float(la[j]*DEgene_factor[i1][j]*batch_factor[i1][j]))
                        l_2.append(float(la[j]*ls7_1[j]*DEgene_factor[i1][j]*batch_factor[i1][j]))
                kkk = 0
                while 1:
                    kkk+=1
                    x = np.random.poisson(tuple(l_2),(1,len(ls7_1)))
                    ls1 = list(map(int,x[0]))
                    if max(ls1) > 0 or kkk > 20:
                        break
                l = ls1

            i1 += 1
        lgamma.append(l_1)
        lgamma_sj.append(l_2)
        lcount.append(l)
    progress(100,width=50)
    print(' ')
    print('Simulation completed')
    #if mod != '-g nocopula':

    #    Mat = np.array(lm)
    #    correlation1 = np.corrcoef(Mat)


    #    sns.set()
    #    ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
    #    plt.show()

    #    ln = []
    #    for i in range(len(lm)):
    #        ln.append(float(np.mean(lm[i])))
    #    ln1 = []
    #    ln1 = list(reversed(list(np.argsort(ln))))
    #    ln2 = []
    #    for i in range(len(lm)):
    #        if i in ln1[:50]:
    #            ln2.append(lm[i])
    #    Mat = np.array(ln2)
    #    correlation1 = np.corrcoef(Mat)

    #    sns.set()
    #    ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
    #    plt.show()




def write_noparam():
    """
    Write file

    Parameters
    ----------
    dp:list
          Simulated dropouted gene
    ldropout:list
          Dropout information
    llisize:list
          Library size
    lcount:list
          Count

    """
    print('Start writing...')
    global ls2,ls3,ls7,ls4p,ls4,ls2_0,ls12c,l_c,correlation,lcount,ltruecount,ld,batch_factor,batch_label,DEgene_factor,Group_label,DEgene_label,Path,Noline_gene

    f = open(path3+'/Simulation information.txt','w')
    if mod == '-c':
        f.write('Method:SimCH-copula-NB')
        f.write('\n')
        f.write('Data:'+path1)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('Copula gene number'+str(copula_number))
        f.write('\n')
        f.write('-1\tGroup ratio:'+str(param_group))
        f.write('\n')
        f.write('-2\tBatch ratio:'+str(param_batch))
        f.write('\n')
        f.write('-3\tBatch variation:'+str(param_batch_var))
        f.write('\n')
        f.write('-4\tPath:'+str(param_path))
        f.write('\n')
        f.write('-5\tDEgene ratio:'+str(param_DEgene))
        f.write('\n')
        f.write('-6\tDEgene variation:'+str(param_DEgene_var))
        f.write('\n')
        f.write('-7\tNon-linear gene ratio:'+str(param_noline))
        f.write('\n')
        f.write('-8\tMarker gene ratio:'+str(param_marker))
        f.write('\n')
        f.write('-9\tLibrary magnification:'+str(param_library))
        f.write('\n')
        f.write('-10\tCell number:'+str(cell_number))
        f.write('\n')
        f.close()
    elif mod == '-c nocopula':
        f.write('Method:SimCH-fit-NB')
        f.write('\n')
        f.write('Data:'+path1)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('-1\tGroup ratio:'+str(param_group))
        f.write('\n')
        f.write('-2\tBatch ratio:'+str(param_batch))
        f.write('\n')
        f.write('-3\tBatch variation:'+str(param_batch_var))
        f.write('\n')
        f.write('-4\tPath:'+str(param_path))
        f.write('\n')
        f.write('-5\tDEgene ratio:'+str(param_DEgene))
        f.write('\n')
        f.write('-6\tDEgene variation:'+str(param_DEgene_var))
        f.write('\n')
        f.write('-7\tNon-linear gene ratio:'+str(param_noline))
        f.write('\n')
        f.write('-8\tMarker gene ratio:'+str(param_marker))
        f.write('\n')
        f.write('-9\tLibrary magnification:'+str(param_library))
        f.write('\n')
        f.write('-10\tCell number:'+str(cell_number))
        f.write('\n')
        f.close()
    elif mod == '-d':
        f.write('Method:SimCH-copula-NBZI')
        f.write('\n')
        f.write('Data:'+path1)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('Copula gene number'+str(copula_number))
        f.write('\n')
        f.write('-1\tGroup ratio:'+str(param_group))
        f.write('\n')
        f.write('-2\tBatch ratio:'+str(param_batch))
        f.write('\n')
        f.write('-3\tBatch variation:'+str(param_batch_var))
        f.write('\n')
        f.write('-4\tPath:'+str(param_path))
        f.write('\n')
        f.write('-5\tDEgene ratio:'+str(param_DEgene))
        f.write('\n')
        f.write('-6\tDEgene variation:'+str(param_DEgene_var))
        f.write('\n')
        f.write('-7\tNon-linear gene ratio:'+str(param_noline))
        f.write('\n')
        f.write('-8\tMarker gene ratio:'+str(param_marker))
        f.write('\n')
        f.write('-9\tLibrary magnification:'+str(param_library))
        f.write('\n')
        f.write('-10\tCell number:'+str(cell_number))
        f.write('\n')
        f.close()
    elif mod == '-d nocopula':
        f.write('Method:SimCH-fit-NBZI')
        f.write('\n')
        f.write('Data:'+path1)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('-1\tGroup ratio:'+str(param_group))
        f.write('\n')
        f.write('-2\tBatch ratio:'+str(param_batch))
        f.write('\n')
        f.write('-3\tBatch variation:'+str(param_batch_var))
        f.write('\n')
        f.write('-4\tPath:'+str(param_path))
        f.write('\n')
        f.write('-5\tDEgene ratio:'+str(param_DEgene))
        f.write('\n')
        f.write('-6\tDEgene variation:'+str(param_DEgene_var))
        f.write('\n')
        f.write('-7\tNon-linear gene ratio:'+str(param_noline))
        f.write('\n')
        f.write('-8\tMarker gene ratio:'+str(param_marker))
        f.write('\n')
        f.write('-9\tLibrary magnification:'+str(param_library))
        f.write('\n')
        f.write('-10\tCell number:'+str(cell_number))
        f.write('\n')
        f.close()
    elif mod == '-g':
        f.write('Method:SimCH-copula-NB test imputation')
        f.write('\n')
        f.write('Data:'+path1)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('Copula gene number'+str(copula_number))
        f.write('\n')
        f.write('-1\tGroup ratio:'+str(param_group))
        f.write('\n')
        f.write('-2\tBatch ratio:'+str(param_batch))
        f.write('\n')
        f.write('-3\tBatch variation:'+str(param_batch_var))
        f.write('\n')
        f.write('-4\tPath:'+str(param_path))
        f.write('\n')
        f.write('-5\tDEgene ratio:'+str(param_DEgene))
        f.write('\n')
        f.write('-6\tDEgene variation:'+str(param_DEgene_var))
        f.write('\n')
        f.write('-7\tNon-linear gene ratio:'+str(param_noline))
        f.write('\n')
        f.write('-8\tMarker gene ratio:'+str(param_marker))
        f.write('\n')
        f.write('-9\tLibrary magnification:'+str(param_library))
        f.write('\n')
        f.write('-10\tCell number:'+str(cell_number))
        f.write('\n')
        f.close()
    elif mod == '-g nocopula':
        f.write('Method:SimCH-fit-NB test imputation')
        f.write('\n')
        f.write('Data:'+path1)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('-1\tGroup ratio:'+str(param_group))
        f.write('\n')
        f.write('-2\tBatch ratio:'+str(param_batch))
        f.write('\n')
        f.write('-3\tBatch variation:'+str(param_batch_var))
        f.write('\n')
        f.write('-4\tPath:'+str(param_path))
        f.write('\n')
        f.write('-5\tDEgene ratio:'+str(param_DEgene))
        f.write('\n')
        f.write('-6\tDEgene variation:'+str(param_DEgene_var))
        f.write('\n')
        f.write('-7\tNon-linear gene ratio:'+str(param_noline))
        f.write('\n')
        f.write('-8\tMarker gene ratio:'+str(param_marker))
        f.write('\n')
        f.write('-9\tLibrary magnification:'+str(param_library))
        f.write('\n')
        f.write('-10\tCell number:'+str(cell_number))
        f.write('\n')
        f.close()




    f = open(path1,'r')
    if mod[-3:] == 'csv':
        reader = csv.reader(f)
        ls = list(reader)
    else:

        f1 = f.readlines()
        ls = []
        for i in f1:
            k = i.strip().split('\t')
            ls.append(k)


    if cell_number == len(ls7):
        f = open(path3+'/Counts.txt','w')
        for i in range(len(ls[0])):
            f.write('\t')
            f.write(str(ls[0][i]))
        f.write('\n')
        for i in range(len(ls4)+len(ls2_0)):
            f.write(ls[i+1][0])
            for j in range(len(ls4p[0])):
                f.write('\t')
                f.write(str(int(lcount[i][j])))
            f.write('\n')
        f.close()
    else:
        f = open(path3+'/Counts.txt','w')
        for i in range(len(lcount[0])):
            f.write('\t')
            f.write('cell'+str(i+1))
        f.write('\n')
        for i in range(len(ls4)+len(ls2_0)):
            f.write(ls[i+1][0])
            for j in range(len(lcount[0])):
                f.write('\t')
                f.write(str(int(lcount[i][j])))
            f.write('\n')
        f.close()

    f = open(path3+'/Gene mean.txt','w')
    i1 = 0
    for i in range(len(ls4)+len(ls2_0)):
        f.write(ls[i+1][0])
        f.write('\t')
        if i in ls2_0:
            f.write(str(0))
        else:
            f.write(str(ls4[i1]))
            i1 += 1
        f.write('\n')
    f.close()
    
    f = open(path3+'/Gene dispersion.txt','w')
    i1 = 0
    for i in range(len(ls4)+len(ls2_0)):
        f.write(ls[i+1][0])
        f.write('\t')
        if i in ls2_0:
            f.write(str(0))
        else:
            f.write(str(ls12c[i1]))
            i1 += 1
        f.write('\n')
    f.close()

    if mod == '-g' or mod == '-g nocopula':
        f = open(path3+'/Cell mean.txt','w')
        if cell_number == len(ls7):
            for i in range(len(ls4p[0])):
                f.write('\t')
                f.write(str(ls[0][i]))
            f.write('\n')
        else:
            for i in range(len(lcount[0])):
                f.write('\t')
                f.write('cell'+str(i+1))
            f.write('\n')
        for i in range(len(ls4)+len(ls2_0)):
            f.write(ls[i+1][0])
            for j in range(len(lcount[0])):
                f.write('\t')
                f.write(str(float(lgamma[i][j])))
            f.write('\n')
        f.close()

        f = open(path3+'/Cell mean size factor.txt','w')
        if cell_number == len(ls7):
            for i in range(len(ls4p[0])):
                f.write('\t')
                f.write(str(ls[0][i]))
            f.write('\n')
        else:
            for i in range(len(lcount[0])):
                f.write('\t')
                f.write('cell'+str(i+1))
            f.write('\n')
        for i in range(len(ls4)+len(ls2_0)):
            f.write(ls[i+1][0])
            for j in range(len(lcount[0])):
                f.write('\t')
                f.write(str(float(lgamma_sj[i][j])))
            f.write('\n')
        f.close()





    if mod == '-d' or mod == '-d nocopula':
        f = open(path3+'/p0.txt','w')
        i1 = 0
        for i in range(len(ls4)+len(ls2_0)):
            f.write(ls[i+1][0])
            f.write('\t')
            if i in ls2_0:
                f.write(str(0))
            else:
                f.write(str(lsk_b[i1]))
                i1 += 1
            f.write('\n')
        f.close()

        f = open(path3+'/ε.txt','w')
        i1 = 0
        for i in range(len(ls4)+len(ls2_0)):
            f.write(ls[i+1][0])
            f.write('\t')
            if i in ls2_0:
                f.write(str(0))
            else:
                f.write(str(lsk_a[i1]))
                i1 += 1
            f.write('\n')
        f.close()


    f = open(path3+'/Library size.txt','w')
    l = []
    for i in range(len(lcount[0])):
        i1 = 0
        for j in range(len(lcount)):
            i1 += lcount[j][i]
        l.append(i1)
    if cell_number == len(ls7):
        for i in range(len(ls4p[0])):
            f.write('\t')
            f.write(str(ls[0][i]))
        f.write('\n')
    else:
        for i in range(len(lcount[0])):
            f.write('\t')
            f.write('cell'+str(i+1))
        f.write('\n')

    for i in range(len(lcount[0])):
        f.write('\t')
        f.write(str(int(l[i])))
    f.close()

    if cell_number == len(ls7):
        f = open(path3+'/Size factor.txt','w')

        for i in range(len(ls4p[0])):
            f.write('\t')
            f.write(str(ls[0][i]))
        f.write('\n')
        for i in range(len(ls4p[0])):
            f.write('\t')
            f.write(str(ls7[i]))
        f.close()
    else:
        f = open(path3+'/Size factor.txt','w')

        for i in range(len(ls7_1)):
            f.write('\t')
            f.write('cell'+str(i+1))
        f.write('\n')
        for i in range(len(ls7_1)):
            f.write('\t')
            f.write(str(ls7_1[i]))
        f.close()

    if param_batch != [1]:
        f = open(path3+'/Batch.txt','w')
        for i in range(len(ls7_1)):
            f.write('\t')
            if cell_number == len(ls7):
                f.write(str(ls[0][i]))
            else:
                f.write('cell'+str(i+1))
        f.write('\n')
        for i in range(len(ls7_1)):
            f.write('\t')
            f.write(str(batch_label[i]))
        f.close()

        f = open(path3+'/Batch factor.txt','w')
        for i in range(len(ls7_1)):
            f.write('\t')
            if cell_number == len(ls7):
                f.write(str(ls[0][i]))
            else:
                f.write('cell'+str(i+1))
        f.write('\n')
        i1 = 0
        for i in range(len(ls)-1):
            f.write(str(ls[i+1][0]))
            if i in ls2_0:
                for j in range(len(ls7_1)):
                    f.write('\t')
                    f.write(str(0))
            else:
                for j in range(len(ls7_1)):
                    f.write('\t')
                    f.write(str(batch_factor[i1][j]/param_library))
                i1 += 1
            f.write('\n')
        f.close()


    if param_group != [1]:
        if param_path == 'No':
            tree = 0
            for i in range(len(param_group)):
                if type(param_group[i]) == type([]):
                    tree = 1
            if tree == 0:
                f = open(path3+'/Group.txt','w')
                for i in range(len(ls7_1)):
                    f.write('\t')
                    if cell_number == len(ls7):
                        f.write(str(ls[0][i]))
                    else:
                        f.write('cell'+str(i+1))

                f.write('\n')
                for i in range(len(ls7_1)):
                    f.write('\t')
                    f.write(str(Group_label[i]))
                f.close()
            else:
                f = open(path3+'/Group.txt','w')
                for i in range(len(ls7_1)):
                    f.write('\t')
                    if cell_number == len(ls7):
                        f.write(str(ls[0][i]))
                    else:
                        f.write('cell'+str(i+1))
                f.write('\n')
                for i in range(len(ls7_1)):
                    f.write('\t')
                    if type(Group_label[i]) != type([]):
                        f.write(str(Group_label[i]))
                    else:
                        f.write(str(Group_label[i][0])+'_'+str(Group_label[i][1]))

                f.close()


            f = open(path3+'/DEgene factor.txt','w')
            for i in range(len(ls7_1)):
                f.write('\t')
                if cell_number == len(ls7):
                    f.write(str(ls[0][i]))
                else:
                    f.write('cell'+str(i+1))
            f.write('\n')
            i1 = 0
            ls1 = []
            for i in range(len(ls4)+len(ls2_0)):
                if i not in ls2_0:
                    ls1.append(i1)
                    i1 += 1
                else:
                    ls1.append(-1)

            for i in range(len(ls4)+len(ls2_0)):
                if ls1[i] != -1 and DEgene_label[ls1[i]] == 1:
                    f.write(ls[i+1][0])
                    for j in range(len(ls7_1)):
                        f.write('\t')
                        f.write(str(DEgene_factor[ls1[i]][j]))
                    f.write('\n')
            f.close()


            f = open(path3+'/Marker gene.txt','w')
            i1 = 0
            ls1 = []
            for i in range(len(ls4)+len(ls2_0)):
                if i not in ls2_0:
                    if i1 in l_marker:
                        f.write(ls[i+1][0])
                        f.write('\n')
                    i1 += 1
            f.close()


        else:
            tree = 0
            for i in range(len(param_group)):
                if type(param_group[i]) == type([]):
                    tree = 1
            if tree == 0:
                f = open(path3+'/Group.txt','w')
                for i in range(len(ls7_1)):
                    f.write('\t')
                    if cell_number == len(ls7):
                        f.write(str(ls[0][i]))
                    else:
                        f.write('cell'+str(i+1))
                f.write('\n')
                for i in range(len(ls7_1)):
                    f.write('\t')
                    f.write(str(Group_label[i]))
                f.close()
            else:
                f = open(path3+'/Group.txt','w')
                for i in range(len(ls7_1)):
                    f.write('\t')
                    if cell_number == len(ls7):
                        f.write(str(ls[0][i]))
                    else:
                        f.write('cell'+str(i+1))
                f.write('\n')
                for i in range(len(ls7_1)):
                    f.write('\t')
                    if type(Group_label[i]) != type([]):
                        f.write(str(Group_label[i]))
                    else:
                        f.write(str(Group_label[i][0])+'_'+str(Group_label[i][1]))
                f.close()

            f = open(path3+'/Path.txt','w')
            for i in range(len(ls7_1)):
                f.write('\t')
                if cell_number == len(ls7):
                    f.write(str(ls[0][i]))
                else:
                    f.write('cell'+str(i+1))
            f.write('\n')
            for i in range(len(ls7_1)):
                f.write('\t')
                f.write(str(Path[i]))
            f.close()

            f = open(path3+'/DEgene factor.txt','w')
            for i in range(len(ls7_1)):
                f.write('\t')
                if cell_number == len(ls7):
                    f.write(str(ls[0][i]))
                else:
                    f.write('cell'+str(i+1))
            f.write('\n')
            i1 = 0
            ls1 = []
            for i in range(len(ls4)+len(ls2_0)):
                if i not in ls2_0:
                    ls1.append(i1)
                    i1 += 1
                else:
                    ls1.append(-1)
            for i in range(len(ls4)+len(ls2_0)):
                if DEgene_label[ls1[i]] == 1:
                    f.write(ls[i+1][0])
                    for j in range(len(ls7_1)):
                        f.write('\t')
                        f.write(str(DEgene_factor[ls1[i]][j]))
                    f.write('\n')
            f.close()

            f = open(path3+'/noline gene factor.txt','w')
            for i in range(len(ls7_1)):
                f.write('\t')
                if cell_number == len(ls7):
                    f.write(str(ls[0][i]))
                else:
                    f.write('cell'+str(i+1))
            f.write('\n')
            i1 = 0
            ls1 = []
            for i in range(len(ls4)+len(ls2_0)):
                if i not in ls2_0:
                    ls1.append(i1)
                    i1 += 1
                else:
                    ls1.append(-1)
            for i in range(len(ls4)+len(ls2_0)):
                if Noline_gene[ls1[i]] == 1:
                    f.write(ls[i+1][0])
                    for j in range(len(ls7_1)):
                        f.write('\t')
                        f.write(str(DEgene_factor[ls1[i]][j]))
                    f.write('\n')
            f.close()





    if mod == '-d nocopula':
        f = open(path3+'/TrueCounts.txt','w')
        for i in range(len(ls7_1)):
            f.write('\t')
            if cell_number == len(ls7):
                f.write(str(ls[0][i]))
            else:
                f.write('cell'+str(i+1))
        f.write('\n')
        for i in range(len(ls4)+len(ls2_0)):
            f.write(ls[i+1][0])
            for j in range(len(ls7_1)):
                f.write('\t')
                f.write(str(int(ltruecount[i][j])))
            f.write('\n')
        f.close()

        f = open(path3+'/Dropout.txt','w')
        for i in range(len(ls7_1)):
            f.write('\t')
            if cell_number == len(ls7):
                f.write(str(ls[0][i]))
            else:
                f.write('cell'+str(i+1))
        f.write('\n')
        for i in range(len(ls4)+len(ls2_0)):
            f.write(ls[i+1][0])
            for j in range(len(ls7_1)):
                f.write('\t')
                f.write(str(ldp[i][j]))
            f.write('\n')
        f.close()
                











































def Remove0_t():
    """
    Read the file and remove the 0 gene expressed

    Parameters
    ----------
    path1:str
          Gene data file path
    path2:str
          Cell type data file path
    param:list
          Parameters used in simulation
          
    Returns
    -------
    ls2:list
          The data of gene
    ls2_a:list
          The type of each cell

    """
    f = open(path1,'r')
    if path1[-3:] == 'csv':
        reader = csv.reader(f)
        ls = list(reader)
    else:

        f1 = f.readlines()

        ls = []
        for i in f1:
            k = i.strip().split('\t')
            ls.append(k)
    if len(ls[0]) < len(ls[1]):
        ls[0].insert(0,' ')
    f.close()


    f = open(path2,'r')
    if path2[-3:] == 'csv':
        reader = csv.reader(f)
        ls1 = list(reader)
        if len(ls) == 1:
            ls_a = ls1[0]
        else:
            ls_a = ls1[1]
    else:
        f2 = f.readlines()
        ls_a = []
        if len(f2) == 1:
            ls_a = f2[0].strip().split('\t')
        else:
            ls_a = f2[1].strip().split('\t')


    print('Number of genes:{}'.format(len(ls)-1))
    global ls2,ls2_a
    ls2 = ls
    ls2_a = ls_a

def nml_t():
    """
    Normalize data with size factor

    Parameters
    ----------
    ls2:list
          The data of gene
    ls2_a:list
          The type of each cell
          
    Returns
    -------
    ls3:list
          The data of nomalized gene
    ls7:list
          Size factor
    ls3_a:list
          Number of cell types
    ls4_a:list
          Data segmented by cell group
    ls7_a:list
          Size factor segmented by cell group

    """
    global ls3,ls7
    ls3 = []
    ls7 = []
    ls1 = []
    for i in range(len(ls2[0])-1):
        sum1 = 0
        for j in range(len(ls2)-1):
            sum1 += float(ls2[j+1][i+1])
        ls1.append(sum1)
    ls1a = sorted(ls1)
    if len(ls1a)%2 == 0:
        k = (ls1a[len(ls1a)//2]+ls1a[len(ls1a)//2-1])/2
    else:
        k = ls1a[len(ls1a)//2]
    #print(k)
    for i in range(len(ls1)):
        ls7.append(ls1[i]/k)

    ls3.append(ls2[0])
    for i in range(len(ls2)-1):
        ls8 = []
        ls8.append(ls2[i+1][0])
        for j in range(len(ls2[i+1])-1):
            ls8.append(float(ls2[i+1][j+1])/ls7[j])
        ls3.append(ls8)

    
    global ls3_a,ls4_a,ls7_a,ls4_b,ls2_a
    ls3_a = []
    ls4_a = []
    ls7_a = []
    ls4_b = []
    
    for i in range(len(ls2_a)):
        if ls2_a[i] not in ls3_a:
            ls3_a.append(ls2_a[i])
    for i in range(len(ls3_a)):
        ls4_a.append([])
        ls4_b.append([])
    for i in range(len(ls3_a)):
        ls7_a.append([])
    ls_cell = []
    for i in range(len(ls3_a)):
        ls = []
        for j in range(len(ls7)):
            if ls2_a[j] == ls3_a[i]:
                ls.append(ls2[0][j])
        ls_cell.append(ls)
    ls_gene = []
    for j in range(len(ls2)-1):
        ls_gene.append(ls2[j+1][0])

    for i in range(len(ls3)-1):
        ls = []
        ls1 = []
        for i1 in range(len(ls3_a)):
            ls.append([])
            ls1.append([])
        for j in range(len(ls3[i])-1):
            for k in range(len(ls3_a)):
                if ls2_a[j] == ls3_a[k]:
                    ls[k].append(ls3[i+1][j+1])
                    ls1[k].append(ls2[i+1][j+1])
        for i2 in range(len(ls)):
            ls4_a[i2].append(ls[i2])
            ls4_b[i2].append(ls1[i2])
    for i in range(len(ls7)):
        for j in range(len(ls3_a)):
            if ls2_a[i] == ls3_a[j]:
                ls7_a[j].append(ls7[i])


    if Subgroup == 'Yes':
        dicta = {}

        for i in range(len(ls4_b)):
            warnings.filterwarnings("ignore")
            l = np.transpose(np.array(ls4_b[i]))
            ann = sc.AnnData(l)
            ann.var_names = ls_gene
            ann.obs_names = ls_cell[i]
            sc.pp.filter_genes(ann, min_cells=3)
            sc.pp.normalize_total(ann, target_sum=1e6)
            sc.pp.log1p(ann)
            sc.pp.highly_variable_genes(ann, min_mean=0.5, max_mean=15,min_disp=0.5)
            sc.tl.pca(ann, svd_solver='arpack')
            sc.pp.neighbors(ann, n_neighbors=20, n_pcs=15)
            sc.tl.leiden(ann,resolution=0.5)
            for j in range(len(ann.obs['leiden'])):
                dicta[ls_cell[i][j]] = str(ls3_a[i])+'_'+str(int(ann.obs['leiden'][j])+1)
        ls2_a = []
        for i in range(len(ls7)):
            ls2_a.append(dicta[ls2[0][i]])

        ls3_a = []
        ls4_a = []
        ls7_a = []
        ls4_b = []

        for i in range(len(ls2_a)):
            if ls2_a[i] not in ls3_a:
                ls3_a.append(ls2_a[i])
        for i in range(len(ls3_a)):
            ls4_a.append([])
            ls4_b.append([])
        for i in range(len(ls3_a)):
            ls7_a.append([])
        ls_cell = []
        for i in range(len(ls7)):
            ls = []
            for j in ls3_a:
                if ls2_a[i] == j:
                    ls.append(ls2[0][i])
            ls_cell.append(ls)
        ls_gene = []
        for j in range(len(ls2)-1):
            ls_gene.append(ls2[j+1][0])

        for i in range(len(ls3)-1):
            ls = []
            ls1 = []
            for i1 in range(len(ls3_a)):
                ls.append([])
                ls1.append([])
            for j in range(len(ls3[i])-1):
                for k in range(len(ls3_a)):
                    if ls2_a[j] == ls3_a[k]:
                        ls[k].append(ls3[i+1][j+1])
                        ls1[k].append(ls2[i+1][j+1])
            for i2 in range(len(ls)):
                ls4_a[i2].append(ls[i2])
                ls4_b[i2].append(ls1[i2])
        for i in range(len(ls7)):
            for j in range(len(ls3_a)):
                if ls2_a[i] == ls3_a[j]:
                    ls7_a[j].append(ls7[i])

    print('Number of cells:{}'.format(len(ls3[1])-1))


def mean_t():
    """
    Calculate mean of each gene

    Parameters
    ----------
    ls2:list
          The data of gene
    ls2_a:list
          The type of each cell

          
    Returns
    -------
    ls4:list
          Mean of each gene


    """
    global ls4
    ls4 = []
    for i in range(len(ls4_a)):
        ls = []
        for j in range(len(ls4_a[i])):
            ls.append(float(np.mean(ls4_a[i][j])))
        ls4.append(ls)

def BCV_t():
    """
    Fit the model of dispersion

    Parameters
    ----------
    ls4:list
          Mean of nor-DEgene
    ls4_b:list
          Mean of DEgene
    ls3:list
          The data of nomalized gene
    ls7:list
          Size factor
    ls3_a:list
          Number of cell types
    ls4_a:list
          Data segmented by cell group
    ls7_a:list
          Size factor segmented by cell group
    ls4:list
          Mean of each gene

    Returns
    -------
    ls12c:list
          Parameter dispersion of gene
    """
    
    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')


    from scipy.special import gammaln
    from scipy.special import digamma
    global ls2,ls3,ls7_a,ls4_a,correlation,dictx,ls12c,ls12a,ls12b,l_c
    ls12 = []
    for i in  range(len(ls4_a)):
        ls12.append([])
    dictx  ={}
    for i in range(len(ls4_a)):
        for j in range(len(ls4_a[i])):
            if ls4[i][j] == 0:
                ls12[i].append(0)
            else:
                ls12[i].append(float(np.var(ls4_a[i][j])/(ls4[i][j]**2)-1/ls4[i][j]))
    ls12d = []
    i1 = 0
    i2 = 0

    ls12a = []
    ls12b = []
    for i in range(len(ls4_a)):
        ls12a.append([])
        ls12b.append([])
    print('Estimate mean and dispersion...')
    warnings.filterwarnings("ignore")
    L = len(ls4_a[0])//100
    for i in range(len(ls4_a[0])):
        if i%L == 0:
            progress(i//L,width=50)
            #print('{}%'.format(i//L))
        for i1 in range(len(ls4_a)):
            if ls4[i1][i] != 0:
                i2 = 0
                i3 = 0
                i4 = 0
                for j in range(len(ls4_a[i1][i])):
                    i2 += (ls7_a[i1][j]*ls4[i1][i])**2
                    i3 += ls7_a[i1][j]
                    i4 += (float(ls4_b[i1][i][j])-ls7_a[i1][j]*ls4[i1][i])**2
                da = (i4-ls4[i1][i]*i3)/i2

                ls12b[i1].append(da)
            else:
                ls12b[i1].append(0)


            if ls12[i1][i] <= 0:
                d1 = 100
            else:
                d1 = ls12[i1][i]
            mu = ls4[i1][i]
            #
            #
            # l = list(np.random.normal(0,d1/3,15))
            # l1 = list(np.random.normal(0,d1/9,15))
            # l.extend(l1)
            #
            # k = 0
            # for j in range(len(ls4_b[i1][i])):
            #     if int(ls4_b[i1][i][j]) == 0:
            #         k += 1
            #
            # k1 = 0
            # for j in range(len(ls7_a[i1])):
            #     k1 += ((d1**(-1)/(ls7[i1][j]*mu+d1**(-1)))**(1/d1))
            # k2 = abs(k1-k)
            # d = d1
            # for j in l:
            #     if d+j > 0:
            #         d2 = d+j
            #     k1 = 0
            #     for j in range(len(ls7)):
            #         k1 += ((d2**(-1)/(ls7_a[i1][j]*mu+d2**(-1)))**(1/d2))
            #     if abs(k1-k) < k2:
            #         d = d2
            #         k2 = abs(k1-k)
            # ls12d.append(d)



            if ls4[i1][i] != 0:
                la = []
                lb = []
                lc = []
                ld = []
                for j in range(len(ls7_a[i1])):
                    la.append(mu*ls7_a[i1][j]*d1)
                    lb.append(1/(mu*ls7_a[i1][j]*d1)+1)
                    lc.append(1/d1+ls4_a[i1][i][j])
                    ld.append(1/d1)
                la1 = np.log(la)
                lb1 = np.log(lb)
                lc1 = digamma(lc)
                ld1 = digamma(ld)
                la1 = pd.Series(la1)
                lb1 = pd.Series(lb1)
                lc1 = pd.Series(lc1)
                ld1 = pd.Series(ld1)
                s1 = pd.Series(ls4_a[i1][i])
                s2 = pd.Series(ls7_a[i1])
                #lnF = 0
                d12 = d1**2
                d13 = d1**3
                s = 1/d12*la1-1/d12+1/d12*lb1+(1+s1*d1)/(d12+mu*s2*d13)-1/d12*lc1+1/d12*ld1
                lnF = float(sum(s))


                #for j in range(len(ls7)):
                #    lnF += 1/d12*la1[j]-1/d12+1/d12*lb1[j]+(1+ls4p[i][j]*d1)/(d12+mu*ls7[j]*d13)-1/d12*lc1[j]+1/d12*ld1[j]

                if lnF < 0:
                    K1 = 0
                    K2 = d1
                else:
                    K1 = d1
                    K2 = 1000
                for j in range(10):
                    d1 = (K1+K2)/2
                    la = []
                    lb = []
                    lc = []
                    ld = []
                    for i2 in range(len(ls7_a[i1])):
                        la.append(mu*ls7_a[i1][i2]*d1)
                        lb.append(1/(mu*ls7_a[i1][i2]*d1)+1)
                        lc.append(1/d1+ls4_a[i1][i][i2])
                        ld.append(1/d1)
                    la1 = np.log(la)
                    lb1 = np.log(lb)
                    lc1 = digamma(lc)
                    ld1 = digamma(ld)
                    la1 = pd.Series(la1)
                    lb1 = pd.Series(lb1)
                    lc1 = pd.Series(lc1)
                    ld1 = pd.Series(ld1)
                    s1 = pd.Series(ls4_a[i1][i])
                    s2 = pd.Series(ls7_a[i1])
                    #lnF = 0
                    d12 = d1**2
                    d13 = d1**3
                    s = 1/d12*la1-1/d12+1/d12*lb1+(1+s1*d1)/(d12+mu*s2*d13)-1/d12*lc1+1/d12*ld1
                    lnF = float(sum(s))

                    #for i1 in range(len(ls7)):
                    #    lnF += 1/d12*la1[i1]-1/d12+1/d12*lb1[i1]+(1+ls4p[i][i1]*d1)/(d12+mu*ls7[i1]*d13)-1/d12*lc1[i1]+1/d12*ld1[i1]
                    if lnF > 0:
                        K1 = d1
                    else:
                        K2 = d1

                ls12a[i1].append((K1+K2)/2)
            else:
                ls12a[i1].append(0)

    ls4_1 = []
    ls12a_1 = []
    ls12b_1 = []
    for i in range(len(ls4)):
        l1 = []
        l2 = []
        l3 = []
        for j in range(len(ls4[i])):
            if ls4[i][j] != 0:
                l1.append(ls4[i][j])
                l2.append(ls12a[i][j])
                l3.append(ls12b[i][j])
        ls4_1.append(l1)
        ls12a_1.append(l2)
        ls12b_1.append(l3)

    ls4_2 = []

    ls4_2 = [list(np.log(i)) for i in ls4_1]

    def funb(x,t1,t2,t3,d0):
        return t1/(t2+np.exp(t3*x))+d0

    ls12c = []
    for i in range(len(ls4_2)):
        ls12c.append([])

    for i in range(len(ls4_2)):

        ls = sorted(ls4_2[i])
        l = (ls[-1]+0.001-ls[0])/20
        l1 = []
        l2 = []
        for i1 in range(20):
            l1.append([])
            l2.append([])

        for i1 in range(len(ls4_2[i])):
            l1[int((ls4_2[i][i1]-ls[0])/l)].append(ls4_2[i][i1])
            l2[int((ls4_2[i][i1]-ls[0])/l)].append(ls12a_1[i][i1])
        l3 = []
        l4 = []
        for i1 in range(len(l1)):
            l2[i1].sort()
            if  len(l1[i1]) > 10:
                l3.append(float(np.mean(l1[i1])))
                l4.append(l2[i1][int(len(l2[i1])/2)])

        l5 = []
        l6 = []
        for i1 in range(20):
            l5.append([])
            l6.append([])
        for i1 in range(len(ls4_2[i])):
            l5[int((ls4_2[i][i1]-ls[0])/l)].append(ls4_2[i][i1])
            l6[int((ls4_2[i][i1]-ls[0])/l)].append(ls12b_1[i][i1])
        l7 = []
        l8 = []
        for i1 in range(len(l5)):
            l6[i1].sort()
            if  len(l1[i1]) > 10:
                l7.append(float(np.mean(l5[i1])))
                l8.append(l6[i1][int(len(l6[i1])/2)])








        ls1 = []
        for i1 in range(len(l4)):
            if l4[i1] != 0:
                ls1.append(l8[i1]/l4[i1])
            else:
                ls1.append(0)
        x = np.array(l3[8:])
        y = np.array(ls1[8:])
        param_bounds=([0,0,-np.inf,0],[np.inf,np.inf,np.inf,np.inf])
        try:
            popt, pcov = curve_fit(funb,x,y,bounds=param_bounds)
        except:
            param_bounds=([0,0,-1000,0],[20000,20000,1000,100])
            popt, pcov = curve_fit(funb,x,y,bounds=param_bounds)
        T1 = popt[0]
        T2 = popt[1]
        T3 = popt[2]
        D0 = popt[3]
        for i1 in range(len(ls4[i])):
            if ls4[i][i1] != 0:
                ls12c[i].append((T1/(T2+np.exp(T3*np.log(ls4[i][i1])))+D0)*ls12a[i][i1])
            else:
                ls12c[i].append(0)




    correlation = []
    dictx = []
    l_c = []
    for i in range(len(ls4)):
        l2 = []
        l3 = []
        l_c1 = []
        dictx1 = {}
        l2 = list(reversed(list(np.argsort(ls4[i]))))
        #print(ls4[i][l2[1999]])
        if mod != '-e nocopula':
            if len(l2) >= copula_number:
                for i1 in range(copula_number):
                    if max(ls4_a[i][l2[i1]]) > 0:
                        l3.append(ls4_a[i][l2[i1]])
                        l_c1.append(l2[i1])
                for i1 in range(len(l_c1)):
                    dictx1[l_c1[i1]] = i1
            else:
                for i1 in range(len(l2)):
                    if max(ls4_a[i][l2[i1]]) > 0:
                        l3.append(ls4_a[i][l2[i1]])
                        l_c1.append(l2[i1])
                for i1 in range(len(l_c1)):
                    dictx1[l_c1[i1]] = i1
            Mat = np.array(l3)
            correlation2 = np.corrcoef(Mat)
            correlation.append(correlation2)
            dictx.append(dictx1)
            l_c.append(l_c1)
            #sns.set()
            #ax = sns.clustermap(correlation2,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
            #plt.show()

            Mat = np.array(l3[:50])
            correlation1 = np.corrcoef(Mat)
            #sns.set()
            #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
            #plt.show()
        else:
            l_c.append([])



def creatcount_noparam_t():
    """

    Simulate count

    Parameters
    ----------
    ls4:list
          Mean of gene
    ls12c:list
          Dispersion of gene
    ls7:list
          Size factor


    Returns
    -------
    llisize:list
          Library size
    lcount:list
          Count

    """
    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')

    from scipy.stats import nbinom
    from scipy.stats import norm
    global lcount,ls4,correlation,ls12c,ls7,dictx,l_c,ls7_a,ls4_a,ls7_1,ls5a

    print(' ')
    print('Start simulation...')
    warnings.filterwarnings("ignore")
    ls5a = []
    ls7_1 = []
    ls1 = [len(i) for i in ls7_a]
    if cell_number == sum(ls1):
        ls7_1 = ls7_a
    else:
        if cell_number < sum(ls1):
            for i in range(len(ls7_a)):
                if i < len(ls7_a)-1:
                    if int(len(ls7_a[i])*cell_number/sum(ls1)) == 0:
                        ls7_1.append([ls7_a[i][0]])
                        ls5a.append(1)
                    else:
                        ls7_1.append(ls7_a[i][:int(len(ls7_a[i])*cell_number/sum(ls1))])
                        ls5a.append(int(len(ls7_a[i])*cell_number/sum(ls1)))
                else:
                    ls1a = [len(j) for j in ls7_1]
                    k = sum(ls1a)
                    k1 = cell_number-k
                    if k1 > 0:
                        ls7_1.append(ls7_a[i][:k1])
                        ls5a.append(k1)
                    else:
                        ls7_1.append(ls7_a[i][0])
                        ls5a.append(1)
                    if len(ls7_a[i]) < k1:
                        while 1:
                            ls7_1[-1].append(ls7_a[i][0])
                            if len(ls7_1[-1]) == k1:
                                break


        else:
            for i in range(len(ls7_a)):
                ls7_2 = np.log(ls7_a[i])
                x_1,x_2 = norm.fit(ls7_2)
                if i < len(ls7_a)-1:
                    lsa = [float(j) for j in ls7_a[i]]
                    while 1:
                        p = float(np.exp(np.random.normal(x_1,x_2,1)))
                        if len(lsa) == int(len(ls7_a[i])*cell_number/sum(ls1)):
                            break
                        lsa.append(p)
                    ls7_1.append(lsa)
                    ls5a.append(len(lsa))
                else:
                    ls1a = [len(j) for j in ls7_1]
                    k = sum(ls1a)
                    k1 = cell_number-k
                    lsa = [j for j in ls7_a[i]]
                    if len(lsa) <= k1:
                        while 1:
                            p = float(np.exp(np.random.normal(x_1,x_2,1)))
                            if len(lsa) == k1:
                                break
                            lsa.append(p)
                    else:
                        lsa = lsa[:k1]
                    ls7_1.append(lsa)
                    ls5a.append(len(lsa))

    if param_library != 1:
        ls7_3 = [i for i in ls7_1]
        ls7_1 = []
        for i in range(len(ls7_3)):
            l = [j*param_library for j in ls7_3[i]]
            ls7_1.append(l)






    if mod != '-e nocopula':
        p = []
        for i in range(len(ls4)):
            l5 = []
            for i1 in range(len(correlation[i])):
                l5.append(0)
            p1 = np.random.multivariate_normal(mean=l5,cov=correlation[i],size=len(ls7_1[i]))


            k  = list(np.transpose(p1))
            p.append(k)
            # k1 = np.corrcoef(np.array(k))
            # sns.set()
            # ax = sns.clustermap(k1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
            # plt.show()
            #
            # k1 = np.corrcoef(np.array(k[:50]))
            # sns.set()
            # ax = sns.clustermap(k1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
            # plt.show()





    lm = []
    lcount = []
    for i in range(len(ls4)):
        lcount.append([])
        lm.append([])

    s = []
    for i in range(len(ls7_a)):
        s.append(np.array(ls7_a[i]))

    LJ = 0
    for i in range(len(ls4_a[0])):
        LJ += ls4_a[0][i][0]

    N = sum(ls1)
    i1 = 0
    i2 = (len(ls4[0]))//100
    for i in range(len(ls4[0])):
        if i%i2 == 0:
            progress(int(i/i2),width=50)
            #print('{}%'.format(int(i/i2)))
        if cell_number == N and param_library == 1:
            for i1 in range(len(ls4)):
                l = []
                if ls4[i1][i] == 0:
                    for j in range(len(ls4_a[i1][0])):
                        l.append(0)
                else:
                    if i in l_c[i1]:
                        l1 = []
                        l2 = []
                        l4 = []
                        l5 = []
                        for j in range(len(ls4_a[i1][0])):
                            l1.append(1/ls12c[i1][i])
                            l2.append(1/(ls12c[i1][i]*ls4[i1][i]*ls7_a[i1][j]+1))
                        l4 = p[i1][dictx[i1][i]]
                        l5 = list(norm.cdf(l4))
                        l = list(nbinom.ppf(l5, l1, l2))

                        l6 = []
                        for i3 in range(len(l)):
                            l6.append(l[i3]/ls7_a[i1][i3])
                        if max(l6) > 0:
                            lm[i1].append(l6)



                    else:
                        ls1 = np.random.gamma(1/ls12c[i1][i],ls4[i1][i]*ls12c[i1][i],len(ls4_a[i1][0]))

                        l1 = []
                        for j in range(len(ls1)):
                            l1.append(float(ls1[j]*ls7_a[i1][j]))
                        kkk = 0
                        while 1:
                            kkk += 1
                            x = np.random.poisson(tuple(l1),(1,len(ls4_a[i1][0])))
                            ls1 = list(map(int,x[0]))
                            if max(ls1) > 0 or kkk > 20:
                                break
                        l = ls1

                    l7 = []
                    l8 = list(np.random.normal(0,ls12c[i1][i]/10,10))
                    l7.extend(l8)
                    l9 = list(np.random.normal(0,ls12c[i1][i]/20,20))
                    l7.extend(l9)
                    mean1 = np.mean(np.log2(np.array(ls4_a[i1][i])*1000000/LJ+1))
                    mean2 = np.mean(np.log2(np.array(l)/s[i1]*1000000/LJ+1))
                    L = abs(mean2-mean1)
                    for i_1 in range(len(l7)):
                        if ls12c[i1][i] + l7[i_1] > 0:
                            if i in l_c[i1]:
                                l1 = []
                                l2 = []
                                l4 = []
                                l5 = []
                                for j in range(len(ls4_a[i1][0])):
                                    l1.append(1/(ls12c[i1][i]+l7[i_1]))
                                    l2.append(1/((ls12c[i1][i]+l7[i_1])*ls4[i1][i]*ls7_a[i1][j]+1))
                                l4 = p[i1][dictx[i1][i]]
                                l5 = list(norm.cdf(l4))
                                la = list(nbinom.ppf(l5, l1, l2))

                                if abs(np.mean(np.log2(np.array(la)/s[i1]*1000000/LJ+1))-mean1) < L:
                                    L = abs(np.mean(np.log2(np.array(la)/s[i1]*1000000/LJ+1))-mean1)
                                    l = la
                                    ls12c[i1][i] = (ls12c[i1][i] + l7[i_1])





                            else:
                                ls1 = np.random.gamma(1/(ls12c[i1][i]+l7[i_1]),ls4[i1][i]*(ls12c[i1][i]+l7[i_1]),len(ls4_a[i1][0]))

                                l1 = []
                                for j in range(len(ls1)):
                                    l1.append(float(ls1[j]*ls7_a[i1][j]))
                                kkk = 0
                                while 1:
                                    kkk += 1
                                    x = np.random.poisson(tuple(l1),(1,len(ls4_a[i1][0])))
                                    ls1 = list(map(int,x[0]))
                                    if max(ls1) > 0 or kkk > 20:
                                        break
                                la = ls1
                                if abs(np.mean(np.log2(np.array(la)/s[i1]*1000000/LJ+1))-mean1) < L:
                                    L = abs(np.mean(np.log2(np.array(la)/s[i1]*1000000/LJ+1))-mean1)
                                    l = la
                                    ls12c[i1][i] = (ls12c[i1][i]+l7[i_1])
                lcount[i1].append(l)
        else:
            for i1 in range(len(ls4)):
                l = []
                if ls4[i1][i] == 0:
                    for j in range(len(ls7_1[i1])):
                        l.append(0)
                else:
                    if i in l_c[i1]:
                        l1 = []
                        l2 = []
                        l4 = []
                        l5 = []
                        for j in range(len(ls7_1[i1])):
                            l1.append(1/ls12c[i1][i])
                            l2.append(1/(ls12c[i1][i]*ls4[i1][i]*ls7_1[i1][j]+1))
                        l4 = p[i1][dictx[i1][i]]
                        l5 = list(norm.cdf(l4))
                        l = list(nbinom.ppf(l5, l1, l2))

                        l6 = []
                        for i3 in range(len(l)):
                            l6.append(l[i3]/ls7_1[i1][i3])
                        if max(l6) > 0:
                            lm[i1].append(l6)



                    else:
                        ls1 = np.random.gamma(1/ls12c[i1][i],ls4[i1][i]*ls12c[i1][i],len(ls7_1[i1]))

                        l1 = []
                        for j in range(len(ls1)):
                            l1.append(float(ls1[j]*ls7_1[i1][j]))
                        kkk = 0
                        while 1:
                            kkk += 1
                            x = np.random.poisson(tuple(l1),(1,len(ls7_1[i1])))
                            ls1 = list(map(int,x[0]))
                            if max(ls1) > 0 or kkk > 20:
                                break
                        l = ls1


                lcount[i1].append(l)




    print(' ')
    print('Simulation completed')
    if mod != '-e nocopula':
        for i in range(len(lm)):
            Mat = np.array(lm[i])
            correlation1 = np.corrcoef(Mat)

            #sns.set()
            #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
            #plt.show()

            ln = []
            for i1 in range(len(lm[i])):
                ln.append(float(np.mean(lm[i][i1])))
            ln1 = []
            ln1 = list(reversed(list(np.argsort(ln))))
            ln2 = []
            for i1 in range(len(lm[i])):
                if i1 in ln1[:50]:
                    ln2.append(lm[i][i1])
            Mat = np.array(ln2)
            correlation1 = np.corrcoef(Mat)

            #sns.set()
            #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
            #plt.show()

def NBzero_mean_dispersion_noparam_t():
    """
    Estimate mean and dispersion of -f or -f nocopula

    Parameters
    ----------
    ls3:list
          The data of nomalized gene

    Returns
    -------
    ls4:list
          Mean of gene
    ls12c:list
          Dispersion of gene

    """
    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')



    global ls4,ls4_a,ls4_b,ls7_a,dictx,ls12c,correlation,l_c,k,k_a,lsk_a,lsk_b

    lsk_a = []
    lsk_b = []
    lme = []
    ldis = []
    for i in range(len(ls4)):
        lsk_a.append([])
        lsk_b.append([])
        lme.append([])
        ldis.append([])

    warnings.filterwarnings("ignore")

    print('Estimate mean,dispersion,γ,ε...')
    for i_1 in range(len(ls4)):
        k3 = len(ls4[i_1])//100
        s = np.array(ls7_a[i_1])
        print(str(i_1+1)+'/'+str(len(ls4)))
        for i in range(len(ls4_a[i_1])):
            if i%k3 == 0:
                progress(int(i/k3),width=50)


            m0 = float(np.mean(ls4_a[i_1][i]))
            p0 = ls4_a[i_1][i].count(0)/len(ls4_a[i_1][i])
            v0 = float(np.var(ls4_a[i_1][i]))
            if m0 != 0:
                d0 = v0/(m0**2)


                ##########################################################################


                l1 = []
                l2 = []
                l3 = []
                l4 = []
                step1 = (1/(1-p0))**(0.2)
                for i1 in range(7):
                    for j1 in range(5):
                        for k1 in range(5):
                            for k2 in range(10):
                                if i1 == 0:
                                    l1.append([m0])
                                    l2.append([0.6*d0+j1*0.1*d0])
                                    l3.append([0.1+0.2*k1])
                                    l4.append([-5+k2*0.56])
                                elif i1 == 1:
                                    if (step1-1)*m0*0.1+m0 > 0:
                                        l1.append([(step1-1)*m0*0.1+m0])
                                        l2.append([0.6*d0+j1*0.1*d0])
                                        l3.append([0.1+0.2*k1])
                                        l4.append([-5+k2*0.56])
                                else:
                                    l1.append([m0*(step1**(i1-1))])
                                    l2.append([0.6*d0+j1*0.1*d0])
                                    l3.append([0.1+0.2*k1])
                                    l4.append([-5+k2*0.56])
                l1 = np.array(l1)
                l2 = np.array(l2)
                l3 = np.array(l3)
                l4 = np.array(np.exp(l4))
                d = l2
                m = l1
                p1 = l3
                k = l4
                g = np.exp(gammaln(1/d+1)-gammaln(1/d))
                s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
                if p0 > 0:
                    s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
                g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
                m1 = np.transpose(np.array(s2.mean(axis=1)))
                ls = []
                for i1 in m1:
                    ls.append([i1])
                m1  = np.array(ls)

                s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
                s4 = s3.mean(axis=1)
                d1 = []
                for i1 in s4:
                    d1.append([i1])
                d1 = np.array(d1)
                d1 = d1-m1**2
                if p0 > 0:
                    s5 = s1.mean(axis=1)
                    ls = []
                    for i1 in s5:
                        ls.append([i1])
                    p1 = np.array(ls)

                if p0 > 0:

                    if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                        l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                    else:
                        l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

                else:
                    l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                r = np.argmin(l)
                m = l1[r][0]
                d = l2[r][0]
                p1 = l3[r][0]
                k = np.log(l4[r][0])



                l1 = []
                l2 = []
                l3 = []
                l4 = []
                step1 = (1/(1-p0))**(0.2)
                for i1 in range(5):
                    for j1 in range(5):
                        for k1 in range(5):
                            for k2 in range(5):
                                if p1+(k1-2)*0.1 >= 0 and p1+(k1-2)*0.1 <=1:
                                    if m == m0:
                                        if (step1-1)*m0*0.1*i1*0.25+m0 > 0:
                                            l1.append([(step1-1)*m0*0.1*i1*0.25+m0])
                                            l2.append([d+(j1-2)*0.05*d0])
                                            l3.append([p1+(k1-2)*0.1])
                                            l4.append([k+(k2-2)*0.28])
                                    elif m == (step1-1)*m0*0.1+m0:
                                        if (step1-1)*m0*0.1*i1*0.5+m0 > 0:
                                            l1.append([(step1-1)*m0*0.1*i1*0.5+m0])
                                            l2.append([d+(j1-2)*0.05*d0])
                                            l3.append([p1+(k1-2)*0.1])
                                            l4.append([k+(k2-2)*0.28])
                                    else:
                                        l1.append([m*(step1**((i1-2)*0.5))])
                                        l2.append([d+(j1-2)*0.05*d0])
                                        l3.append([p1+(k1-2)*0.1])
                                        l4.append([k+(k2-2)*0.28])
                l1 = np.array(l1)
                l2 = np.array(l2)
                l3 = np.array(l3)
                l4 = np.array(np.exp(l4))
                d = l2
                m = l1
                p1 = l3
                k = l4
                g = np.exp(gammaln(1/d+1)-gammaln(1/d))
                s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
                if p0 > 0:
                    s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
                g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
                m1 = np.transpose(np.array(s2.mean(axis=1)))
                ls = []
                for i1 in m1:
                    ls.append([i1])
                m1  = np.array(ls)

                s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
                s4 = s3.mean(axis=1)
                d1 = []
                for i1 in s4:
                    d1.append([i1])
                d1 = np.array(d1)
                d1 = d1-m1**2
                if p0 > 0:
                    s5 = s1.mean(axis=1)
                    ls = []
                    for i1 in s5:
                        ls.append([i1])
                    p1 = np.array(ls)

                if p0 > 0:

                    if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                        l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                    else:
                        l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

                else:
                    l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                r = np.argmin(l)
                m = l1[r][0]
                d = l2[r][0]
                p1 = l3[r][0]
                k = np.log(l4[r][0])





                l1 = []
                l2 = []
                l3 = []
                l4 = []
                step1 = (1/(1-p0))**0.2
                for i1 in range(5):
                    for j1 in range(5):
                        for k1 in range(5):
                            for k2 in range(5):
                                if p1+(k1-2)*0.1 >= 0 and p1+(k1-2)*0.1 <=1:
                                    if m == m0:
                                        if (step1-1)*m0*0.1*i1*0.125+m0 > 0:
                                            l1.append([(step1-1)*m0*0.1*i1*0.125+m0])
                                            l2.append([d+(j1-2)*0.025*d0])
                                            l3.append([p1+(k1-2)*0.05])
                                            l4.append([k+(k2-2)*0.14])
                                    elif m <= (step1-1)*m0*0.2+m0:
                                        if (step1-1)*m0*0.05*(i1-2)*0.25+m > 0:
                                            l1.append([(step1-1)*m0*0.05*(i1-2)*0.25+m])
                                            l2.append([d+(j1-2)*0.025*d0])
                                            l3.append([p1+(k1-2)*0.05])
                                            l4.append([k+(k2-2)*0.14])
                                    else:
                                        l1.append([m*(step1**((i1-2)*0.25))])
                                        l2.append([d+(j1-2)*0.025*d0])
                                        l3.append([p1+(k1-2)*0.05])
                                        l4.append([k+(k2-2)*0.14])
                l1 = np.array(l1)
                l2 = np.array(l2)
                l3 = np.array(l3)
                l4 = np.array(np.exp(l4))
                d = l2
                m = l1
                p1 = l3
                k = l4
                g = np.exp(gammaln(1/d+1)-gammaln(1/d))
                s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
                if p0 > 0:
                    s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
                g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
                m1 = np.transpose(np.array(s2.mean(axis=1)))
                ls = []
                for i1 in m1:
                    ls.append([i1])
                m1  = np.array(ls)

                s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
                s4 = s3.mean(axis=1)
                d1 = []
                for i1 in s4:
                    d1.append([i1])
                d1 = np.array(d1)
                d1 = d1-m1**2
                if p0 > 0:
                    s5 = s1.mean(axis=1)
                    ls = []
                    for i1 in s5:
                        ls.append([i1])
                    p1 = np.array(ls)

                if p0 > 0:

                    if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                        l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                    else:
                        l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

                else:
                    l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                r = np.argmin(l)
                m = l1[r][0]
                d = l2[r][0]
                p1 = l3[r][0]
                k = np.log(l4[r][0])



                l1 = []
                l2 = []
                l3 = []
                l4 = []
                step1 = (1/(1-p0))**0.2
                for i1 in range(5):
                    for j1 in range(5):
                        for k1 in range(5):
                            for k2 in range(5):
                                if p1+(k1-2)*0.1 >= 0 and p1+(k1-2)*0.1 <=1:
                                    if m == m0:
                                        if (step1-1)*m0*0.1*i1*0.0625+m0 > 0:
                                            l1.append([(step1-1)*m0*0.1*i1*0.0625+m0])
                                            l2.append([d+(j1-2)*0.0125*d0])
                                            l3.append([p1+(k1-2)*0.025])
                                            l4.append([k+(k2-2)*0.07])
                                    elif m <= (step1-1)*m0*0.2+m0:
                                        if (step1-1)*m0*0.05*(i1-2)*0.125+m > 0:
                                            l1.append([(step1-1)*m0*0.05*(i1-2)*0.125+m])
                                            l2.append([d+(j1-2)*0.0125*d0])
                                            l3.append([p1+(k1-2)*0.025])
                                            l4.append([k+(k2-2)*0.07])
                                    else:
                                        l1.append([m*(step1**((i1-2)*0.125))])
                                        l2.append([d+(j1-2)*0.0125*d0])
                                        l3.append([p1+(k1-2)*0.025])
                                        l4.append([k+(k2-2)*0.07])
                l1 = np.array(l1)
                l2 = np.array(l2)
                l3 = np.array(l3)
                l4 = np.array(np.exp(l4))
                d = l2
                m = l1
                p1 = l3
                k = l4
                g = np.exp(gammaln(1/d+1)-gammaln(1/d))
                s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))
                if p0 > 0:
                    s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
                g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
                m1 = np.transpose(np.array(s2.mean(axis=1)))
                ls = []
                for i1 in m1:
                    ls.append([i1])
                m1  = np.array(ls)

                s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
                s4 = s3.mean(axis=1)
                d1 = []
                for i1 in s4:
                    d1.append([i1])
                d1 = np.array(d1)
                d1 = d1-m1**2
                if p0 > 0:
                    s5 = s1.mean(axis=1)
                    ls = []
                    for i1 in s5:
                        ls.append([i1])
                    p1 = np.array(ls)

                if p0 > 0:

                    if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                        l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                    else:
                        l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

                else:
                    l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                r = np.argmin(l)












                M = abs(m0-m1)[r]
                D = (abs(v0-d1)/v0)[r]
                if p0 > 0:
                    P = (abs(p0-p1)/p0)[r]
                else:
                    P = 0.05
                r_1 = r
                m1_1 = m1
                d1_1 = d1
                p1_1 = p1

                m = l1[r][0]
                d = l2[r][0]
                p1 = l3[r][0]
                k = l4[r][0]


                L = min(l)[0]
                L1 = L


                m_result = m
                d_result = d
                p1_result = p1
                k_result = k


                for i2 in range(20):

                    if D < 0.05:
                        if M > 0:
                            l1 = np.random.normal(0,M,64)
                            l2 = np.random.normal(0,0.025*d_result,64)
                            if P*p0 > 0:
                                l3 = np.random.normal(0,P*p0,64)
                            else:
                                l3 = np.random.normal(0,p1_result*0.001,64)
                            if (m_result*p1_result) <= 0 or P*p0 <= 0:
                                l4 = np.random.normal(0,k_result*0.1,64)
                            else:
                                if P*p0/(m_result*p1_result) > 0.001 and P*p0/(m_result*p1_result) < 10:
                                    l4 = np.random.normal(0,P*p0/(m_result*p1_result),64)
                                else:
                                    l4 = np.random.normal(0,k_result*0.1,64)
                        else:
                            l1 = np.random.normal(0,m_result*0.001,64)
                            l2 = np.random.normal(0,0.025*d_result,64)
                            if P*p0 > 0:
                                l3 = np.random.normal(0,P*p0,64)
                            else:
                                l3 = np.random.normal(0,p1_result*0.001,64)
                            if (m_result*p1_result) <= 0 or P*p0 <= 0:
                                l4 = np.random.normal(0,k_result*0.1,64)
                            else:
                                if P*p0/(m_result*p1_result) > 0.001 and P*p0/(m_result*p1_result) < 10:
                                    l4 = np.random.normal(0,P*p0/(m_result*p1_result),64)
                                else:
                                    l4 = np.random.normal(0,k_result*0.1,64)
                    else:
                        l1 = np.random.normal(0,0.01*m_result,64)
                        l2 = np.random.normal(0,D*d_result,64)
                        if P*p0 > 0:
                            l3 = np.random.normal(0,P*p0,64)
                        else:
                            l3 = np.random.normal(0,p1_result*0.001,64)
                        if (m_result*p1_result) <= 0 or P*p0 <= 0:
                            l4 = np.random.normal(0,k_result*0.1,64)
                        else:
                            if P*p0/(m_result*p1_result) > 0.001 and P*p0/(m_result*p1_result) < 10:
                                l4 = np.random.normal(0,P*p0/(m_result*p1_result),64)
                            else:
                                l4 = np.random.normal(0,k_result*0.1,64)


                    l1a = []
                    l2a = []
                    l3a = []
                    l4a = []

                    for j in range(64):
                        if m_result+l1[j] > 0 and d_result+l2[j] > 0 and p1_result+l3[j] > 0 and p1_result+l3[j] <= 1 and k_result+l4[j] > 0:
                            l1a.append([m_result+l1[j]])
                            l2a.append([d_result+l2[j]])
                            l3a.append([p1_result+l3[j]])
                            l4a.append([k_result+l4[j]])
                    l1 = np.array(l1a)
                    l2 = np.array(l2a)
                    l3 = np.array(l3a)
                    l4 = np.array(l4a)

                    if len(l1a) > 0:

                        d = l2
                        m = l1
                        p1 = l3
                        k = l4
                        g = np.exp(gammaln(1/d+1)-gammaln(1/d))


                        s2 = d*m*g*(1-p1*(1+k*m*s*d)**(-1/d-1))


                        if p0 > 0:
                            s1 = (1+d*m*s)**(-1/d)*(1-p1*(k/(1/(d*m*s)+1)+1)**(-1/d))+(1+k*d*m*s)**(-1/d)*p1
                        g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
                        m1 = np.transpose(np.array(s2.mean(axis=1)))
                        ls = []
                        for i1 in m1:
                            ls.append([i1])
                        m1  = np.array(ls)

                        s3 = ((d*m)**2*g1*(1-(1+k*d*m*s)**(-1/d-2)*p1)+(d*m/s)*g*(1-(1+k*d*m*s)**(-1/d-1)*p1))
                        s4 = s3.mean(axis=1)
                        d1 = []
                        for i1 in s4:
                            d1.append([i1])
                        d1 = np.array(d1)
                        d1 = d1-m1**2
                        if p0 > 0:
                            s5 = s1.mean(axis=1)
                            ls = []
                            for i1 in s5:
                                ls.append([i1])
                            p1 = np.array(ls)



                        if p0 > 0:

                            if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.25+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                                l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                            else:
                                l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)

                        else:
                            l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)



                        if min(l)[0] < L:
                            r = np.argmin(l)
                            da = d1[r]
                            pa = p1[r]
                            ma = m1[r]

                            M = abs(m0-m1)[r]
                            D = (abs(v0-d1)/v0)[r]
                            if p0 > 0:
                                P = (abs(p0-p1)/p0)[r]
                            else:
                                P = 0.05
                            m_result = l1[r][0]
                            d_result = l2[r][0]
                            p1_result = l3[r][0]
                            k_result = l4[r][0]
                            L = min(l)[0]
                # print(D,M,P)
                # if L1 != L:
                #     print('Yes')
                #     if p0 > 0:
                #         print(L,m0,ma,v0,da,p0,pa)
                #         print(p1_result,k_result,m_result,d_result)
                #     else:
                #         print(L,m0,ma,v0,da)
                #         print(p1_result,k_result,m_result,d_result)
                # else:
                #     print('No')
                #     if p0 > 0:
                #         print(L1,m0,m1_1[r_1],v0,d1_1[r_1],p0,p1_1[r_1])
                #         print(p1_result,k_result,m_result,d_result)
                #     else:
                #         print(L1,m0,m1_1[r_1],v0,d1_1[r_1])
                #         print(p1_result,k_result,m_result,d_result)




                l1 = []
                l2 = []

                step1 = (1/(1-p0))**0.2
                for i1 in range(11):
                    for j1 in range(11):
                        if m_result == m0:
                            if (step1-1)*m0*0.1*i1*0.03125+m0 > 0 and d_result+(j1-5)*0.00625*d0 > 0:
                                l1.append([(step1-1)*m0*0.1*i1*0.03125+m0])
                                l2.append([d_result+(j1-5)*0.00625*d0])

                        elif m_result <= (step1-1)*m0*0.2+m0:
                            if (step1-1)*m0*0.025*(i1-5)*0.125+m_result > 0 and d_result+(j1-5)*0.00625*d0 > 0:
                                l1.append([(step1-1)*m0*0.025*(i1-5)*0.125+m_result])
                                l2.append([d_result+(j1-5)*0.00625*d0])

                        else:
                            if m_result*(step1**((i1-5)*0.0625)) > 0 and d_result+(j1-5)*0.00625*d0 > 0:
                                l1.append([m_result*(step1**((i1-5)*0.0625))])
                                l2.append([d_result+(j1-5)*0.00625*d0])
                if len(l1) > 0:
                    l1 = np.array(l1)
                    l2 = np.array(l2)

                    d = l2
                    m = l1

                    g = np.exp(gammaln(1/d+1)-gammaln(1/d))
                    s2 = d*m*g*(1-p1_result*(1+k_result*m*s*d)**(-1/d-1))
                    if p0 > 0:
                        s1 = (1+d*m*s)**(-1/d)*(1-p1_result*(k_result/(1/(d*m*s)+1)+1)**(-1/d))+(1+k_result*d*m*s)**(-1/d)*p1_result
                    g1 = np.exp(gammaln(1/d+2)-gammaln(1/d))
                    m1 = np.transpose(np.array(s2.mean(axis=1)))
                    ls = []
                    for i1 in m1:
                        ls.append([i1])
                    m1  = np.array(ls)

                    s3 = ((d*m)**2*g1*(1-(1+k_result*d*m*s)**(-1/d-2)*p1_result)+(d*m/s)*g*(1-(1+k_result*d*m*s)**(-1/d-1)*p1_result))
                    s4 = s3.mean(axis=1)
                    d1 = []
                    for i1 in s4:
                        d1.append([i1])
                    d1 = np.array(d1)
                    d1 = d1-m1**2
                    if p0 > 0:
                        s5 = s1.mean(axis=1)
                        ls = []
                        for i1 in s5:
                            ls.append([i1])
                        p1 = np.array(ls)

                    if p0 > 0:

                        if (abs(p0-p1)/p0)[np.argmin((abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0))] < 1:
                            l = (abs(p0-p1)/p0)**0.2+(abs(m0-m1)/m0)**0.5
                        else:
                            l = (abs(p0-p1)/p0)**4+(abs(m0-m1)/m0)**0.5

                    else:
                        l = (abs(m0-m1)/m0)**0.5+(abs(v0-d1)/v0)
                    r = np.argmin(l)
                    m_result = l1[r][0]
                    d_result = l2[r][0]




                lme[i_1].append(m_result)
                ldis[i_1].append(d_result)
                lsk_a[i_1].append(k_result)
                lsk_b[i_1].append(p1_result)
            else:
                lme[i_1].append(0)
                ldis[i_1].append(0)
                lsk_a[i_1].append(0)
                lsk_b[i_1].append(0)

    ls4 = lme
    ls12c = ldis

    progress(100,width=50)

    dictx = []
    correlation = []
    l_c = []
    for i1 in range(len(ls4_a)):
        dictx1 = {}

        l2 = []
        l3 = []

        l2 = list(reversed(list(np.argsort(lme[i1]))))
        y = len(l2)
        for j in range(len(l2)):
            if lme[i1][l2[j]] == 0:
                y = j
                break

        lc = []
        if mod != '-f nocopula':
            if y > copula_number:
                for i in range(copula_number):
                    if max(ls4_a[i1][l2[i]]) > 0:
                        l3.append(ls4_a[i1][l2[i]])
                        lc.append(l2[i])
                for i in range(len(lc)):
                    dictx1[lc[i]] = i
                dictx.append(dictx1)
            else:
                for i in range(y):
                    if max(ls4_a[i1][l2[i]]) > 0:
                        l3.append(ls4_a[i1][l2[i]])
                        lc.append(l2[i])
                for i in range(len(lc)):
                    dictx1[lc[i]] = i
                dictx.append(dictx1)
            Mat = np.array(l3)
            correlation1 = np.corrcoef(Mat)
            correlation.append(correlation1)
            #sns.set()
            #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
            #plt.show()

            Mat = np.array(l3[:50])
            correlation1 = np.corrcoef(Mat)
            #sns.set()
            #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
            #plt.show()
        l_c.append(lc)


def creatcount_NBzero_noparam_t():
    """
    Simulate count of -f or -f nocopula

    Parameters
    ----------
    ls4:list
          Mean of gene
    ls12c:list
          Dispersion of gene
    ls7:list
          Size factor

    Returns
    -------
    lcount:list
          Count

    """
    from scipy.stats import nbinom
    from scipy.stats import norm
    global lcount,ls4,correlation,ls12c,ls7_a,dictx,ltruecount,ldp,k,l_c,k_a,ls7_1,ls5a
    def progress(percent,width=50):
        if percent >= 100:
            percent=100

        show_str=('[%%-%ds]' %width) %(int(width * percent/100)*"#")
        print('\r%s %d%%' %(show_str,percent),end='')


    ls5a = []
    ls7_1 = []
    ls1 = [len(i) for i in ls7_a]
    if cell_number == sum(ls1):
        ls7_1 = ls7_a
    else:
        if cell_number < sum(ls1):
            for i in range(len(ls7_a)):
                if i < len(ls7_a)-1:
                    ls7_1.append(ls7_a[i][:int(len(ls7_a[i])*cell_number/sum(ls1))])
                    ls5a.append(int(len(ls7_a[i])*cell_number/sum(ls1)))
                else:
                    ls1a = [len(j) for j in ls7_1]
                    k = sum(ls1a)
                    k1 = cell_number-k
                    ls7_1.append(ls7_a[i][:k1])
                    ls5a.append(k1)
        else:
            for i in range(len(ls7_a)):
                ls7_2 = np.log(ls7_a[i])
                x_1,x_2 = norm.fit(ls7_2)
                if i < len(ls7_a)-1:
                    lsa = [float(j) for j in ls7_a[i]]
                    while 1:
                        p = float(np.exp(np.random.normal(x_1,x_2,1)))
                        if len(lsa) == int(len(ls7_a[i])*cell_number/sum(ls1)):
                            break
                        lsa.append(p)
                    ls7_1.append(lsa)
                    ls5a.append(len(lsa))
                else:
                    ls1a = [len(j) for j in ls7_1]
                    k = sum(ls1a)
                    k1 = cell_number-k
                    lsa = [j for j in ls7_a[i]]
                    if len(lsa) <= k1:
                        while 1:
                            p = float(np.exp(np.random.normal(x_1,x_2,1)))
                            if len(lsa) == k1:
                                break
                            lsa.append(p)
                    else:
                        lsa = lsa[:k1]
                    ls7_1.append(lsa)
                    ls5a.append(len(lsa))

    if param_library != 1:
        ls7_3 = [i for i in ls7_1]
        ls7_1 = []
        for i in range(len(ls7_3)):
            l = [j*param_library for j in ls7_3[i]]
            ls7_1.append(l)


    N = sum(ls1)


    lcount = []
    ltruecount = []
    ldp = []
    print(' ')
    print('Start simulation...')
    def randomm(m):
        k = random.random()
        if k < m:
            return 0
        else:
            return 1
    def integral(m,d,p,k,k_a):
        s = np.array(p)
        return list(np.exp((-1/d)*np.log(m*d)+(-1/d-s)*np.log(1/(d*m)+1)+gammaln(1/d+s)-gammaln(1/d)-gammaln(1+s))*(1-np.exp(-k_a)*(1+k/(1/(d*m)+1))**(-1/d-s)))

    if mod != '-f nocopula':
        k2 = []
        for i in range(len(correlation)):
            l5 = []
            for j in range(len(correlation[i])):
                l5.append(0)
            p = np.random.multivariate_normal(mean=l5,cov=correlation[i],size=len(ls7_1[i]))
            k2.append(list(np.transpose(p)))






    lm = []
    ltruecount = []
    ldp = []
    for i in range(len(ls4)):
        lm.append([])
        ltruecount.append([])
        ldp.append([])
        lcount.append([])
    k3 = (len(ls4[0]))//100
    for i1 in range(len(ls4[0])):
        if i1%k3 == 0:
            progress(int(i1/k3),width=50)
        if cell_number == N and param_library == 1:
            for j1 in range(len(ls4)):
                l = []
                if ls4[j1][i1] == 0:
                    for j in range(len(ls4_a[j1][0])):
                        l.append(0)
                    if mod == '-f nocopula':
                        l2 = []
                        l3 = []
                        for i in range(len(l)):
                            l2.append(int(0))
                            l3.append('uncertain')
                        ltruecount[j1].append(l2)
                        ldp[j1].append(l3)
                else:
                    if i1 in l_c[j1]:
                        l4 = []
                        l5 = []
                        l4 = k2[j1][dictx[j1][i1]]
                        l5 = list(norm.cdf(l4))
                        for j in range(len(ls4_a[j1][0])):
                            j2 = 0
                            p1 = float((1+ls12c[j1][i1]*ls4[j1][i1]*ls7_a[j1][j])**(-1/ls12c[j1][i1])*(1-(lsk_a[j1][i1]/(1/(ls12c[j1][i1]*ls4[j1][i1]*ls7_a[j1][j])+1)+1)**(-1/ls12c[j1][i1]))+(1+lsk_a[j1][i1]*ls12c[j1][i1]*ls4[j1][i1]*ls7_a[j1][j])**(-1/ls12c[j1][i1]))
                            j3 = p1
                            if l5[j] <= p1:
                                i4 = 0
                            else:
                                while 1:
                                    ls = []
                                    for j4 in range(500):
                                        ls.append(int(j4+1+500*j2))
                                    j2 += 1
                                    ls1 = integral(ls4[j1][i1]*ls7_a[j1][j],ls12c[j1][i1],ls,lsk_a[j1][i1],-np.log(lsk_b[j1][i1]))
                                    for j4 in range(len(ls1)):
                                        j5 = j3
                                        j3 += ls1[j4]
                                        if j3 > l5[j]  or (j3 > 0.99 and j5 == j3):
                                            i4 = ls[j4]
                                            break
                                    else:
                                        continue
                                    break
                            l.append(i4)

                        l6 = []
                        for i3 in range(len(l)):
                            l6.append(l[i3]/ls7_a[j1][i3])
                        if max(l6) > 0:
                            lm[j1].append(l6)



                    else:
                        ls1 = np.random.gamma(1/ls12c[j1][i1],ls4[j1][i1]*ls12c[j1][i1],len(ls4_a[j1][0]))
                        l1 = []
                        for j in range(len(ls1)):
                            l1.append(float(ls1[j]*ls7_a[j1][j]))
                        c = 0
                        while 1:
                            c += 1
                            x = np.random.poisson(tuple(l1),(1,len(ls4_a[j1][0])))
                            ls1 = list(map(int,x[0]))
                            if max(ls1) > 0 or c > 10:
                                break
                        l2 = ls1

                        if mod == '-f nocopula':
                            c = 0
                            while 1:
                                c += 1
                                l = []
                                l3 = []
                                for i2 in range(len(l1)):
                                    m = randomm(lsk_b[j1][i1]*np.exp(-lsk_a[j1][i1]*l1[i2]))
                                    l.append(int(l2[i2]*m))
                                    if m == 0:
                                        l3.append('true')
                                    else:
                                        l3.append('false')
                                if max(l) > 0 or c >10:
                                    break
                            ltruecount[j1].append(l2)
                            ldp[j1].append(l3)
                        else:
                            c = 0
                            while 1:
                                c += 1
                                l = []
                                for i2 in range(len(l1)):
                                    m = randomm(lsk_b[j1][i1]*np.exp(-lsk_a[j1][i1]*l1[i2]))
                                    l.append(int(l2[i2]*m))
                                if max(l) > 0 or c > 10:
                                    break
                lcount[j1].append(l)
        else:
            for j1 in range(len(ls4)):
                l = []
                if ls4[j1][i1] == 0:
                    for j in range(len(ls7_1[j1])):
                        l.append(0)
                    if mod == '-f nocopula':
                        l2 = []
                        l3 = []
                        for i in range(len(l)):
                            l2.append(int(0))
                            l3.append('uncertain')
                        ltruecount[j1].append(l2)
                        ldp[j1].append(l3)
                else:
                    if i1 in l_c[j1]:
                        l4 = []
                        l5 = []
                        l4 = k2[j1][dictx[j1][i1]]
                        l5 = list(norm.cdf(l4))
                        for j in range(len(ls7_1[j1])):
                            j2 = 0
                            p1 = float((1+ls12c[j1][i1]*ls4[j1][i1]*ls7_1[j1][j])**(-1/ls12c[j1][i1])*(1-(lsk_a[j1][i1]/(1/(ls12c[j1][i1]*ls4[j1][i1]*ls7_1[j1][j])+1)+1)**(-1/ls12c[j1][i1]))+(1+lsk_a[j1][i1]*ls12c[j1][i1]*ls4[j1][i1]*ls7_1[j1][j])**(-1/ls12c[j1][i1]))
                            j3 = p1
                            if l5[j] <= p1:
                                i4 = 0
                            else:
                                while 1:
                                    ls = []
                                    for j4 in range(500):
                                        ls.append(int(j4+1+500*j2))
                                    j2 += 1
                                    ls1 = integral(ls4[j1][i1]*ls7_1[j1][j],ls12c[j1][i1],ls,lsk_a[j1][i1],-np.log(lsk_b[j1][i1]))
                                    for j4 in range(len(ls1)):
                                        j5 = j3
                                        j3 += ls1[j4]
                                        if j3 > l5[j]  or (j3 > 0.99 and j5 == j3):
                                            i4 = ls[j4]
                                            break
                                    else:
                                        continue
                                    break
                            l.append(i4)

                        l6 = []
                        for i3 in range(len(l)):
                            l6.append(l[i3]/ls7_1[j1][i3])
                        if max(l6) > 0:
                            lm[j1].append(l6)



                    else:
                        ls1 = np.random.gamma(1/ls12c[j1][i1],ls4[j1][i1]*ls12c[j1][i1],len(ls7_1[j1]))
                        l1 = []
                        for j in range(len(ls1)):
                            l1.append(float(ls1[j]*ls7_1[j1][j]))
                        c = 0
                        while 1:
                            c += 1
                            x = np.random.poisson(tuple(l1),(1,len(ls7_1[j1])))
                            ls1 = list(map(int,x[0]))
                            if max(ls1) > 0 or c > 10:
                                break
                        l2 = ls1

                        if mod == '-f nocopula':
                            c = 0
                            while 1:
                                c += 1
                                l = []
                                l3 = []
                                for i2 in range(len(l1)):
                                    m = randomm(lsk_b[j1][i1]*np.exp(-lsk_a[j1][i1]*l1[i2]))
                                    l.append(int(l2[i2]*m))
                                    if m == 0:
                                        l3.append('true')
                                    else:
                                        l3.append('false')
                                if max(l) > 0 or c >10:
                                    break
                            ltruecount[j1].append(l2)
                            ldp[j1].append(l3)
                        else:
                            c = 0
                            while 1:
                                c += 1
                                l = []
                                for i2 in range(len(l1)):
                                    m = randomm(lsk_b[j1][i1]*np.exp(-lsk_a[j1][i1]*l1[i2]))
                                    l.append(int(l2[i2]*m))
                                if max(l) > 0 or c > 10:
                                    break
                lcount[j1].append(l)

    progress(100,width=50)
    print(' ')
    print('Simulation completed')
    if mod != '-f nocopula':
        for i in range(len(lm)):
            Mat = np.array(lm[i])
            correlation1 = np.corrcoef(Mat)

            #sns.set()
            #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
            #plt.show()

            ln = []
            for i1 in range(len(lm[i])):
                ln.append(float(np.mean(lm[i][i1])))
            ln1 = []
            ln1 = list(reversed(list(np.argsort(ln))))
            ln2 = []
            for i1 in range(len(lm[i])):
                if i1 in ln1[:50]:
                    ln2.append(lm[i][i1])
            Mat = np.array(ln2)
            correlation1 = np.corrcoef(Mat)

            #sns.set()
            #ax = sns.clustermap(correlation1,metric = 'euclidean',vmax=1,vmin=-1,cmap='bwr')
            #plt.show()





            
def write_t():
    """
    Write file

    Parameters
    ----------
    dp:list
          Simulated dropouted gene
    ldropout:list
          Dropout information
    llisize:list
          Library size
    lcount:list
          Count

    """
    print('Start writing...')
    global ls2,ls2_a,ls3_a,ls7_a,ls4,ls12c,lcount,ldp,ltruecount,lsk_a,lsk_b,ls7_1,ls5a
    ls = []
    for i in range(len(ls3_a)):
        ls.append([])
    for i in range(len(ls2_a)):
        for j in range(len(ls3_a)):
            if ls3_a[j] == ls2_a[i]:
                ls[j].append(i)
    ls1 = [len(j) for j in ls7_a]
    N = sum(ls1)


    f  = open(path3+'/Simulation information.txt','w')
    if mod == '-e':
        f.write('Method:SimCH-copula-NB(multiple groups)')
        f.write('\n')
        f.write('Data:'+path1+' '+path2)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('Copula gene number'+str(copula_number))
        f.write('\n')
        f.write('Subgroup:'+str(Subgroup))
        f.write('\n')
        f.write('-1\tCell number:'+str(cell_number))
        f.write('\n')
        f.write('-2\tLibrary magnification:'+str(param_library))
        f.write('\n')
        f.close()
    elif mod == '-e nocopula':
        f.write('Method:SimCH-fit-NB(multiple groups)')
        f.write('\n')
        f.write('Data:'+path1+' '+path2)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('Subgroup:'+str(Subgroup))
        f.write('\n')
        f.write('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
        f.write('\n')
        f.write('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')
        f.write('\n')
        f.close()
    elif mod == '-f':
        f.write('Method:SimCH-copula-NBZI(multiple groups)')
        f.write('\n')
        f.write('Data:'+path1+' '+path2)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('Copula gene number'+str(copula_number))
        f.write('\n')
        f.write('Subgroup:'+str(Subgroup))
        f.write('\n')
        f.write('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
        f.write('\n')
        f.write('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')
        f.write('\n')
        f.close()
    elif mod == '-f nocopula':
        f.write('Method:SimCH-fit-NBZI(multiple groups)')
        f.write('\n')
        f.write('Data:'+path1+' '+path2)
        f.write('\n')
        f.write('Seed1:'+str(seed_1))
        f.write('\n')
        f.write('Seed2:'+str(seed_2))
        f.write('\n')
        f.write('Subgroup:'+str(Subgroup))
        f.write('\n')
        f.write('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
        f.write('\n')
        f.write('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')
        f.write('\n')
        f.close()





    if cell_number == N:
        f = open(path3+'/Counts.txt','w')
        for i in range(len(ls)):
            for j in range(len(ls[i])):
                f.write('\t')
                f.write(str(ls2[0][ls[i][j]+1]))
        f.write('\n')
    else:
        f = open(path3+'/Counts.txt','w')
        k = 1
        for i in range(len(ls7_1)):
            for j in range(len(ls7_1[i])):
                f.write('\t')
                f.write('cell'+str(k))
                k+=1
        f.write('\n')

    for i in range(len(lcount[0])):
        f.write(ls2[i+1][0])
        l = []
        for j in range(len(lcount)):
            l.extend(lcount[j][i])
        for j in range(len(l)):
            f.write('\t')
            f.write(str(int(l[j])))
        f.write('\n')
    f.close()

    f = open(path3+'/Gene mean.txt','w')
    for i in range(len(ls3_a)):
        f.write('\t')
        f.write(str(ls3_a[i]))
    f.write('\n')
    for i in range(len(ls4[0])):
        f.write(str(ls2[i+1][0]))

        l = []
        for j in range(len(ls4)):
            l.append(ls4[j][i])
        for j in range(len(l)):
            f.write('\t')
            f.write(str(l[j]))
        f.write('\n')
    f.close()
    if mod == '-f' or mod == '-f nocopula':
        f = open(path3+'/γ.txt','w')
        for i in range(len(ls3_a)):
            f.write('\t')
            f.write(str(ls3_a[i]))
        f.write('\n')
        for i in range(len(ls4[0])):
            f.write(str(ls2[i+1][0]))
            for j in range(len(lsk_b)):
                f.write('\t')
                f.write(str(lsk_b[j][i]))
            f.write('\n')
        f.close()

        f = open(path3+'/ε.txt','w')
        for i in range(len(ls3_a)):
            f.write('\t')
            f.write(str(ls3_a[i]))
        f.write('\n')
        for i in range(len(ls4[0])):
            f.write(str(ls2[i+1][0]))
            for j in range(len(lsk_a)):
                f.write('\t')
                f.write(str(lsk_a[j][i]))
            f.write('\n')
        f.close()

    f = open(path3+'/Gene dispersion.txt','w')
    for i in range(len(ls3_a)):
        f.write('\t')
        f.write(str(ls3_a[i]))
    f.write('\n')
    for i in range(len(ls12c[0])):
        f.write(str(ls2[i+1][0]))

        l = []
        for j in range(len(ls12c)):
            l.append(ls12c[j][i])
        for j in range(len(l)):
            f.write('\t')
            f.write(str(l[j]))
        f.write('\n')
    f.close()

    f = open(path3+'/Cell type.txt','w')
    if cell_number == N:
        for i in range(len(ls)):
            for j in range(len(ls[i])):
                f.write('\t')
                f.write(str(ls2[0][ls[i][j]+1]))
        f.write('\n')
        for i in range(len(ls)):
            for j in range(len(ls[i])):
                f.write('\t')
                f.write(str(ls3_a[i]))
        f.close()
    else:
        k = 1
        for i in range(len(ls7_1)):
            for j in range(len(ls7_1[i])):
                f.write('\t')
                f.write('cell'+str(k))
                k+=1
        f.write('\n')
        for i in range(len(ls5a)):
            for j in range(ls5a[i]):
                f.write('\t')
                f.write(str(ls3_a[i]))
        f.close()


    f = open(path3+'/Library size.txt','w')
    if cell_number == N:
        for i in range(len(ls)):
            for j in range(len(ls[i])):
                f.write('\t')
                f.write(str(ls2[0][ls[i][j]+1]))
        f.write('\n')
        for i in range(len(lcount)):
            for j in range(len(lcount[i][0])):
                l = 0
                for j1 in range(len(lcount[i])):
                    l += int(lcount[i][j1][j])
                f.write('\t')
                f.write(str(l))
        f.close()
    else:
        k = 1
        for i in range(len(ls7_1)):
            for j in range(len(ls7_1[i])):
                f.write('\t')
                f.write('cell'+str(k))
                k+=1
        f.write('\n')
        for i in range(len(lcount)):
            for j in range(len(lcount[i][0])):
                l = 0
                for j1 in range(len(lcount[i])):
                    l += int(lcount[i][j1][j])
                f.write('\t')
                f.write(str(l))
        f.close()




    if mod == '-f nocopula':
        f = open(path3+'/TrueCounts.txt','w')
        for i in range(len(ls2[0])):
            f.write('\t')
            f.write(str(ls2[0][i]))
        f.write('\n')
        for i in range(len(lcount[0])):
            f.write(ls2[i+1][0])
            l = []
            for j in range(len(lcount)):
                l.extend(ltruecount[j][i])
            for j in range(len(l)):
                f.write('\t')
                f.write(str(int(l[j])))
            f.write('\n')
        f.close()



        f = open(path3+'/Dropout.txt','w')
        for i in range(len(ls2[0])):
            f.write('\t')
            f.write(str(ls2[0][i]))
        f.write('\n')

        for i in range(len(lcount[0])):
            f.write(ls2[i+1][0])
            l = []
            for j in range(len(lcount)):
                l.extend(ldp[j][i])
            for j in range(len(l)):
                f.write('\t')
                f.write(str(l[j]))
            f.write('\n')
        f.close()




    #lmean = []
    #lmean_t = []
    #lmean_z = []
    number_t = []
    lbcv = []
    ldp = []
    lpoi = []
    lisizef = []
    ldropout = []
    llisize = []
    




        
import time       
import os
import numpy as np
import random
import math
from scipy.optimize import curve_fit
from scipy.special import gamma
import scipy.stats as stats
from scipy.stats import f_oneway
import scipy
#import tkinter as tk
#from tkinter import filedialog
import matplotlib.pyplot as plt
from random import shuffle
import pandas as pd
#import seaborn as sns
from sympy import *
from scipy.special import gamma
from scipy.special import gammaln
from scipy.special import digamma
import warnings
from scipy.optimize import fsolve
import csv
import scanpy as sc
random.seed(960329)
np.random.seed(990623)
param = []
seed_1 = 960329
seed_2 = 990623
while 1:
    print("\nSimCH is a simulation tool for evaluating scRNA-seq computational methods.\n")
    print("Version:\tv0.1")
    print("Authors:\tWang Gongming, Sun Lei")
    print("Contacts:\tsunlei(at)yzu.edu.cn")
    print("License:\tGNU GENERAL PUBLIC LICENSE -- Version 3")
    print("Copyright (C) 2022 Yangzhou University\n")
    print("It provides simulation modes as follows:\n")
    print('1\tsimulation based on SimCH-flex-NB'+'\n'+'2\tsimulation based on SimCH-flex-NBZI'+'\n'\
      +'3\tsimulation based on SimCH-fit-NB'+'\n'+'4\tsimulation based on SimCH-fit-NBZI' +'\n' \
      +'5\tsimulation based on SimCH-copula-NB'+'\n'+'6\tsimulation based on SimCH-copula-NBZI'+'\n'\
      +'7\tsimulation of independent multiple groups based on SimCH-copula-NB'+'\n'+'8\tsimulation of independent multiple groups based on SimCH-fit-NB'\
      +'\n'+'9\tsimulation of independent multiple groups based on SimCH-copula-NBZI'+'\n'+'10\tsimulation of independent multiple groups based on SimCH-fit-NBZI'+'\n'\
          +'11\tsimulation based on SimCH-copula-NB for evaluating imputation methods\n\tin terms of gene co-expression preservation'+'\n'+'-q\tquit the program\n')

    while(1):
        input_mod = input('Please input option number: ')
        if input_mod == "1" :
            mod = "-a"
            break;
        elif input_mod == "2" :
            mod = "-b"
            break;
        elif input_mod == "5" :
            mod = "-c"
            break;
        elif input_mod == "3" :
            mod = "-c nocopula"
            break;
        elif input_mod == "6" :
            mod = "-d"
            break;
        elif input_mod == "4" :
            mod = "-d nocopula"
            break;
        elif input_mod == "7" :
            mod = "-e"
            break;
        elif input_mod == "8" :
            mod = "-e nocopula"
            break;
        elif input_mod == "" :
            mod = "9"
            break;
        elif input_mod == "10" :
            mod = "-f nocopula"
            break;
        elif input_mod == "11" :
            mod = "-g"
            break;
        elif input_mod == " " :
            mod = " "
            break;
        elif input_mod == "-q":
            import sys
            sys.exit(0);
        else:
            print("Invalid input!!!\n")

    if mod == '-a' or mod == '-b' or mod == '-c' or mod == '-c nocopula' or mod == '-d' or mod == '-d nocopula' or mod == '-e' or mod == '-e nocopula' or mod == '-f' or '-f nocopula' or ' ' or '-g':

        if mod == '-a':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package: ')
                    seed2 = input('Please enter new seed2 for running the numpy package: ')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')
            print('\n'+'Estimating parameters for SimCH-flex-NB simulation...')
            Remove0()
            nml()
            sizef()
            mean()

            time_start=time.time()

            fitmean()
            BCV()
            time_end=time.time()
            #print('totally cost',time_end-time_start)

            param.append([1])
            param.append([0,0.2])
            param.append([1])
            param.append(0.2)
            param.append([0,0.8])
            param.append(0)
            param.append(0.01)
            param.append(0.01)
            param.append(1)

            print('-1     '+'Number of genes:{}'.format(param[0]))
            print('-2     '+'Number of non-zero genes:{}'.format(param[1]))
            print('-3     '+'Number of cells:{}'.format(param[2]))
            print('-4     '+'Ln(Sj) GMM3:{}'.format(param[3]))
            print('-5     '+'Ln(mean) GMM5:{}'.format(param[6][2:]))
            print('-6     '+'Relationship parameter of mean and dispersion:{} '.format(param[19:23]))
            print('-7     '+'Degree of freedom:{}'.format(param[23]))
            print('-8     '+'Group ratio (for extended simulation):{}'.format(param[27]))
            print('-9     '+'Batch ratio (for extended simulation):{}'.format(param[25]))
            print('-10    '+'Batch variation (for extended simulation):{}'.format(param[26]))
            print('-11    '+'Number of path (for extended simulation):{}'.format(param[30]))
            print('-12    '+'DEgene ratio (for extended simulation):{}'.format(param[28]))
            print('-13    '+'DEgene variation (for extended simulation):{}'.format(param[29]))
            print('-14    '+'Non-linear gene ratio (for extended simulation):{}'.format(param[31]))
            print('-15    '+'Marker gene ratio (for extended simulation):{}'.format(param[32]))
            print('-16    '+'Library magnification (for extended simulation):{}'.format(param[33]))

            print('\n')
            print('If you want to modify the parameters above, please enter a number (e.g. -2)')
            print('Otherwise continue to run')
            print('-v means view current parameters')
            print('-s means start simulation')

            while 1:
                mod1 = input('Please input: ')
                if mod1 == '-s' or mod1 == '-v' or mod1 == '-1' or mod1 == '-2' or mod1 == '-3' or mod1 == '-4' or mod1 == '-5' or mod1 == '-6' or mod1 == '-7' or mod1 == '-8' or mod1 == '-9' or mod1 == '-10' or mod1 == '-11' or mod1 == '-12' or mod1 == '-13' or mod1 == '14':

                    if mod1 == '-s':
                        while 1:
                            path2 = input('Please enter a path for saving output files: ')
                            if os.path.exists(path2) == True:
                                break
                            else:
                                print('Path does not exist')
                        time_start=time.time()
                        print('\n'+'Start simulation...')
                        createmean()
                        print('Mean done！')
                        creatBCV()
                        print('BCV done！')
                        creatsizef()
                        print('Size factor done！')
                        creatpoi()
                        print('Poisson done！')
                        time_end=time.time()
                        print('Total cost', time_end-time_start)
                        write()
                        print('Done!\n')
                        print('-m means simulate on this real data again, otherwise go back to the top menu')
                        a = input('Please input: ')
                        if a != '-m':
                            break
                        else:
                            print('-1     '+'Number of genes:{}'.format(param[0]))
                            print('-2     '+'Number of non-zero genes:{}'.format(param[1]))
                            print('-3     '+'Number of cells:{}'.format(param[2]))
                            print('-4     '+'Ln(Sj) GMM3:{}'.format(param[3]))
                            print('-5     '+'Ln(mean) GMM5:{}'.format(param[6][2:]))
                            print('-6     '+'Relationship parameter of mean and dispersion:{} '.format(param[19:23]))
                            print('-7     '+'Degree of freedom:{}'.format(param[23]))
                            print('-8     '+'Group ratio (for extended simulation):{}'.format(param[27]))
                            print('-9     '+'Batch ratio (for extended simulation):{}'.format(param[25]))
                            print('-10    '+'Batch variation (for extended simulation):{}'.format(param[26]))
                            print('-11    '+'Number of path (for extended simulation):{}'.format(param[30]))
                            print('-12    '+'DEgene ratio (for extended simulation):{}'.format(param[28]))
                            print('-13    '+'DEgene variation (for extended simulation):{}'.format(param[29]))
                            print('-14    '+'Non-linear gene ratio (for extended simulation):{}'.format(param[31]))
                            print('-15    '+'Marker gene ratio (for extended simulation):{}'.format(param[32]))
                            print('-16    '+'Library magnification (for extended simulation):{}'.format(param[33]))
                            print('\n')
                            print('If you want to modify the above parameters, please enter a number')
                            print('-v means view current parameters')
                            print('-s means start simulation')
                    elif mod1 == '-v':
                        print('-1     '+'Number of genes:{}'.format(param[0]))
                        print('-2     '+'Number of non-zero genes:{}'.format(param[1]))
                        print('-3     '+'Number of cells:{}'.format(param[2]))
                        print('-4     '+'Ln(Sj) GMM3:{}'.format(param[3]))
                        print('-5     '+'Ln(mean) GMM5:{}'.format(param[6][2:]))
                        print('-6     '+'Relationship parameter of mean and dispersion:{} '.format(param[19:23]))
                        print('-7     '+'Degree of freedom:{}'.format(param[23]))
                        print('-8     '+'Group ratio (for extended simulation):{}'.format(param[27]))
                        print('-9     '+'Batch ratio (for extended simulation):{}'.format(param[25]))
                        print('-10    '+'Batch variation (for extended simulation):{}'.format(param[26]))
                        print('-11    '+'Number of path (for extended simulation):{}'.format(param[30]))
                        print('-12    '+'DEgene ratio (for extended simulation):{}'.format(param[28]))
                        print('-13    '+'DEgene variation (for extended simulation):{}'.format(param[29]))
                        print('-14    '+'Non-linear gene ratio (for extended simulation):{}'.format(param[31]))
                        print('-15    '+'Marker gene ratio (for extended simulation):{}'.format(param[32]))
                        print('-16    '+'Library magnification (for extended simulation):{}'.format(param[33]))
                        print('\n')
                        print('If you want to modify the above parameters, please enter a number')
                        print('-v means view current parameters')
                        print('-s means start simulation')
                    elif mod1 == '-1':
                        mod2 = input('Please input new number of genes: ')
                        try:
                            if int(mod2) > 0:
                                param[0] = int(mod2)
                        except:
                            print('Format error')

                    elif mod1 == '-2':
                        mod2 = input('Please input new number of non-zero genes: ')
                        try:
                            if int(mod2) > 0:
                                param[1] = int(mod2)
                        except:
                            print('Format error')
                    elif mod1 == '-3':
                        mod2 = input('Please input new number of cells: ')
                        try:
                            if int(mod2) > 0:
                                param[2] = int(mod2)
                        except:
                            print('Format error')
                    elif mod1 == '-4':
                        mod2 = input('Please input new ln(Sj) GMM3: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            if sum(b[0:3]) == 1:
                                if b[4] > 0 and b[6] > 0 and b[8] > 0:
                                    param[3] = b
                                else:
                                    print('Standard deviation must be greater than 0')
                            else:
                                print('The sum of the proportions must be 1')
                        except:
                            print('Format error')
                    elif mod1 == '-5':
                        mod2 = input('Please input new ln(mean) GMM5: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            if sum(b[0:5]) == 1:
                                if b[6] > 0 and b[8] > 0 and b[10] > 0 and b[12] > 0 and b[14] > 0:
                                    param[6][2:] = b
                                else:
                                    print('Standard deviation must be greater than 0')
                            else:
                                print('The sum of the proportions must be 1')
                        except:
                            print('Format error')
                    elif mod1 == '-6':
                        mod2 = input('Please input new relationship parameter of mean and dispersion: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            param[19:23] = b
                        except:
                            print('Format error')
                    elif mod1 == '-7':
                        mod2 = input('Please input new degree of freedom: ')
                        try:
                            a = int(mod2)
                            param[23] = a
                        except:
                            print('Format error')
                    elif mod1 == '-8':
                        while 1:
                            a1 = input('Please input new group ratio: ')
                            b = a1[1:-1].split(',')
                            try:
                                c = [float(i) for i in b]
                                if sum(c) == 1:
                                    param_group = c
                                    break
                                else:
                                    print('Format error, the sum of the proportions must be 1')
                            except:
                                c = []
                                d = []
                                for i in range(len(b)):
                                    if '[' in b[i]:
                                        c.append(i)
                                    if  ']' in b[i]:
                                        d.append(i)
                                f = []
                                g = []
                                j = 1
                                if (len(c) + len(d))%2 == 0:
                                    e = a1[1:-1].split(',')
                                    for i in range(len(e)):
                                        if '[' in e[i] or ']' in e[i] or j == -1:
                                            if '[' in e[i]:
                                                j *= -1
                                                try:
                                                    g.append(float(e[i].replace('[','').replace(']','')))
                                                except:
                                                    print('Format error, the input data should be a floating point number')
                                            elif ']' in e[i]:
                                                try:
                                                    g.append(float(e[i].replace(']','').replace('[','')))
                                                except:
                                                    print('Format error, the input data should be a floating point number')
                                                j *= -1
                                                f.append(g)
                                                g = []
                                            else:
                                                try:
                                                    g.append(float(e[i]))
                                                except:
                                                    print('Format error, the input data should be a floating point number')
                                        else:
                                            try:
                                                f.append(float(e[i]))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                    i1 = 0
                                    for i in range(len(f)):
                                        if type(f[i]) == type(0.1):
                                            i1 += f[i]
                                        elif type(f[i]) == type([]):
                                            i1 += sum(f[i])
                                    if i1 == 1:
                                        param[27] = f
                                        break
                                    else:
                                        print('Format error, the sum of the proportions must be 1')
                                else:
                                    print('Format error, the input data should be a floating point number')



                    elif mod1 == '-9':
                        mod2 = input('Please input new batch ratio: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            if sum(b) == 1:
                                param[25] = b
                            else:
                                print('The sum of the proportions must be 1')
                        except:
                            print('Format error')
                    elif mod1 == '-10':
                        mod2 = input('Please input new batch variation: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            param[26] = b
                        except:
                            print('Format error')
                    elif mod1 == '-11':
                        mod2 = input('Please input new number of path: ')
                        try:
                            a = int(mod2)
                            param[30] = a
                        except:
                            print('Format error')
                    elif mod1 == '-12':
                        mod2 = input('Please input new DEgene ratio: ')
                        try:
                            a = float(mod2)
                            param[28] = a
                        except:
                            print('Format error')
                    elif mod1 == '-13':
                        mod2 = input('Please input new DEgene variation: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            param[29] = b
                        except:
                            print('Format error')
                    elif mod1 == '-14':
                        mod2 = input('Please input new Non-linear gene ratio: ')
                        try:
                            a = float(mod2)
                            param[31] = a
                        except:
                            print('Format error')
                    elif mod1 == '-15':
                        mod2 = input('Please input new marker gene ratio: ')
                        try:
                            a = float(mod2)
                            if a < param[28] and a > 0:
                                param[32] = a
                            else:
                                print('Format error')
                        except:
                            print('Format error')
                    elif mod1 == '-16':
                        mod2 = input('Please input new library magnification value: ')
                        try:
                            a = float(mod2)
                            if a > 0:
                                param[33] = a
                            else:
                                print('Format error')
                        except:
                            print('Format error')
                else:
                    print('Format error')
        elif mod == '-b':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package:')
                    seed2 = input('Please enter new seed2 for running the numpy package:')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')
            print('\n'+'Estimating parameters for SimCH-flex-NBZI simulation...')
            Remove0()
            nml()
            sizef()
            time_start=time.time()
            NBzero_mean_dispersion()
            time_end=time.time()
            #print('totally cost',time_end-time_start)
            param.append([1])
            param.append([0,0.2])
            param.append([1])
            param.append(0.2)
            param.append([0,0.8])
            param.append(0)
            param.append(0.01)
            param.append(0.01)
            param.append(1)
            print('-1     '+'Number of genes:{}'.format(param[0]))
            print('-2     '+'Number of non-zero genes:{}'.format(param[1]))
            print('-3     '+'Number of cells:{}'.format(param[2]))
            print('-4     '+'Ln(Sj) GMM3:{}'.format(param[3]))
            print('-5     '+'γ and ε:{} '.format(param[5]))
            print('-6     '+'Ln(mean) GMM5:{}'.format(param[6][2:]))
            print('-7     '+'Relationship parameter of mean and dispersion:{} '.format(param[19:23]))
            print('-8     '+'Degree of freedom:{}'.format(param[23]))
            print('-9     '+'Group ratio (for extended simulation):{}'.format(param[27]))
            print('-10    '+'Batch (for extended simulation):{}'.format(param[25]))
            print('-11    '+'Batch variation (for extended simulation):{}'.format(param[26]))
            print('-12    '+'Number of path (for extended simulation):{}'.format(param[30]))
            print('-13    '+'DEgene ratio (for extended simulation):{}'.format(param[28]))
            print('-14    '+'DEgene variation (for extended simulation):{}'.format(param[29]))
            print('-15    '+'Non-linear gene ratio (for extended simulation):{}'.format(param[31]))
            print('-16    '+'Marker gene ratio (for extended simulation):{}'.format(param[32]))
            print('-17    '+'Library magnification (for extended simulation):{}'.format(param[33]))

            print('\n')
            print('If you want to modify the above parameters, please enter a number')
            print('Otherwise continue to run')
            print('-v means view current parameters')
            print('-s means start simulation')
            while 1:
                mod1 = input('Please input: ')
                if mod1 == '-s' or mod1 == '-v' or mod1 == '-1' or mod1 == '-2' or mod1 == '-3' or mod1 == '-4' or mod1 == '-5' or mod1 == '-6' or mod1 == '-7' or mod1 == '-8' or mod1 == '-9' or mod1 == '-10' or mod1 == '-11' or mod1 == '-12' or mod1 == '-13' or mod1 == '-14' or mod1 == '-15':
                    if mod1 == '-s':
                        while 1:
                            path2 = input('Please enter the saving path: ')
                            if os.path.exists(path2) == True:
                                break
                            else:
                                print('Path does not exist')
                        print('\n'+'Start simulation...')
                        time_start=time.time()
                        createmean()
                        print('Mean done！')
                        creatBCV()
                        print('BCV done！')
                        creatsizef()
                        print('Size factor done！')
                        creatpoi()
                        print('Poisson done！')
                        creatdropout()
                        print('Dropout done！')
                        time_end=time.time()
                        print('Total cost',time_end-time_start)
                        write()
                        print('Done!\n')
                        print('-m means simulate on this data again, otherwise exit the simulation')
                        a = input('Please input: ')
                        if a != '-m':
                            break
                        else:
                            print('-1     '+'Number of genes:{}'.format(param[0]))
                            print('-2     '+'Number of non-zero genes:{}'.format(param[1]))
                            print('-3     '+'Number of cells:{}'.format(param[2]))
                            print('-4     '+'Ln(Sj) GMM3:{}'.format(param[3]))
                            print('-5     '+'P0 and ε:{} '.format(param[5]))
                            print('-6     '+'Ln(mean) GMM5:{}'.format(param[6][2:]))
                            print('-7     '+'Relationship parameter of mean and dispersion:{} '.format(param[19:23]))
                            print('-8     '+'Degree of freedom:{}'.format(param[23]))
                            print('-9     '+'Group ratio (for extended simulation):{}'.format(param[27]))
                            print('-10    '+'Batch (for extended simulation):{}'.format(param[25]))
                            print('-11    '+'Batch variation (for extended simulation):{}'.format(param[26]))
                            print('-12    '+'Number of path (for extended simulation):{}'.format(param[30]))
                            print('-13    '+'DEgene ratio (for extended simulation):{}'.format(param[28]))
                            print('-14    '+'DEgene variation (for extended simulation):{}'.format(param[29]))
                            print('-15    '+'Non-linear gene ratio (for extended simulation):{}'.format(param[31]))
                            print('-16    '+'Marker gene ratio (for extended simulation):{}'.format(param[32]))
                            print('-17    '+'Library magnification (for extended simulation):{}'.format(param[33]))

                            print('\n')
                            print('If you want to modify the above parameters, please enter a number')
                            print('-v means view current parameters')
                            print('-s means start simulation')
                    elif mod1 == '-v':
                        print('-1     '+'Number of genes:{}'.format(param[0]))
                        print('-2     '+'Number of non-zero genes:{}'.format(param[1]))
                        print('-3     '+'Number of cells:{}'.format(param[2]))
                        print('-4     '+'Ln(Sj) GMM3:{}'.format(param[3]))
                        print('-5     '+'P0 and ε:{} '.format(param[5]))
                        print('-6     '+'Ln(mean) GMM5:{}'.format(param[6][2:]))
                        print('-7     '+'Relationship parameter of mean and dispersion:{} '.format(param[19:23]))
                        print('-8     '+'Degree of freedom:{}'.format(param[23]))
                        print('-9     '+'Group ratio (for extended simulation):{}'.format(param[27]))
                        print('-10    '+'Batch (for extended simulation):{}'.format(param[25]))
                        print('-11    '+'Batch variation (for extended simulation):{}'.format(param[26]))
                        print('-12    '+'Number of path (for extended simulation):{}'.format(param[30]))
                        print('-13    '+'DEgene ratio (for extended simulation):{}'.format(param[28]))
                        print('-14    '+'DEgene variation (for extended simulation):{}'.format(param[29]))
                        print('-15    '+'Non-linear gene ratio (for extended simulation):{}'.format(param[31]))
                        print('-16    '+'Marker gene ratio (for extended simulation):{}'.format(param[32]))
                        print('-17    '+'Library magnification (for extended simulation):{}'.format(param[33]))

                        print('\n')
                        print('If you want to modify the above parameters, please enter a number')
                        print('-v means view current parameters')
                        print('-s means start simulation')
                    elif mod1 == '-1':
                        mod2 = input('Please input new number of genes: ')
                        try:
                            if int(mod2) > 0:
                                param[0] = int(mod2)
                        except:
                            print('Format error')

                    elif mod1 == '-2':
                        mod2 = input('Please input new number of non-zero genes: ')
                        try:
                            if int(mod2) > 0:
                                param[1] = int(mod2)
                        except:
                            print('Format error')
                    elif mod1 == '-3':
                        mod2 = input('Please input new number of cells: ')
                        try:
                            if int(mod2) > 0:
                                param[2] = int(mod2)
                        except:
                            print('Format error')
                    elif mod1 == '-4':
                        mod2 = input('Please input new ln(Sj) GMM3: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            if sum(b[0:3]) == 1:
                                if b[4] > 0 and b[6] > 0 and b[8] > 0:
                                    param[3] = b
                                else:
                                    print('Standard deviation must be greater than 0')
                            else:
                                print('The sum of the proportions must be 1')
                        except:
                            print('Format error')
                    elif mod1 == '-5':
                        mod2 = input('Please input new p0 and ε: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            if b[0] > 0 and b[0] < 1 and b[1] > 0 and len(b) == 2:
                                param[3] = b
                            else:
                                    print('0 < p0 < 1 and ε > 0')
                        except:
                            print('Format error')

                    elif mod1 == '-6':
                        mod2 = input('Please input new ln(mean) GMM5: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            if sum(b[0:5]) == 1:
                                if b[6] > 0 and b[8] > 0 and b[10] > 0 and b[12] > 0 and b[14] > 0:
                                    param[6][2:] = b
                                else:
                                    print('Standard deviation must be greater than 0')
                            else:
                                print('The sum of the proportions must be 1')
                        except:
                            print('Format error')
                    elif mod1 == '-7':
                        mod2 = input('Please input new relationship parameter of mean and dispersion: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            param[19:23] = b
                        except:
                            print('Format error')
                    elif mod1 == '-8':
                        mod2 = input('Please input new degree of freedom: ')
                        try:
                            a = int(mod2)
                            param[23] = a
                        except:
                            print('Format error')
                    elif mod1 == '-9':
                        while 1:
                            a1 = input('Please input new group ratio: ')
                            b = a1[1:-1].split(',')
                            try:
                                c = [float(i) for i in b]
                                if sum(c) == 1:
                                    param_group = c
                                    break
                                else:
                                    print('Format error, the sum of the proportions must be 1')
                            except:
                                c = []
                                d = []
                                for i in range(len(b)):
                                    if '[' in b[i]:
                                        c.append(i)
                                    if  ']' in b[i]:
                                        d.append(i)
                                f = []
                                g = []
                                j = 1
                                if (len(c) + len(d))%2 == 0:
                                    e = a1[1:-1].split(',')
                                    for i in range(len(e)):
                                        if '[' in e[i] or ']' in e[i] or j == -1:
                                            if '[' in e[i]:
                                                j *= -1
                                                try:
                                                    g.append(float(e[i].replace('[','').replace(']','')))
                                                except:
                                                    print('Format error, the input data should be a floating point number')
                                            elif ']' in e[i]:
                                                try:
                                                    g.append(float(e[i].replace(']','').replace('[','')))
                                                except:
                                                    print('Format error, the input data should be a floating point number')
                                                j *= -1
                                                f.append(g)
                                                g = []
                                            else:
                                                try:
                                                    g.append(float(e[i]))
                                                except:
                                                    print('Format error, the input data should be a floating point number')
                                        else:
                                            try:
                                                f.append(float(e[i]))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                    i1 = 0
                                    for i in range(len(f)):
                                        if type(f[i]) == type(0.1):
                                            i1 += f[i]
                                        elif type(f[i]) == type([]):
                                            i1 += sum(f[i])
                                    if i1 == 1:
                                        param[27] = f
                                        break
                                    else:
                                        print('Format error, the sum of the proportions must be 1')
                                else:
                                    print('Format error, the input data should be a floating point number')


                    elif mod1 == '-10':
                        mod2 = input('Please input new batch ratio: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            if sum(b) == 1:
                                param[25] = b
                            else:
                                print('The sum of the proportions must be 1')
                        except:
                            print('Format error')
                    elif mod1 == '-11':
                        mod2 = input('Please input new batch variation: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            param[26] = b
                        except:
                            print('Format error')
                    elif mod1 == '-12':
                        mod2 = input('Please input new number of path:')
                        try:
                            a = int(mod2)
                            param[30] = a
                        except:
                            print('Format error')
                    elif mod1 == '-13':
                        mod2 = input('Please input new DEgene ratio: ')
                        try:
                            a = float(mod2)
                            param[28] = a
                        except:
                            print('Format error')
                    elif mod1 == '-14':
                        mod2 = input('Please input new DEgene variation: ')
                        try:
                            a = mod2[1:-1].strip().split(',')
                            b = [float(i) for i in a]
                            param[29] = b
                        except:
                            print('Format error')
                    elif mod1 == '-15':
                        mod2 = input('Please input new non-linear gene ratio: ')
                        try:
                            a = float(mod2)
                            param[31] = a
                        except:
                            print('Format error')
                    elif mod1 == '-16':
                        mod2 = input('Please input new marker gene ratio: ')
                        try:
                            a = float(mod2)
                            if a < param[28] and a > 0:
                                param[32] = a
                            else:
                                print('Format error')
                        except:
                            print('Format error')
                    elif mod1 == '-17':
                        mod2 = input('Please input new library magnification value: ')
                        try:
                            a = float(mod2)
                            if a > 0:
                                param[33] = a
                            else:
                                print('Format error')
                        except:
                            print('Format error')

                else:
                    print('Format error')

        elif mod == '-c':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package: ')
                    seed2 = input('Please enter new seed2 for running the numpy package: ')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')


            copula_number = 2000

            print('-1\tCopula genes number:'+str(copula_number)+int(20-len(str(copula_number)))*' ')
            print('If you want to modify the above parameters, please enter a number')
            print('Enter other to start the simulation')
            a = input('Please input: ')
            if a == '-1':
                while 1:
                    try:
                        a1 = input('Please input new copula genes number:')
                        b = int(a1)
                        if b > 0 :
                            copula_number = b
                            break
                        else:
                            print('Copula genes number should be greater than 0')
                    except:
                        print('Format error, the input data should be a positive integer')
            print('\n'+'Estimating parameters for SimCH-copula-NB simulation...')
            Remove0_noparam()
            nml_noparam()
            time_start=time.time()
            mean_noparam()
            BCV_noparam()
            time_end=time.time()
            #print('totally cost',time_end-time_start)

            param_group = [1]
            param_batch = [1]
            param_batch_var = [0,0.2]
            param_path = 'No'
            param_DEgene = 0.2
            param_DEgene_var = [0,0.5]
            param_noline = 0.01
            param_marker = 0.01
            param_library = 1
            cell_number = len(ls7)

            print('\n')
            print('-1\tGroup ratio (for extended simulation):'+str(param_group)+int(28-len(str(param_group)))*' ')
            print('-2\tBatch ratio (for extended simulation):'+str(param_batch)+int(28-len(str(param_batch)))*' ')
            print('-3\tBatch variation (for extended simulation):'+str(param_batch_var)+int(24-len(str(param_batch_var)))*' ')
            print('-4\tPath (for extended simulation):'+str(param_path)+int(35-len(str(param_path)))*' ')
            print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene)+int(27-len(str(param_DEgene)))*' ')
            print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var)+int(23-len(str(param_DEgene_var)))*' ')
            print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline)+int(18-len(str(param_noline)))*' ')
            print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker)+int(22-len(str(param_marker)))*' ')
            print('-9\tLibrary magnification (for extended simulation):'+str(param_library)+int(18-len(str(param_library)))*' ')
            print('-10\tCell number (for extended simulation):'+str(cell_number)+int(28-len(str(cell_number)))*' ')
            print('\n')
            print('If you want to modify the above parameters, please enter a number')
            print('-v means view current parameters')
            print('-s means start simulation')

            while 1:
                a = input('Please input: ')
                if a == '-1':
                    while 1:
                        a1 = input('Please input new group ratio: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            if sum(c) == 1:
                                param_group = c
                                break
                            else:
                                print('Format error, the sum of the proportions must be 1')
                        except:
                            c = []
                            d = []
                            for i in range(len(b)):
                                if '[' in b[i]:
                                    c.append(i)
                                if  ']' in b[i]:
                                    d.append(i)
                            f = []
                            g = []
                            j = 1
                            if (len(c) + len(d))%2 == 0:
                                e = a1[1:-1].split(',')
                                for i in range(len(e)):
                                    if '[' in e[i] or ']' in e[i] or j == -1:
                                        if '[' in e[i]:
                                            j *= -1
                                            try:
                                                g.append(float(e[i].replace('[','').replace(']','')))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                        elif ']' in e[i]:
                                            try:
                                                g.append(float(e[i].replace(']','').replace('[','')))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                            j *= -1
                                            f.append(g)
                                            g = []
                                        else:
                                            try:
                                                g.append(float(e[i]))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                    else:
                                        try:
                                            f.append(float(e[i]))
                                        except:
                                            print('Format error, the input data should be a floating point number')
                                i1 = 0
                                for i in range(len(f)):
                                    if type(f[i]) == type(0.1):
                                        i1 += f[i]
                                    elif type(f[i]) == type([]):
                                        i1 += sum(f[i])
                                if i1 == 1:
                                    param_group = f
                                    break
                                else:
                                    print('Format error, the sum of the proportions must be 1')
                            else:
                                print('Format error, the input data should be a floating point number')




                elif a == '-2':
                    while 1:
                        a1 = input('Please input new batch ratio: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            if sum(c) == 1:
                                param_batch = c
                                break
                            else:
                                print('Format error, the sum of the proportions must be 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-3':
                    while 1:
                        a1 = input('Please input new batch variation: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            param_batch_var = c
                            break
                        except:
                            print('Format error, the input data should be a floating point number')


                elif a == '-4':
                    while 1:
                        try:
                            a1 = input('If you wan to input path (Yes or No)? ')
                            if a1 == 'Yes' or a1 == 'No':
                                param_path = a1
                                break
                            else:
                                print('Enter Yes or No')
                        except:
                            print('Format error')
                elif a == '-5':
                    while 1:
                        try:
                            a1 = input('Please input new DEgene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1:
                                param_DEgene = b
                                break
                            else:
                                print('DEgene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-6':
                    while 1:
                        a1 = input('Please input new DEgene variation:')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            param_DEgene_var = c
                            break
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-7':
                    while 1:
                        try:
                            a1 = input('Please input new Non-linear gene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1 and b+param_DEgene <= 1:
                                param_noline = b
                                break
                            else:
                                print('Non-linear gene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')

                elif a == '-8':
                    while 1:
                        try:
                            a1 = input('Please input new marker gene ratio:')
                            b = float(a1)
                            if b >= 0 and b <= 1 and b < param_DEgene:
                                param_marker = b
                                break
                            else:
                                print('Marker gene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-9':
                    while 1:
                        try:
                            a1 = input('Please input new library magnification: ')
                            b = float(a1)
                            if b > 0 :
                                param_library = b
                                break
                            else:
                                print('Library magnification should be greater than 0')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-10':
                    while 1:
                        try:
                            a1 = input('Please input new cell number:')
                            b = int(a1)
                            if b > 0 :
                                cell_number = b
                                break
                            else:
                                print('Cell number should be greater than 0')
                        except:
                            print('Format error, the input data should be a positive integer')
                elif a == '-v':
                    print('\n')
                    print('-1\tGroup ratio (for extended simulation):'+str(param_group)+int(28-len(str(param_group)))*' ')
                    print('-2\tBatch ratio (for extended simulation):'+str(param_batch)+int(28-len(str(param_batch)))*' ')
                    print('-3\tBatch variation (for extended simulation):'+str(param_batch_var)+int(24-len(str(param_batch_var)))*' ')
                    print('-4\tPath (for extended simulation):'+str(param_path)+int(35-len(str(param_path)))*' ')
                    print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene)+int(27-len(str(param_DEgene)))*' ')
                    print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var)+int(23-len(str(param_DEgene_var)))*' ')
                    print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline)+int(18-len(str(param_noline)))*' ')
                    print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker)+int(22-len(str(param_marker)))*' ')
                    print('-9\tLibrary magnification (for extended simulation):'+str(param_library)+int(18-len(str(param_library)))*' ')
                    print('-10\tCell number (for extended simulation):'+str(cell_number)+int(28-len(str(cell_number)))*' ')
                    print('\n')
                elif a == '-s':
                    time_start=time.time()

                    creatcount_noparam()
                    time_end=time.time()
                    #print('totally cost',time_end-time_start)
                    while 1:
                        path3 = input('Please enter the saving path: ')
                        if os.path.exists(path3) == True:
                            break
                        else:
                            print('Path does not exist')
                    write_noparam()
                    print('Done!')
                    print('\n')
                    print('-m means simulate on this data again, otherwise exit the simulation')
                    a1 = input('Please input: ')
                    if a1 != '-m':
                        break
                    else:
                        print('\n')
                        print('-1\tGroup ratio:'+str(param_group))
                        print('-2\tBatch ratio:'+str(param_batch))
                        print('-3\tBatch variation:'+str(param_batch_var))
                        print('-4\tPath:'+str(param_path))
                        print('-5\tDEgene ratio:'+str(param_DEgene))
                        print('-6\tDEgene variation:'+str(param_DEgene_var))
                        print('-7\tNon-linear gene ratio:'+str(param_noline))
                        print('-8\tMarker gene ratio:'+str(param_marker))
                        print('-9\tLibrary magnification:'+str(param_library))
                        print('-10\tCell number:'+str(cell_number))
                        print('If you want to modify the above parameters, please enter a number')
                        print('-v means view current parameters')
                        print('-s means start simulation')
                        print('\n')
                else:
                    print('Format error')

        elif mod == '-c nocopula':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package:')
                    seed2 = input('Please enter new seed2 for running the numpy package:')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')
            copula_number = 2000
            print('-1\tCopula genes number:'+str(copula_number)+int(20-len(str(copula_number)))*' ')
            print('If you want to modify the above parameters, please enter a number')
            print('Enter other to start the simulation')
            a = input('Please input: ')
            if a == '-1':
                while 1:
                    try:
                        a1 = input('Please input new copula genes number:')
                        b = int(a1)
                        if b > 0 :
                            copula_number = b
                            break
                        else:
                            print('Copula genes number should be greater than 0')
                    except:
                        print('Format error, the input data should be a positive integer')
            print('\n'+'Estimating parameters for SimCH-fit-NB simulation...')
            Remove0_noparam()
            nml_noparam()
            time_start=time.time()
            mean_noparam()
            BCV_noparam()
            time_end=time.time()
            #print('totally cost',time_end-time_start)

            param_group = [1]
            param_batch = [1]
            param_batch_var = [0,0.2]
            param_path = 'No'
            param_DEgene = 0.2
            param_DEgene_var = [0,0.5]
            param_noline = 0.01
            param_marker = 0.01
            param_library = 1
            cell_number = len(ls7)

            print('\n')
            print('-1\tGroup ratio (for extended simulation):'+str(param_group)+int(28-len(str(param_group)))*' ')
            print('-2\tBatch ratio (for extended simulation):'+str(param_batch)+int(28-len(str(param_batch)))*' ')
            print('-3\tBatch variation (for extended simulation):'+str(param_batch_var)+int(24-len(str(param_batch_var)))*' ')
            print('-4\tPath (for extended simulation):'+str(param_path)+int(35-len(str(param_path)))*' ')
            print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene)+int(27-len(str(param_DEgene)))*' ')
            print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var)+int(23-len(str(param_DEgene_var)))*' ')
            print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline)+int(18-len(str(param_noline)))*' ')
            print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker)+int(22-len(str(param_marker)))*' ')
            print('-9\tLibrary magnification (for extended simulation):'+str(param_library)+int(18-len(str(param_library)))*' ')
            print('-10\tCell number (for extended simulation):'+str(cell_number)+int(28-len(str(cell_number)))*' ')
            print('\n')
            print('If you want to modify the above parameters, please enter a number')
            print('Otherwise continue to run')
            print('-v means view current parameters')
            print('-s means start simulation')
            while 1:
                a = input('Please input: ')
                if a == '-1':
                    while 1:
                        a1 = input('Please input new group ratio:')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            if sum(c) == 1:
                                param_group = c
                                break
                            else:
                                print('Format error, the sum of the proportions must be 1')
                        except:
                            c = []
                            d = []
                            for i in range(len(b)):
                                if '[' in b[i]:
                                    c.append(i)
                                if  ']' in b[i]:
                                    d.append(i)
                            f = []
                            g = []
                            j = 1
                            if (len(c) + len(d))%2 == 0:
                                e = a1[1:-1].split(',')
                                for i in range(len(e)):
                                    if '[' in e[i] or ']' in e[i] or j == -1:
                                        if '[' in e[i]:
                                            j *= -1
                                            try:
                                                g.append(float(e[i].replace('[','').replace(']','')))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                        elif ']' in e[i]:
                                            try:
                                                g.append(float(e[i].replace(']','').replace('[','')))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                            j *= -1
                                            f.append(g)
                                            g = []
                                        else:
                                            try:
                                                g.append(float(e[i]))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                    else:
                                        try:
                                            f.append(float(e[i]))
                                        except:
                                            print('Format error, the input data should be a floating point number')
                                i1 = 0
                                for i in range(len(f)):
                                    if type(f[i]) == type(0.1):
                                        i1 += f[i]
                                    elif type(f[i]) == type([]):
                                        i1 += sum(f[i])
                                if i1 == 1:
                                    param_group = f
                                    break
                                else:
                                    print('Format error, the sum of the proportions must be 1')
                            else:
                                print('Format error, the input data should be a floating point number')


                elif a == '-2':
                    while 1:
                        a1 = input('Please input new batch ratio: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            if sum(c) == 1:
                                param_batch = c
                                break
                            else:
                                print('Format error, the sum of the proportions must be 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-3':
                    while 1:
                        a1 = input('Please input new batch variation:')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            param_batch_var = c
                            break
                        except:
                            print('Format error, the input data should be a floating point number')


                elif a == '-4':
                    while 1:
                        try:
                            a1 = input('If you want to input path (Yes or No)? ')
                            if a1 == 'Yes' or a1 == 'No':
                                param_path = a1
                                break
                            else:
                                print('Enter Yes or No')
                        except:
                            print('Format error')
                elif a == '-5':
                    while 1:
                        try:
                            a1 = input('Please input new DEgene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1:
                                param_DEgene = b
                                break
                            else:
                                print('DEgene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-6':
                    while 1:
                        a1 = input('Please input new DEgene variation:')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            param_DEgene_var = c
                            break
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-7':
                    while 1:
                        try:
                            a1 = input('Please input new Non-linear gene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1 and b+param_DEgene <= 1:
                                param_DEgene_var = b
                                break
                            else:
                                print('Non-linear gene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-8':
                    while 1:
                        try:
                            a1 = input('Please input new marker gene ratio:')
                            b = float(a1)
                            if b >= 0 and b <= 1 and b < param_DEgene:
                                param_marker = b
                                break
                            else:
                                print('Marker gene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-9':
                    while 1:
                        try:
                            a1 = input('Please input new library magnification: ')
                            b = float(a1)
                            if b > 0 :
                                param_library = b
                                break
                            else:
                                print('Library magnification should be greater than 0')
                        except:
                            print('Format error, the input data should be  a floating point number')
                elif a == '-10':
                    while 1:
                        try:
                            a1 = input('Please input new cell number:')
                            b = int(a1)
                            if b > 0 :
                                cell_number = b
                                break
                            else:
                                print('Cell number should be greater than 0')
                        except:
                            print('Format error, the input data should be a positive integer')
                elif a == '-v':
                    print('-1\tGroup ratio (for extended simulation):'+str(param_group)+int(28-len(str(param_group)))*' ')
                    print('-2\tBatch ratio (for extended simulation):'+str(param_batch)+int(28-len(str(param_batch)))*' ')
                    print('-3\tBatch variation (for extended simulation):'+str(param_batch_var)+int(24-len(str(param_batch_var)))*' ')
                    print('-4\tPath (for extended simulation):'+str(param_path)+int(35-len(str(param_path)))*' ')
                    print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene)+int(27-len(str(param_DEgene)))*' ')
                    print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var)+int(23-len(str(param_DEgene_var)))*' ')
                    print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline)+int(18-len(str(param_noline)))*' ')
                    print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker)+int(22-len(str(param_marker)))*' ')
                    print('-9\tLibrary magnification (for extended simulation):'+str(param_library)+int(18-len(str(param_library)))*' ')
                    print('-10\tCell number (for extended simulation):'+str(cell_number)+int(28-len(str(cell_number)))*' ')
                elif a == '-s':
                    time_start=time.time()
                    creatcount_noparam()
                    time_end=time.time()
                    #print('totally cost',time_end-time_start)

                    while 1:
                        path3 = input('Please enter the saving path: ')
                        if os.path.exists(path3) == True:
                            break
                        else:
                            print('Path does not exist')
                    write_noparam()
                    print('Done!')
                    print('\n')
                    print('-m means simulate on this data again, otherwise exit the simulation')
                    a1 = input('Please input: ')
                    if a1 != '-m':
                        break
                    else:
                        print('-1\tGroup ratio (for extended simulation):'+str(param_group)+int(28-len(str(param_group)))*' ')
                        print('-2\tBatch ratio (for extended simulation):'+str(param_batch)+int(28-len(str(param_batch)))*' ')
                        print('-3\tBatch variation (for extended simulation):'+str(param_batch_var)+int(24-len(str(param_batch_var)))*' ')
                        print('-4\tPath (for extended simulation):'+str(param_path)+int(35-len(str(param_path)))*' ')
                        print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene)+int(27-len(str(param_DEgene)))*' ')
                        print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var)+int(23-len(str(param_DEgene_var)))*' ')
                        print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline)+int(18-len(str(param_noline)))*' ')
                        print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker)+int(22-len(str(param_marker)))*' ')
                        print('-9\tLibrary magnification (for extended simulation):'+str(param_library)+int(18-len(str(param_library)))*' ')
                        print('-10\tCell number (for extended simulation):'+str(cell_number)+int(28-len(str(cell_number)))*' ')
                        print('\n')
                        print('If you want to modify the above parameters, please enter a number')
                        print('-v means view current parameters')
                        print('-s means start simulation')
                else:
                    print('Format error')


        elif mod == '-d':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package: ')
                    seed2 = input('Please enter new seed2 for running the numpy package: ')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')

            copula_number = 2000
            print('-1\tCopula genes number:'+str(copula_number)+int(20-len(str(copula_number)))*' ')
            print('If you want to modify the above parameters, please enter a number')
            print('Enter other to start the simulation')
            a = input('Please input: ')
            if a == '-1':
                while 1:
                    try:
                        a1 = input('Please input new copula genes number: ')
                        b = int(a1)
                        if b > 0 :
                            copula_number = b
                            break
                        else:
                            print('Copula genes number should be greater than 0')
                    except:
                        print('Format error, the input data should be a positive integer')
            print('\n'+'Estimating parameters for SimCH-copula-NBZI simulation...')
            Remove0_noparam()
            nml_noparam()
            time_start=time.time()

            NBzero_mean_dispersion_noparam()
            time_end=time.time()
            #print('totally cost',time_end-time_start)


            param_group = [1]
            param_batch = [1]
            param_batch_var = [0,0.2]
            param_path = 'No'
            param_DEgene = 0.2
            param_DEgene_var = [0,0.5]
            param_noline = 0.01
            param_marker = 0.01
            param_library = 1
            cell_number = len(ls7)
            print('\n')
            print('-1\tGroup ratio (for extended simulation):'+str(param_group))
            print('-2\tBatch ratio (for extended simulation):'+str(param_batch))
            print('-3\tBatch variation (for extended simulation):'+str(param_batch_var))
            print('-4\tPath (for extended simulation):'+str(param_path))
            print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene))
            print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var))
            print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline))
            print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker))
            print('-9\tLibrary magnification (for extended simulation):'+str(param_library))
            print('-10\tCell number (for extended simulation):'+str(cell_number))
            print('\n')
            print('If you want to modify the above parameters, please enter a number')
            print('-v means view current parameters')
            print('-s means start simulation')
            while 1:
                a = input('Please input: ')
                if a == '-1':
                    while 1:
                        a1 = input('Please input new group ratio: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            if sum(c) == 1:
                                param_group = c
                                break
                            else:
                                print('Format error, the sum of the proportions must be 1')
                        except:
                            c = []
                            d = []
                            for i in range(len(b)):
                                if '[' in b[i]:
                                    c.append(i)
                                if  ']' in b[i]:
                                    d.append(i)
                            f = []
                            g = []
                            j = 1
                            if (len(c) + len(d))%2 == 0:
                                e = a1[1:-1].split(',')
                                for i in range(len(e)):
                                    if '[' in e[i] or ']' in e[i] or j == -1:
                                        if '[' in e[i]:
                                            j *= -1
                                            try:
                                                g.append(float(e[i].replace('[','').replace(']','')))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                        elif ']' in e[i]:
                                            try:
                                                g.append(float(e[i].replace(']','').replace('[','')))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                            j *= -1
                                            f.append(g)
                                            g = []
                                        else:
                                            try:
                                                g.append(float(e[i]))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                    else:
                                        try:
                                            f.append(float(e[i]))
                                        except:
                                            print('Format error, the input data should be a floating point number')
                                i1 = 0
                                for i in range(len(f)):
                                    if type(f[i]) == type(0.1):
                                        i1 += f[i]
                                    elif type(f[i]) == type([]):
                                        i1 += sum(f[i])
                                if i1 == 1:
                                    param_group = f
                                    break
                                else:
                                    print('Format error, the sum of the proportions must be 1')
                            else:
                                print('Format error, the input data should be a floating point number')


                elif a == '-2':
                    while 1:
                        a1 = input('Please input new batch ratio: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            if sum(c) == 1:
                                param_batch = c
                                break
                            else:
                                print('Format error, the sum of the proportions must be 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-3':
                    while 1:
                        a1 = input('Please input new batch variation: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            param_batch_var = c
                            break
                        except:
                            print('Format error, the input data should be a floating point number')


                elif a == '-4':
                    while 1:
                        try:
                            a1 = input('If you want to input path (Yes or No)? ')
                            if a1 == 'Yes' or a1 == 'No':
                                param_path = a1
                                break
                            else:
                                print('Enter Yes or No')
                        except:
                            print('Format error')
                elif a == '-5':
                    while 1:
                        try:
                            a1 = input('Please input new DEgene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1:
                                param_DEgene = b
                                break
                            else:
                                print('DEgene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-6':
                    while 1:
                        a1 = input('Please input new DEgene variation: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            param_DEgene_var = c
                            break
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-7':
                    while 1:
                        try:
                            a1 = input('Please input new Non-linear gene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1 and b+param_DEgene <= 1:
                                param_DEgene_var = b
                                break
                            else:
                                print('Non-linear gene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-8':
                    while 1:
                        try:
                            a1 = input('Please input new marker gene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1 and b < param_DEgene:
                                param_marker = b
                                break
                            else:
                                print('Marker gene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-9':
                    while 1:
                        try:
                            a1 = input('Please input new library magnification: ')
                            b = float(a1)
                            if b > 0 :
                                param_library = b
                                break
                            else:
                                print('Library magnification should be greater than 0')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-10':
                    while 1:
                        try:
                            a1 = input('Please input new cell number: ')
                            b = int(a1)
                            if b > 0 :
                                cell_number = b
                                break
                            else:
                                print('Cell number should be greater than 0')
                        except:
                            print('Format error, the input data should be a positive integer')
                elif a == '-v':
                    print('-1\tGroup ratio (for extended simulation):'+str(param_group))
                    print('-2\tBatch ratio (for extended simulation):'+str(param_batch))
                    print('-3\tBatch variation (for extended simulation):'+str(param_batch_var))
                    print('-4\tPath (for extended simulation):'+str(param_path))
                    print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene))
                    print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var))
                    print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline))
                    print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker))
                    print('-9\tLibrary magnification:'+str(param_library))
                    print('-10\tCell number (for extended simulation):'+str(cell_number))
                elif a == '-s':
                    time_start=time.time()
                    creatcount_NBzero_noparam()
                    time_end=time.time()
                    print('Total cost',time_end-time_start)
                    while 1:
                        path3 = input('Please enter the saving path: ')
                        if os.path.exists(path3) == True:
                            break
                        else:
                            print('Path does not exist')
                    write_noparam()
                    print('Done!')
                    print('\n')
                    print('-m means simulate on this data again, otherwise exit the simulation')
                    a1 = input('Please input: ')
                    if a1 != '-m':
                        break
                    else:
                        print('-1\tGroup ratio (for extended simulation):'+str(param_group))
                        print('-2\tBatch ratio (for extended simulation):'+str(param_batch))
                        print('-3\tBatch variation (for extended simulation):'+str(param_batch_var))
                        print('-4\tPath (for extended simulation):'+str(param_path))
                        print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene))
                        print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var))
                        print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline))
                        print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker))
                        print('-9\tLibrary magnification (for extended simulation):'+str(param_library))
                        print('-10\tCell number (for extended simulation):'+str(cell_number))
                        print('\n')
                        print('If you want to modify the above parameters, please enter a number')
                        print('-v means view current parameters')
                        print('-s means start simulation')
                else:
                    print('Format error')



        elif mod =='-d nocopula':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package: ')
                    seed2 = input('Please enter new seed2 for running the numpy package: ')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')

            copula_number = 2000
            print('-1\tCopula genes number:'+str(copula_number)+int(20-len(str(copula_number)))*' ')
            print('If you want to modify the above parameters, please enter a number')
            print('Enter other to start the simulation')
            a = input('Please input: ')
            if a == '-1':
                while 1:
                    try:
                        a1 = input('Please input new copula genes number: ')
                        b = int(a1)
                        if b > 0 :
                            copula_number = b
                            break
                        else:
                            print('Copula genes number should be greater than 0')
                    except:
                        print('Format error, the input data should be a positive integer')
            print('\n'+'Estimating parameters for SimCH-fit-NBZI simulation...')
            Remove0_noparam()
            nml_noparam()
            time_start=time.time()
            NBzero_mean_dispersion_noparam()
            time_end=time.time()
            #print('totally cost',time_end-time_start)

            param_group = [1]
            param_batch = [1]
            param_batch_var = [0,0.2]
            param_path = 'No'
            param_DEgene = 0.2
            param_DEgene_var = [0,0.5]
            param_noline = 0.01
            param_marker = 0.01
            param_library = 1
            cell_number = len(ls7)
            print('\n')
            print('-1\tGroup ratio (for extended simulation):'+str(param_group)+int(28-len(str(param_group)))*' ')
            print('-2\tBatch ratio (for extended simulation):'+str(param_batch)+int(28-len(str(param_batch)))*' ')
            print('-3\tBatch variation (for extended simulation):'+str(param_batch_var)+int(24-len(str(param_batch_var)))*' ')
            print('-4\tPath (for extended simulation):'+str(param_path)+int(35-len(str(param_path)))*' ')
            print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene)+int(27-len(str(param_DEgene)))*' ')
            print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var)+int(23-len(str(param_DEgene_var)))*' ')
            print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline)+int(18-len(str(param_noline)))*' ')
            print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker)+int(22-len(str(param_marker)))*' ')
            print('-9\tLibrary magnification (for extended simulation):'+str(param_library)+int(18-len(str(param_library)))*' ')
            print('-10\tCell number (for extended simulation):'+str(cell_number)+int(28-len(str(cell_number)))*' ')
            print('\n')
            print('If you want to modify the above parameters, please enter a number')
            print('-v means view current parameters')
            print('-s means start simulation')
            while 1:
                a = input('Please input: ')
                if a == '-1':
                    while 1:
                        a1 = input('Please input new group ratio: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            if sum(c) == 1:
                                param_group = c
                                break
                            else:
                                print('Format error, the sum of the proportions must be 1')
                        except:
                            c = []
                            d = []
                            for i in range(len(b)):
                                if '[' in b[i]:
                                    c.append(i)
                                if  ']' in b[i]:
                                    d.append(i)
                            f = []
                            g = []
                            j = 1
                            if (len(c) + len(d))%2 == 0:
                                e = a1[1:-1].split(',')
                                for i in range(len(e)):
                                    if '[' in e[i] or ']' in e[i] or j == -1:
                                        if '[' in e[i]:
                                            j *= -1
                                            try:
                                                g.append(float(e[i].replace('[','').replace(']','')))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                        elif ']' in e[i]:
                                            try:
                                                g.append(float(e[i].replace(']','').replace('[','')))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                            j *= -1
                                            f.append(g)
                                            g = []
                                        else:
                                            try:
                                                g.append(float(e[i]))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                    else:
                                        try:
                                            f.append(float(e[i]))
                                        except:
                                            print('Format error, the input data should be a floating point number')
                                i1 = 0
                                for i in range(len(f)):
                                    if type(f[i]) == type(0.1):
                                        i1 += f[i]
                                    elif type(f[i]) == type([]):
                                        i1 += sum(f[i])
                                if i1 == 1:
                                    param_group = f
                                    break
                                else:
                                    print('Format error, the sum of the proportions must be 1')
                            else:
                                print('Format error, the input data should be a floating point number')


                elif a == '-2':
                    while 1:
                        a1 = input('Please input new batch ratio: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            if sum(c) == 1:
                                param_batch = c
                                break
                            else:
                                print('Format error, the sum of the proportions must be 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-3':
                    while 1:
                        a1 = input('Please input new batch variation: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            param_batch_var = c
                            break
                        except:
                            print('Format error, the input data should be a floating point number')


                elif a == '-4':
                    while 1:
                        try:
                            a1 = input('Please input path(Yes or No): ')
                            if a1 == 'Yes' or a1 == 'No':
                                param_path = a1
                                break
                            else:
                                print('Enter Yes or No')
                        except:
                            print('Format error')
                elif a == '-5':
                    while 1:
                        try:
                            a1 = input('Please input new DEgene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1:
                                param_DEgene = b
                                break
                            else:
                                print('DEgene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-6':
                    while 1:
                        a1 = input('Please input new DEgene variation: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            param_DEgene_var = c
                            break
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-7':
                    while 1:
                        try:
                            a1 = input('Please input new Non-linear gene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1 and b+param_DEgene <= 1:
                                param_DEgene_var = b
                                break
                            else:
                                print('Non-linear gene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-8':
                    while 1:
                        try:
                            a1 = input('Please input new marker gene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1 and b < param_DEgene:
                                param_marker = b
                                break
                            else:
                                print('Marker gene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-9':
                    while 1:
                        try:
                            a1 = input('Please input new library magnification: ')
                            b = float(a1)
                            if b > 0 :
                                param_library = b
                                break
                            else:
                                print('Library magnification should be greater than 0')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-10':
                    while 1:
                        try:
                            a1 = input('Please input new cell number: ')
                            b = int(a1)
                            if b > 0 :
                                cell_number = b
                                break
                            else:
                                print('Cell number should be greater than 0')
                        except:
                            print('Format error, the input data should be a positive integer')
                elif a == '-v':
                    print('-1\tGroup ratio (for extended simulation):'+str(param_group)+int(28-len(str(param_group)))*' ')
                    print('-2\tBatch ratio (for extended simulation):'+str(param_batch)+int(28-len(str(param_batch)))*' ')
                    print('-3\tBatch variation (for extended simulation):'+str(param_batch_var)+int(24-len(str(param_batch_var)))*' ')
                    print('-4\tPath (for extended simulation):'+str(param_path)+int(35-len(str(param_path)))*' ')
                    print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene)+int(27-len(str(param_DEgene)))*' ')
                    print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var)+int(23-len(str(param_DEgene_var)))*' ')
                    print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline)+int(18-len(str(param_noline)))*' ')
                    print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker)+int(22-len(str(param_marker)))*' ')
                    print('-9\tLibrary magnification (for extended simulation):'+str(param_library)+int(18-len(str(param_library)))*' ')
                    print('-10\tCell number (for extended simulation):'+str(cell_number)+int(28-len(str(cell_number)))*' ')
                elif a == '-s':
                    time_start=time.time()
                    creatcount_NBzero_noparam()
                    time_end=time.time()
                    print('Total cost',time_end-time_start)
                    while 1:
                        path3 = input('Please enter the saving path: ')
                        if os.path.exists(path3) == True:
                            break
                        else:
                            print('Path does not exist')
                    write_noparam()
                    print('Done!')
                    print('\n')
                    print('-m means simulate on this data again, otherwise exit the simulation')
                    a1 = input('Please input: ')
                    if a1 != '-m':
                        break
                    else:
                        print('-1\tGroup ratio (for extended simulation):'+str(param_group)+int(28-len(str(param_group)))*' ')
                        print('-2\tBatch ratio (for extended simulation):'+str(param_batch)+int(28-len(str(param_batch)))*' ')
                        print('-3\tBatch variation (for extended simulation):'+str(param_batch_var)+int(24-len(str(param_batch_var)))*' ')
                        print('-4\tPath (for extended simulation):'+str(param_path)+int(35-len(str(param_path)))*' ')
                        print('-5\tDEgene ratio (for extended simulation):'+str(param_DEgene)+int(27-len(str(param_DEgene)))*' ')
                        print('-6\tDEgene variation (for extended simulation):'+str(param_DEgene_var)+int(23-len(str(param_DEgene_var)))*' ')
                        print('-7\tNon-linear gene ratio (for extended simulation):'+str(param_noline)+int(18-len(str(param_noline)))*' ')
                        print('-8\tMarker gene ratio (for extended simulation):'+str(param_marker)+int(22-len(str(param_marker)))*' ')
                        print('-9\tLibrary magnification (for extended simulation):'+str(param_library)+int(18-len(str(param_library)))*' ')
                        print('-10\tCell number (for extended simulation):'+str(cell_number)+int(28-len(str(cell_number)))*' ')
                        print('\n')
                        print('If you want to modify the above parameters, please enter a number')
                        print('-v means view current parameters')
                        print('-s means start simulation')
                else:
                    print('Format error')




        elif mod == '-e':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')
            while 1:
                path2 = input('Please give the path to the cell group file: ')
                if os.path.exists(path2) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package: ')
                    seed2 = input('Please enter new seed2 for running the numpy package: ')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')




            copula_number = 2000
            Subgroup = 'No'

            print('-1\tCopula genes number:'+str(copula_number)+int(20-len(str(copula_number)))*' ')
            print('-2\tSubgroup:'+str(Subgroup))
            print('If you want to modify the above parameters, please enter a number')
            print('Enter other to start the simulation')
            while 1:
                a = input('Please input: ')
                if a == '-1':
                    while 1:
                        try:
                            a1 = input('Please input new copula genes number: ')
                            b = int(a1)
                            if b > 0 :
                                copula_number = b
                                break
                            else:
                                print('Copula genes number should be greater than 0')
                        except:
                            print('Format error, the input data should be a positive integer')
                if a == '-2':
                    while 1:
                        try:
                            a1 = input('If you want to input subgroup (Yes or No)? ')
                            if a1 == 'Yes' or a1 == 'No':
                                Subgroup = a1
                                break
                            else:
                                print('Input should be Yes or No')
                        except:
                            print('Format error')
                if a != '-1' and a != '-2':
                    break

            print('\n'+'Estimating parameters for SimCH-copula-NB simulation of independent multiple groups...')
            Remove0_t()
            nml_t()
            mean_t()
            BCV_t()
            ls1 = [len(i) for i in ls7_a]
            cell_number = sum(ls1)
            param_library = 1

            print('\n')
            print('If you want to modify the above parameters, please enter a number')
            print('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
            print('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')

            print('-v means view current parameters')
            print('-s means start simulation')
            while 1:
                a = input('Please input: ')
                if a == '-1':
                    while 1:
                        try:
                            a1 = input('Please input new cell number: ')
                            b = int(a1)
                            if b > 0 :
                                cell_number = b
                                break
                            else:
                                print('Cell number should be greater than 0')
                        except:
                            print('Format error, the input data should be a positive integer')
                elif a == '-2':
                    while 1:
                        try:
                            a1 = input('Please input new library magnification value: ')
                            b = float(a1)
                            if b > 0 :
                                param_library = b
                                break
                            else:
                                print('Library magnification should be greater than 0')
                        except:
                            print('Format error, the input data should be a floating point number')

                elif a == '-v':
                    print('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
                    print('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')
                    print('-v means view current parameters')
                    print('-s means start simulation')
                elif a == '-s':
                    while 1:
                        path3 = input('Please enter the saving path: ')
                        if os.path.exists(path3) == True:
                            break
                        else:
                            print('Path does not exist')
                    creatcount_noparam_t()
                    write_t()
                    print('Done！')
                else:
                    print('Format error')

        elif mod == '-e nocopula':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')

            while 1:
                path2 = input('Please give the path to the cell group file: ')
                if os.path.exists(path2) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package: ')
                    seed2 = input('Please enter new seed2 for running the numpy package: ')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')


            copula_number = 2000
            Subgroup = 'No'

            print('-1\tCopula genes number:'+str(copula_number))
            print('-2\tSubgroup:'+str(Subgroup))
            print('If you want to modify the above parameters, please enter a number')
            print('Enter other to start the simulation')
            a = input('Please input: ')
            if a == '-1':
                while 1:
                    try:
                        a1 = input('Please input new copula genes number: ')
                        b = int(a1)
                        if b > 0 :
                            copula_number = b
                            break
                        else:
                            print('Copula genes number should be greater than 0')
                    except:
                        print('Format error, the input data should be a positive integer')
            if a == '-2':
                while 1:
                    try:
                        a1 = input('If you want to input subgroup (Yes or No)? ')
                        if a1 == 'Yes' or a1 == 'No':
                            Subgroup = a1
                            break
                        else:
                            print('Input should be Yes or No')
                    except:
                        print('Format error')



            print('\n'+'Estimating parameters for SimCH-fit-NB simulation of independent multiple groups...')
            Remove0_t()
            nml_t()
            mean_t()
            BCV_t()
            ls1 = [len(i) for i in ls7_a]
            cell_number = sum(ls1)
            param_library = 1
            print('\n')
            print('If you want to modify the above parameters, please enter a number')
            print('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
            print('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')


            print('-v means view current parameters')
            print('-s means start simulation')
            while 1:
                a = input('Please input: ')
                if a == '-1':
                    while 1:
                        try:
                            a1 = input('Please input new cell number: ')
                            b = int(a1)
                            if b > 0 :
                                cell_number = b
                                break
                            else:
                                print('Cell number should be greater than 0')
                        except:
                            print('Format error, the input data should be a positive integer')
                elif a == '-2':
                    while 1:
                        try:
                            a1 = input('Please input new library magnification value: ')
                            b = float(a1)
                            if b > 0 :
                                param_library = b
                                break
                            else:
                                print('Library magnification value should be greater than 0')
                        except:
                            print('Format error, the input data should be a floating point number')

                elif a == '-v':
                    print('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
                    print('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')
                    print('-v means view current parameters')
                    print('-s means start simulation')
                elif a == '-s':
                    while 1:
                        path3 = input('Please enter the saving path: ')
                        if os.path.exists(path3) == True:
                            break
                        else:
                            print('Path does not exist')
                    creatcount_noparam_t()
                    write_t()
                    print('Done！')
                else:
                    print('Format error')
        elif mod == '-f':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')

            while 1:
                path1 = input('Please give the path to the cell group file: ')
                if os.path.exists(path2) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package: ')
                    seed2 = input('Please enter new seed2 for running the numpy package: ')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')

            copula_number = 2000
            Subgroup = 'No'

            print('-1\tCopula genes number:'+str(copula_number))
            print('-2\tSubgroup:'+str(Subgroup))
            print('If you want to modify the above parameters, please enter a number')
            print('Enter other to start the simulation')
            a = input('Please input: ')
            if a == '-1':
                while 1:
                    try:
                        a1 = input('Please input new copula genes number:')
                        b = int(a1)
                        if b > 0 :
                            copula_number = b
                            break
                        else:
                            print('Copula genes number should be greater than 0')
                    except:
                        print('Format error, the input data should be a positive integer')
            if a == '-2':
                while 1:
                    try:
                        a1 = input('If you want to input subgroup (Yes or No): ')
                        if a1 == 'Yes' or a1 == 'No':
                            Subgroup = a1
                            break
                        else:
                            print('Input should be Yes or No')
                    except:
                        print('Format error')


            print('\n'+'Estimating parameters for SimCH-copula-NBZI simulation of independent multiple groups...')
            Remove0_t()
            nml_t()
            mean_t()
            NBzero_mean_dispersion_noparam_t()
            ls1 = [len(i) for i in ls7_a]
            cell_number = sum(ls1)
            param_library = 1
            print('\n')
            print('If you want to modify the above parameters, please enter a number')
            print('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
            print('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')



            print('-v means view current parameters')
            print('-s means start simulation')
            while 1:
                a = input('Please input: ')
                if a == '-1':
                    while 1:
                        try:
                            a1 = input('Please input new cell number: ')
                            b = int(a1)
                            if b > 0 :
                                cell_number = b
                                break
                            else:
                                print('Cell number should be greater than 0')
                        except:
                            print('Format error, the input data should be a positive integer')
                elif a == '-2':
                    while 1:
                        try:
                            a1 = input('Please input new library magnification value: ')
                            b = float(a1)
                            if b > 0 :
                                param_library = b
                                break
                            else:
                                print('Library magnification value should be greater than 0')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-v':
                    print('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
                    print('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')
                    print('Otherwise continue to run')
                    print('-v means view current parameters')
                    print('-s means start simulation')
                elif a == '-s':
                    creatcount_NBzero_noparam_t()
                    while 1:
                        path3 = input('Please enter the saving path: ')
                        if os.path.exists(path3) == True:
                            break
                        else:
                            print('Path does not exist')
                    write_t()
                    print('Done！')
                else:
                    print('Format error')

        elif mod == '-f nocopula':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')

            while 1:
                path2 = input('Please give the path to the cell group file: ')
                if os.path.exists(path2) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package: ')
                    seed2 = input('Please enter new seed2 for running the numpy package: ')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')

            copula_number = 2000
            Subgroup = 'No'

            print('-1\tCopula genes number:'+str(copula_number))
            print('-2\tSubgroup:'+str(Subgroup))
            print('If you want to modify the above parameters, please enter a number')
            print('Enter other to start the simulation')
            a = input('Please input: ')
            if a == '-1':
                while 1:
                    try:
                        a1 = input('Please input new number of copula genes: ')
                        b = int(a1)
                        if b > 0 :
                            copula_number = b
                            break
                        else:
                            print('Copula genes number should be greater than 0')
                    except:
                        print('Format error, the input data should be a positive integer')
            if a == '-2':
                while 1:
                    try:
                        a1 = input('If you want to input subgroup? (Yes or No):')
                        if a1 == 'Yes' or a1 == 'No':
                            Subgroup = a1
                            break
                        else:
                            print('Input should be Yes or No')
                    except:
                        print('Format error')



            print('\n'+'Estimating parameters for SimCH-fit-NBZI simulation of independent multiple groups...')
            Remove0_t()
            nml_t()
            mean_t()
            NBzero_mean_dispersion_noparam_t()
            ls1 = [len(i) for i in ls7_a]
            cell_number = sum(ls1)
            param_library = 1
            print('\n')
            print('If you want to modify the above parameters, please enter a number (e.g. -2)')
            print('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
            print('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')
            print('Otherwise continue to run')

            print('-v means view current parameters')
            print('-s means start simulation')
            while 1:
                a = input('Please input: ')
                if a == '-1':
                    while 1:
                        try:
                            a1 = input('Please input new cell number: ')
                            b = int(a1)
                            if b > 0 :
                                cell_number = b
                                break
                            else:
                                print('Cell number should be greater than 0')
                        except:
                            print('Format error, the input data should be a positive integer')
                elif a == '-2':
                    while 1:
                        try:
                            a1 = input('Please input new library magnification: ')
                            b = float(a1)
                            if b > 0 :
                                param_library = b
                                break
                            else:
                                print('Library magnification should be greater than 0')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-v':
                    print('Cell number:'+str(cell_number)+int(28-len(str(cell_number)))*' '+'-1')
                    print('Library magnification:'+str(param_library)+int(18-len(str(param_library)))*' '+'-2')
                    print('-v means view current parameters')
                    print('-s means start simulation')
                elif a == '-s':
                    creatcount_NBzero_noparam_t()
                    while 1:
                        path3 = input('Please enter the saving path: ')
                        if os.path.exists(path3) == True:
                            break
                        else:
                            print('Path does not exist')
                    write_t()
                    print('Done！')
                else:
                    print('Format error')



        elif mod == '-g':
            while 1:
                path1 = input('Please give the path to the count matrix file: ')
                if os.path.exists(path1) == True:
                    break
                else:
                    print('File path does not exist')
            print('Change seeds for generating new series of random numbers.')
            mod4 = input('Please input "seed" to make change or just skip this: ')
            if mod4 == 'seed':
                while 1:
                    seed1 = input('Please enter new seed1 for running the random package: ')
                    seed2 = input('Please enter new seed2 for running the numpy package: ')
                    try:
                        seed_1 = int(seed1)
                        seed_2 = int(seed2)
                        random.seed(seed_1)
                        np.random.seed(seed_2)
                        break
                    except:
                        print('Input format error')
            copula_number = 2000
            print('-1\tCopula genes number:'+str(copula_number)+int(20-len(str(copula_number)))*' ')
            print('If you want to modify the above parameters, please enter a number')
            print('Enter other to start the simulation')
            a = input('Please input: ')
            if a == '-1':
                while 1:
                    try:
                        a1 = input('Please input new copula genes number: ')
                        b = int(a1)
                        if b > 0 :
                            copula_number = b
                            break
                        else:
                            print('Copula genes number should be greater than 0')
                    except:
                        print('Format error, the input data should be a positive integer')
            print('\n'+'Estimating parameters for SimCH-copula-NB simulation for testing imputation methods')
            Remove0_noparam()
            nml_noparam()
            mean_noparam()
            BCV_noparam()

            param_group = [1]
            param_batch = [1]
            param_batch_var = [0,0.2]
            param_path = 'No'
            param_DEgene = 0.2
            param_DEgene_var = [0,0.5]
            param_noline = 0.01
            param_marker = 0.01
            param_library = 1
            cell_number = len(ls7)

            print('\n')
            print('-1\tGroup ratio:'+str(param_group))
            print('-2\tBatch ratio:'+str(param_batch))
            print('-3\tBatch variation:'+str(param_batch_var))
            print('-4\tPath:'+str(param_path))
            print('-5\tDEgene ratio:'+str(param_DEgene))
            print('-6\tDEgene variation:'+str(param_DEgene_var))
            print('-7\tNon-linear gene ratio:'+str(param_noline))
            print('-8\tMarker gene ratio:'+str(param_marker))
            print('-9\tLibrary magnification:'+str(param_library))
            print('-10\tCell number:'+str(cell_number))

            print('\n')
            print('If you want to modify the above parameters, please enter a number')
            print('-v means view current parameters')
            print('-s means start simulation')

            while 1:
                a = input('Please input: ')
                if a == '-1':
                    while 1:
                        a1 = input('Please input new group ratio:')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            if sum(c) == 1:
                                param_group = c
                                break
                            else:
                                print('Format error, the sum of the proportions must be 1')
                        except:
                            c = []
                            d = []
                            for i in range(len(b)):
                                if '[' in b[i]:
                                    c.append(i)
                                if  ']' in b[i]:
                                    d.append(i)
                            f = []
                            g = []
                            j = 1
                            if (len(c) + len(d))%2 == 0:
                                e = a1[1:-1].split(',')
                                for i in range(len(e)):
                                    if '[' in e[i] or ']' in e[i] or j == -1:
                                        if '[' in e[i]:
                                            j *= -1
                                            try:
                                                g.append(float(e[i].replace('[','').replace(']','')))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                        elif ']' in e[i]:
                                            try:
                                                g.append(float(e[i].replace(']','').replace('[','')))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                            j *= -1
                                            f.append(g)
                                            g = []
                                        else:
                                            try:
                                                g.append(float(e[i]))
                                            except:
                                                print('Format error, the input data should be a floating point number')
                                    else:
                                        try:
                                            f.append(float(e[i]))
                                        except:
                                            print('Format error, the input data should be a floating point number')
                                i1 = 0
                                for i in range(len(f)):
                                    if type(f[i]) == type(0.1):
                                        i1 += f[i]
                                    elif type(f[i]) == type([]):
                                        i1 += sum(f[i])
                                if i1 == 1:
                                    param_group = f
                                    break
                                else:
                                    print('Format error, the sum of the proportions must be 1')
                            else:
                                print('Format error, the input data should be a floating point number')


                elif a == '-2':
                    while 1:
                        a1 = input('Please input new batch ratio: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            if sum(c) == 1:
                                param_batch = c
                                break
                            else:
                                print('Format error, the sum of the proportions must be 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-3':
                    while 1:
                        a1 = input('Please input new batch variation: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            param_batch_var = c
                            break
                        except:
                            print('Format error, the input data should be a floating point number')


                elif a == '-4':
                    while 1:
                        try:
                            a1 = input('Please input path(Yes or No): ')
                            if a1 == 'Yes' or a1 == 'No':
                                param_path = a1
                                break
                            else:
                                print('Enter Yes or No')
                        except:
                            print('Format error')
                elif a == '-5':
                    while 1:
                        try:
                            a1 = input('Please input new DEgene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1:
                                param_DEgene = b
                                break
                            else:
                                print('DEgene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-6':
                    while 1:
                        a1 = input('Please input new DEgene variation: ')
                        b = a1[1:-1].split(',')
                        try:
                            c = [float(i) for i in b]
                            param_DEgene_var = c
                            break
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-7':
                    while 1:
                        try:
                            a1 = input('Please input new Non-linear gene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1 and b+param_DEgene <= 1:
                                param_noline = b
                                break
                            else:
                                print('Non-linear gene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-8':
                    while 1:
                        try:
                            a1 = input('Please input new marker gene ratio: ')
                            b = float(a1)
                            if b >= 0 and b <= 1 and b < param_DEgene:
                                param_marker = b
                                break
                            else:
                                print('Marker gene ratio should be greater than 0 and less than 1')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-9':
                    while 1:
                        try:
                            a1 = input('Please input new library magnification: ')
                            b = float(a1)
                            if b > 0 :
                                param_library = b
                                break
                            else:
                                print('Library magnification should be greater than 0')
                        except:
                            print('Format error, the input data should be a floating point number')
                elif a == '-10':
                    while 1:
                        try:
                            a1 = input('Please input new cell number: ')
                            b = int(a1)
                            if b > 0 :
                                cell_number = b
                                break
                            else:
                                print('Cell number should be greater than 0')
                        except:
                            print('Format error, the input data should be a positive integer')
                elif a == '-v':
                    print('\n')
                    print('-1\tGroup ratio:'+str(param_group))
                    print('-2\tBatch ratio:'+str(param_batch))
                    print('-3\tBatch variation:'+str(param_batch_var))
                    print('-4\tPath:'+str(param_path))
                    print('-5\tDEgene ratio:'+str(param_DEgene))
                    print('-6\tDEgene variation:'+str(param_DEgene_var))
                    print('-7\tNon-linear gene ratio:'+str(param_noline))
                    print('-8\tMarker gene ratio:'+str(param_marker))
                    print('-9\tLibrary magnification:'+str(param_library))
                    print('-10\tCell number:'+str(cell_number))
                    print('\n')
                elif a == '-s':
                    creatcount_noparam_test_imputation_cor()
                    while 1:
                        path3 = input('Please enter the saving path: ')
                        if os.path.exists(path3) == True:
                            break
                        else:
                            print('Path does not exist')
                    write_noparam()
                    print('Done!')
                    print('\n')
                    print('-m means simulate on this data again, otherwise exit the simulation')
                    a1 = input('Please input: ')
                    if a1 != '-m':
                        break
                    else:
                        print('\n')
                        print('-1\tGroup ratio:'+str(param_group))
                        print('-2\tBatch ratio:'+str(param_batch))
                        print('-3\tBatch variation:'+str(param_batch_var))
                        print('-4\tPath:'+str(param_path))
                        print('-5\tDEgene ratio:'+str(param_DEgene))
                        print('-6\tDEgene variation:'+str(param_DEgene_var))
                        print('-7\tNon-linear gene ratio:'+str(param_noline))
                        print('-8\tMarker gene ratio:'+str(param_marker))
                        print('-9\tLibrary magnification:'+str(param_library))
                        print('-10\tCell number:'+str(cell_number))
                        print('If you want to modify the above parameters, please enter a number')
                        print('-v means view current parameters')
                        print('-s means start simulation')
                        print('\n')
                else:
                    print('Format error')
        elif mod == ' ':
            break
        else:
            print('Format error')
            print('\n')


