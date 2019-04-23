import itertools
import math
#from pylab import *

def getAccuracy_PHFAWMM(PHFAWMMandPHFAMM):
    sum = 0.0;
    U = PHFAWMMandPHFAMM[0][0]
    V = PHFAWMMandPHFAMM[0][1]
    for i in range(0,len(U)):
        sum += pow(U[i], 2) - pow(V[i], 2)
    result=sum / len(U)
    return result

def getScore_PHFAWMM(Accuracy):
    result=(1+Accuracy)/2.0
    return result

def factorial(n):   #表示求n的阶乘的函数
    result = n
    for i in range(1,n):
        result *= i
    return result

def WPHFAMM_Hamacher_5(HPs,w,Q,Theta):#Mλp,#HPs表示多个PHFNs，是一个多维数组，每一个元素代表一个HP
    UV = []
    hpu = []
    hpv = []
    #HPs是一个数组，有多个元素，每一个元素代表一个PH,假设HPs中有不超过5个HP
    U = []
    V = []
    eta = Theta-1
    qq = 0
    for qj in Q:
        qq += qj
    lamda = float(1 / qq)
    n = len(HPs)
    fn1 = 1 / factorial(n)
    for i in range(len(HPs)): #将UV分离出来，分别进行处理
        ui=[]
        vi=[]
        for j in range(len(HPs[i])):
            ui.append(HPs[i][j][0])
            vi.append(HPs[i][j][1])
        U.append(ui)
        V.append(vi)
    for hp0ui in U[0]: #
        for hp1ui in U[1]:  #
            if (len(HPs)>2):
                for hp2ui in U[2]:  #
                    if(len(HPs)>3):
                        for hp3ui in U[3]:  #
                            if (len(HPs) > 4):
                                for hp4ui in U[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui), str(hp4ui)]
                                    n = len(HPs)
                                    S1 = 1
                                    S2 = 1
                                    # 排列数总数为n!个
                                    for i in itertools.permutations(Ott, n):
                                        # 每一个i代表一个排列数，每一个排列数进行U计算
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            m10 = pow(Theta-eta*(1-pow(j[k],2)),n*w[k])
                                            m11 = pow(1-pow(j[k],2),n*w[k])
                                            M1 *= pow(m10+(pow(Theta,2)-1)*m11,Q[k])
                                            M2 *= pow(m10-m11,Q[k])
                                        S1 *= pow(M1+(pow(Theta,2)-1)*M2,fn1)
                                        S2 *= pow(M1-M2,fn1)
                                    ro1 = pow(S1-S2,lamda)
                                    ro2 = pow(S1+(pow(Theta,2)-1)*S2,lamda)
                                    u = math.sqrt((Theta*ro1)/(ro2+eta*ro1))  # log(b,a):表示以a底b的对数
                                    hpu.append(u)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui)]
                                n = len(HPs)
                                S1 = 1
                                S2 = 1
                                # 排列数总数为n!个
                                for i in itertools.permutations(Ott, n):
                                    # 每一个i代表一个排列数，每一个排列数进行U计算
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        m10 = pow(Theta - eta * (1 - pow(j[k], 2)), n * w[k])
                                        m11 = pow(1 - pow(j[k], 2), n * w[k])
                                        M1 *= pow(m10 + (pow(Theta, 2) - 1) * m11, Q[k])
                                        M2 *= pow(m10 - m11, Q[k])
                                    S1 *= pow(M1 + (pow(Theta, 2) - 1) * M2, fn1)
                                    S2 *= pow(M1 - M2, fn1)
                                ro1 = pow(S1 - S2, lamda)
                                ro2 = pow(S1 + (pow(Theta, 2) - 1) * S2, lamda)
                                u = math.sqrt((Theta * ro1) / (ro2 + eta * ro1))  # log(b,a):表示以a底b的对数
                                hpu.append(u)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0ui), str(hp1ui), str(hp2ui)]
                        n = len(HPs)
                        S1 = 1
                        S2 = 1
                        # 排列数总数为n!个
                        for i in itertools.permutations(Ott, n):
                            # 每一个i代表一个排列数，每一个排列数进行U计算
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                m10 = pow(Theta - eta * (1 - pow(j[k], 2)), n * w[k])
                                m11 = pow(1 - pow(j[k], 2), n * w[k])
                                M1 *= pow(m10 + (pow(Theta, 2) - 1) * m11, Q[k])
                                M2 *= pow(m10 - m11, Q[k])
                            S1 *= pow(M1 + (pow(Theta, 2) - 1) * M2, fn1)
                            S2 *= pow(M1 - M2, fn1)
                        ro1 = pow(S1 - S2, lamda)
                        ro2 = pow(S1 + (pow(Theta, 2) - 1) * S2, lamda)
                        u = math.sqrt((Theta * ro1) / (ro2 + eta * ro1))  # log(b,a):表示以a底b的对数
                        hpu.append(u)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0ui), str(hp1ui)]
                n = len(HPs)
                S1 = 1
                S2 = 1
                # 排列数总数为n!个
                for i in itertools.permutations(Ott, n):
                    # 每一个i代表一个排列数，每一个排列数进行U计算
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        m10 = pow(Theta - eta * (1 - pow(j[k], 2)), n * w[k])
                        m11 = pow(1 - pow(j[k], 2), n * w[k])
                        M1 *= pow(m10 + (pow(Theta, 2) - 1) * m11, Q[k])
                        M2 *= pow(m10 - m11, Q[k])
                    S1 *= pow(M1 + (pow(Theta, 2) - 1) * M2, fn1)
                    S2 *= pow(M1 - M2, fn1)
                ro1 = pow(S1 - S2, lamda)
                ro2 = pow(S1 + (pow(Theta, 2) - 1) * S2, lamda)
                u = math.sqrt((Theta * ro1) / (ro2 + eta * ro1))  # log(b,a):表示以a底b的对数
                hpu.append(u)
    # #计算V
    for hp0vi in V[0]:  #
        for hp1vi in V[1]:  #
            if (len(HPs) > 2):
                for hp2vi in V[2]:  #
                    if(len(HPs) > 3):
                        for hp3vi in V[3]:  #
                            if (len(HPs) > 4):
                                for hp4vi in V[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi), str(hp4vi)]
                                    S1 = 1
                                    S2 = 1
                                    for i in itertools.permutations(Ott, len(HPs)):
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            m10 = pow(Theta - eta * pow(j[k], 2), n * w[k])
                                            m11 = pow(j[k], 2 * n * w[k])
                                            M1 *= pow(m10 + (pow(Theta, 2) - 1) * m11, Q[k])
                                            M2 *= pow(m10 - m11, Q[k])
                                        S1 *= pow(M1 + (pow(Theta, 2) - 1) * M2, fn1)
                                        S2 *= pow(M1 - M2, fn1)
                                    ro1 = pow(S1 - S2, lamda)
                                    ro2 = pow(S1 + (pow(Theta, 2) - 1) * S2, lamda)
                                    v = math.sqrt((ro2 - ro1) / (ro2 + eta * ro1))
                                    hpv.append(v)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi)]
                                S1 = 1
                                S2 = 1
                                for i in itertools.permutations(Ott, len(HPs)):
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        m10 = pow(Theta - eta * pow(j[k], 2), n * w[k])
                                        m11 = pow( j[k], 2*n * w[k])
                                        M1 *= pow(m10 + (pow(Theta, 2) - 1) * m11, Q[k])
                                        M2 *= pow(m10 - m11, Q[k])
                                    S1 *= pow(M1 + (pow(Theta, 2) - 1) * M2, fn1)
                                    S2 *= pow(M1 - M2, fn1)
                                ro1 = pow(S1 - S2, lamda)
                                ro2 = pow(S1 + (pow(Theta, 2) - 1) * S2, lamda)
                                v = math.sqrt((ro2 - ro1) / (ro2 + eta*ro1))
                                hpv.append(v)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0vi), str(hp1vi), str(hp2vi)]
                        S1 = 1
                        S2 = 1
                        for i in itertools.permutations(Ott, len(HPs)):
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                m10 = pow(Theta - eta * pow(j[k], 2), n * w[k])
                                m11 = pow(j[k], 2 * n * w[k])
                                M1 *= pow(m10 + (pow(Theta, 2) - 1) * m11, Q[k])
                                M2 *= pow(m10 - m11, Q[k])
                            S1 *= pow(M1 + (pow(Theta, 2) - 1) * M2, fn1)
                            S2 *= pow(M1 - M2, fn1)
                        ro1 = pow(S1 - S2, lamda)
                        ro2 = pow(S1 + (pow(Theta, 2) - 1) * S2, lamda)
                        v = math.sqrt((ro2 - ro1) / (ro2 + eta * ro1))
                        hpv.append(v)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0vi), str(hp1vi)]
                S1 = 1
                S2 = 1
                for i in itertools.permutations(Ott, len(HPs)):
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        m10 = pow(Theta - eta * pow(j[k], 2), n * w[k])
                        m11 = pow(j[k], 2 * n * w[k])
                        M1 *= pow(m10 + (pow(Theta, 2) - 1) * m11, Q[k])
                        M2 *= pow(m10 - m11, Q[k])
                    S1 *= pow(M1 + (pow(Theta, 2) - 1) * M2, fn1)
                    S2 *= pow(M1 - M2, fn1)
                ro1 = pow(S1 - S2, lamda)
                ro2 = pow(S1 + (pow(Theta, 2) - 1) * S2, lamda)
                v = math.sqrt((ro2 - ro1) / (ro2 + eta * ro1))
                hpv.append(v)
    UV.append([hpu,hpv])
    return UV # 返回这个计算后的结果

def PHFAMM_Hamacher_5(HPs,w,Q,Theta):#Mλp,#HPs表示多个PHFNs，是一个多维数组，每一个元素代表一个HP
    UV = []
    hpu = []
    hpv = []
    #HPs是一个数组，有多个元素，每一个元素代表一个PH,假设HPs中有不超过5个HP
    U = []
    V = []
    eta = Theta-1
    qq = 0
    for qj in Q:
        qq += qj
    lamda = float(1 / qq)
    n = len(HPs)
    fn1 = 1 / factorial(n)
    for i in range(len(HPs)): #将UV分离出来，分别进行处理
        ui=[]
        vi=[]
        for j in range(len(HPs[i])):
            ui.append(HPs[i][j][0])
            vi.append(HPs[i][j][1])
        U.append(ui)
        V.append(vi)
    for hp0ui in U[0]: #
        for hp1ui in U[1]:  #
            if (len(HPs)>2):
                for hp2ui in U[2]:  #
                    if(len(HPs)>3):
                        for hp3ui in U[3]:  #
                            if (len(HPs) > 4):
                                for hp4ui in U[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui), str(hp4ui)]
                                    n = len(HPs)
                                    # 排列数总数为n!个
                                    S1 = 1
                                    S2 = 1
                                    for i in itertools.permutations(Ott, n):
                                        # 每一个i代表一个排列数，每一个排列数进行U计算
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            M1 *= pow(Theta-eta*pow(j[k],2),Q[k])
                                            M2 *= pow(j[k],2*Q[k])
                                        S1 *= pow(M1+(pow(eta,2)-1)*M2,fn1)
                                        S2 *= pow(M1-M2,fn1)
                                    ro1 = pow(S1-S2,lamda)
                                    ro2 = pow(S1+(pow(eta,2)-1)*S2,lamda)
                                    u = math.sqrt((Theta*ro1)/(ro2+eta*ro1))  # log(b,a):表示以a底b的对数
                                    hpu.append(u)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui)]
                                n = len(HPs)
                                # 排列数总数为n!个
                                S1 = 1
                                S2 = 1
                                for i in itertools.permutations(Ott, n):
                                    # 每一个i代表一个排列数，每一个排列数进行U计算
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        M1 *= pow(Theta - eta * pow(j[k], 2), Q[k])
                                        M2 *= pow(j[k], 2 * Q[k])
                                    S1 *= pow(M1 + (pow(eta, 2) - 1) * M2, fn1)
                                    S2 *= pow(M1 - M2, fn1)
                                ro1 = pow(S1 - S2, lamda)
                                ro2 = pow(S1 + (pow(eta, 2) - 1) * S2, lamda)
                                u = math.sqrt((Theta * ro1) / (ro2 + eta * ro1))  # log(b,a):表示以a底b的对数
                                hpu.append(u)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0ui), str(hp1ui), str(hp2ui)]
                        n = len(HPs)
                        # 排列数总数为n!个
                        S1 = 1
                        S2 = 1
                        for i in itertools.permutations(Ott, n):
                            # 每一个i代表一个排列数，每一个排列数进行U计算
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                M1 *= pow(Theta - eta * pow(j[k], 2), Q[k])
                                M2 *= pow(j[k], 2 * Q[k])
                            S1 *= pow(M1 + (pow(eta, 2) - 1) * M2, fn1)
                            S2 *= pow(M1 - M2, fn1)
                        ro1 = pow(S1 - S2, lamda)
                        ro2 = pow(S1 + (pow(eta, 2) - 1) * S2, lamda)
                        u = math.sqrt((Theta * ro1) / (ro2 + eta * ro1))  # log(b,a):表示以a底b的对数
                        hpu.append(u)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0ui), str(hp1ui)]
                n = len(HPs)
                # 排列数总数为n!个
                S1 = 1
                S2 = 1
                for i in itertools.permutations(Ott, n):
                    # 每一个i代表一个排列数，每一个排列数进行U计算
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        M1 *= pow(Theta - eta * pow(j[k], 2), Q[k])
                        M2 *= pow(j[k], 2 * Q[k])
                    S1 *= pow(M1 + (pow(eta, 2) - 1) * M2, fn1)
                    S2 *= pow(M1 - M2, fn1)
                ro1 = pow(S1 - S2, lamda)
                ro2 = pow(S1 + (pow(eta, 2) - 1) * S2, lamda)
                u = math.sqrt((Theta * ro1) / (ro2 + eta * ro1))  # log(b,a):表示以a底b的对数
                hpu.append(u)
    # #计算V
    for hp0vi in V[0]:  #
        for hp1vi in V[1]:  #
            if (len(HPs) > 2):
                for hp2vi in V[2]:  #
                    if(len(HPs) > 3):
                        for hp3vi in V[3]:  #
                            if (len(HPs) > 4):
                                for hp4vi in V[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi), str(hp4vi)]
                                    S1 = 1
                                    S2 = 1
                                    for i in itertools.permutations(Ott, len(HPs)):
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            M1 *= pow(Theta - eta * (1-pow(j[k], 2)), Q[k])
                                            M2 *= pow(1-pow(j[k],2),Q[k])
                                        S1 *= pow(M1-M2,fn1)
                                        S2 *= pow(M1+(pow(Theta,2)-1)*M2,fn1)
                                    ro1 = pow((pow(Theta,2)-1)*S1+S2,lamda)
                                    ro2 = pow(S2-S1,lamda)
                                    v = math.sqrt((ro1-ro2)/(ro1+eta*ro2))
                                    hpv.append(v)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi)]
                                S1 = 1
                                S2 = 1
                                for i in itertools.permutations(Ott, len(HPs)):
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        M1 *= pow(Theta - eta * (1 - pow(j[k], 2)), Q[k])
                                        M2 *= pow(1 - pow(j[k], 2), Q[k])
                                    S1 *= pow(M1 - M2, fn1)
                                    S2 *= pow(M1 + (pow(Theta, 2) - 1) * M2, fn1)
                                ro1 = pow((pow(Theta, 2) - 1) * S1 + S2, lamda)
                                ro2 = pow(S2 - S1, lamda)
                                v = math.sqrt((ro1 - ro2) / (ro1 + eta * ro2))
                                hpv.append(v)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0vi), str(hp1vi), str(hp2vi)]
                        S1 = 1
                        S2 = 1
                        for i in itertools.permutations(Ott, len(HPs)):
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                M1 *= pow(Theta - eta * (1 - pow(j[k], 2)), Q[k])
                                M2 *= pow(1 - pow(j[k], 2), Q[k])
                            S1 *= pow(M1 - M2, fn1)
                            S2 *= pow(M1 + (pow(Theta, 2) - 1) * M2, fn1)
                        ro1 = pow((pow(Theta, 2) - 1) * S1 + S2, lamda)
                        ro2 = pow(S2 - S1, lamda)
                        v = math.sqrt((ro1 - ro2) / (ro1 + eta * ro2))
                        hpv.append(v)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0vi), str(hp1vi)]
                S1 = 1
                S2 = 1
                for i in itertools.permutations(Ott, len(HPs)):
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        M1 *= pow(Theta - eta * (1 - pow(j[k], 2)), Q[k])
                        M2 *= pow(1 - pow(j[k], 2), Q[k])
                    S1 *= pow(M1 - M2, fn1)
                    S2 *= pow(M1 + (pow(Theta, 2) - 1) * M2, fn1)
                ro1 = pow((pow(Theta, 2) - 1) * S1 + S2, lamda)
                ro2 = pow(S2 - S1, lamda)
                v = math.sqrt((ro1 - ro2) / (ro1 + eta * ro2))
                hpv.append(v)
    UV.append([hpu,hpv])
    return UV # 返回这个计算后的结果

def WPHFAMM_Aglebraic_5(HPs,w,Q,Theta):#Mλp,#HPs表示多个PHFNs，是一个多维数组，每一个元素代表一个HP
    UV = []
    hpu = []
    hpv = []
    #HPs是一个数组，有多个元素，每一个元素代表一个PH,假设HPs中有不超过5个HP
    #计算U
    #print(HPs)
    U = []
    V = []
    for i in range(len(HPs)):
        ui=[]
        vi=[]
        for j in range(len(HPs[i])):
            ui.append(HPs[i][j][0])
            vi.append(HPs[i][j][1])
        U.append(ui)
        V.append(vi)
    for hp0ui in U[0]: #
        for hp1ui in U[1]:  #
            if (len(HPs)>2):
                for hp2ui in U[2]:  #
                    if(len(HPs)>3):
                        for hp3ui in U[3]:  #
                            if (len(HPs) > 4):
                                for hp4ui in U[4]:  #
                                    #在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui), str(hp4ui)]
                                    BB = 1;
                                    n = len(HPs)
                                    for i in itertools.permutations(Ott, n):
                                        #每一个i代表一个排列数
                                        #将i中元素进行A运算
                                        m0 = pow(pow(1 - pow(1 - pow(float(i[0]), 2 ),(n*w[0])),0.5),Q[0])
                                        m1 = pow(pow(1 - pow(1 - pow(float(i[1]), 2), (n * w[1])), 0.5), Q[1])
                                        m2 = pow(pow(1 - pow(1 - pow(float(i[2]), 2), (n * w[2])), 0.5), Q[2])
                                        m3 = pow(pow(1 - pow(1 - pow(float(i[3]), 2), (n * w[3])), 0.5), Q[3])
                                        m4 = pow(pow(1 - pow(1 - pow(float(i[4]), 2), (n * w[4])), 0.5), Q[4])
                                        a=1-pow(m0*m1*m2*m3*m4,2)
                                        #将AA中所有元素进行B运算
                                        BB *= a;
                                    qq = 0.0
                                    for qj in Q:
                                        qq +=qj

                                    u = pow(pow(1-pow(BB,float(1 / factorial(n))),float(1 / qq)),0.5)#factorial表示求n的阶乘的函数
                                    #u = pow(pow(1 - pow(BB,float(1 / factorial(n))),2),float(1 / qq))#factorial表示求n的阶乘的函数
                                    hpu.append(u)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui)]
                                BB = 1;
                                n = len(HPs)
                                for i in itertools.permutations(Ott, len(Ott)):
                                    # 每一个i代表一个排列数
                                    # 将i中元素进行A运算
                                    m0 = pow(pow(1 - pow(1 - pow(float(i[0]), 2), (n * w[0])), 0.5), Q[0])
                                    m1 = pow(pow(1 - pow(1 - pow(float(i[1]), 2), (n * w[1])), 0.5), Q[1])
                                    m2 = pow(pow(1 - pow(1 - pow(float(i[2]), 2), (n * w[2])), 0.5), Q[2])
                                    m3 = pow(pow(1 - pow(1 - pow(float(i[3]), 2), (n * w[3])), 0.5), Q[3])
                                    a = 1 - pow(m0 * m1 * m2 * m3, 2)
                                    # 将AA中所有元素进行B运算
                                    BB *= a;
                                qq = 0
                                for qj in Q:
                                    qq += qj
                                #u = pow(pow(1 - pow(BB, float(1 / factorial(n))), 2),float(1 / qq))  # factorial表示求n的阶乘的函数
                                u = pow(pow(1 - pow(BB, float(1 / factorial(n))), float(1 / qq)),
                                        0.5)  # factorial表示求n的阶乘的函数
                                hpu.append(u)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0ui), str(hp1ui), str(hp2ui)]
                        BB = 1;
                        n = len(HPs)
                        for i in itertools.permutations(Ott, n - 1):
                            # 每一个i代表一个排列数
                            # 将i中元素进行A运算
                            m0 = pow(pow(1 - pow(1 - pow(float(i[0], 2)), (n * w[0])), 0.5), Q[0])
                            m1 = pow(pow(1 - pow(1 - pow(float(i[1], 2)), (n * w[1])), 0.5), Q[1])
                            m2 = pow(pow(1 - pow(1 - pow(float(i[2], 2)), (n * w[2])), 0.5), Q[2])
                            a = 1 - pow(m0 * m1 * m2 , 2)
                            # 将AA中所有元素进行B运算
                            BB *= a;
                        qq = 0
                        for qj in Q:
                            qq += qj
                        #u = pow(pow(1 - pow(BB, float(1 / factorial(n))), 2), float(1 / qq))  # factorial表示求n的阶乘的函数
                        u = pow(pow(1 - pow(BB, float(1 / factorial(n))), float(1 / qq)), 0.5)  # factorial表示求n的阶乘的函数
                        hpu.append(u)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0ui), str(hp1ui)]
                BB = 1;
                n = len(HPs)
                for i in itertools.permutations(Ott, n - 1):
                    # 每一个i代表一个排列数
                    # 将i中元素进行A运算
                    m0 = pow(pow(1 - pow(1 - pow(float(i[0], 2)), (n * w[0])), 0.5), Q[0])
                    m1 = pow(pow(1 - pow(1 - pow(float(i[1], 2)), (n * w[1])), 0.5), Q[1])
                    a = 1 - pow(m0 * m1, 2)
                    # 将AA中所有元素进行B运算
                    BB *= a;
                qq = 0
                for qj in Q:
                    qq += qj
                #u = pow(pow(1 - pow(BB, float(1 / factorial(n))), 2), float(1 / qq))  # factorial表示求n的阶乘的函数
                u = pow(pow(1 - pow(BB, float(1 / factorial(n))), float(1 / qq)), 0.5)  # factorial表示求n的阶乘的函数
                hpu.append(u)
    # #计算V
    for hp0vi in V[0]:  #
        for hp1vi in V[1]:  #
            if (len(HPs) > 2):
                for hp2vi in V[2]:  #
                    if(len(HPs) > 3):
                        for hp3vi in V[3]:  #
                            if (len(HPs) > 4):
                                for hp4vi in V[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi), str(hp4vi)]
                                    BB = 1;
                                    n = len(HPs)
                                    for i in itertools.permutations(Ott, n - 1):
                                        # 每一个i代表一个排列数
                                        # 将i中元素进行A运算
                                        m0 = pow(1 - pow(float(i[0]), 2 * n * w[0]), Q[0])
                                        m1 = pow(1 - pow(float(i[1]), 2 * n * w[1]), Q[1])
                                        m2 = pow(1 - pow(float(i[2]), 2 * n * w[2]), Q[2])
                                        m3 = pow(1 - pow(float(i[3]), 2 * n * w[3]), Q[3])
                                        m4 = pow(1 - pow(float(i[4]), 2 * n * w[4]), Q[4])
                                        a=pow(1-m0*m1*m2*m3*m4,0.5)
                                        # 将元素进行B运算
                                        BB *=a
                                    qq = 0
                                    for qj in Q:
                                        qq += qj
                                    v=pow(1-pow(1-pow(pow(BB,float(1/factorial(n))),2),float(1/qq)),0.5)# factorial表示求n的阶乘的函数
                                    hpv.append(v)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi)]
                                BB = 1;
                                n = len(HPs)
                                for i in itertools.permutations(Ott, len(Ott)):
                                    # 每一个i代表一个排列数
                                    # 将i中元素进行A运算
                                    m0 = pow(1 - pow(float(i[0]), 2 * n * w[0]), Q[0])
                                    m1 = pow(1 - pow(float(i[1]), 2 * n * w[1]), Q[1])
                                    m2 = pow(1 - pow(float(i[2]), 2 * n * w[2]), Q[2])
                                    m3 = pow(1 - pow(float(i[3]), 2 * n * w[3]), Q[3])
                                    a = pow(1 - m0 * m1 * m2 * m3, 0.5)
                                    # 将元素进行B运算
                                    BB *= a
                                qq = 0
                                for qj in Q:
                                    qq += qj
                                v = pow(1 - pow(1 - pow(pow(BB, float(1 / factorial(n))), 2), float(1 / qq)),0.5)  # factorial表示求n的阶乘的函数
                                hpv.append(v)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0vi), str(hp1vi), str(hp2vi)]
                        BB = 1;
                        n = len(HPs)
                        for i in itertools.permutations(Ott, n - 1):
                            # 每一个i代表一个排列数
                            # 将i中元素进行A运算
                            m0 = pow(1 - pow(float(i[0]), 2 * n * w[0]), Q[0])
                            m1 = pow(1 - pow(float(i[1]), 2 * n * w[1]), Q[1])
                            m2 = pow(1 - pow(float(i[2]), 2 * n * w[2]), Q[2])
                            a = pow(1 - m0 * m1 * m2, 0.5)
                            # 将元素进行B运算
                            BB *= a
                        qq = 0
                        for qj in Q:
                            qq += qj
                        v = pow(1 - pow(1 - pow(pow(BB, float(1 / factorial(n))), 2), float(1 / qq)),0.5)  # factorial表示求n的阶乘的函数
                        hpv.append(v)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0vi), str(hp1vi)]
                BB = 1;
                n = len(HPs)
                for i in itertools.permutations(Ott, n - 1):
                    # 每一个i代表一个排列数
                    # 将i中元素进行A运算
                    m0 = pow(1 - pow(float(i[0]), 2 * n * w[0]), Q[0])
                    m1 = pow(1 - pow(float(i[1]), 2 * n * w[1]), Q[1])
                    a = pow(1 - m0 * m1, 0.5)
                    # 将元素进行B运算
                    BB *= a
                qq = 0
                for qj in Q:
                    qq += qj
                v = pow(1 - pow(1 - pow(pow(BB, float(1 / factorial(n))), 2), float(1 / qq)),0.5)  # factorial表示求n的阶乘的函数
                hpv.append(v)
    UV.append([hpu,hpv])
    return UV # 返回这个计算后的结果

def PHFAMM_Aglebraic_5(HPs,w,Q,Theta):#Mλp,#HPs表示多个PHFNs，是一个多维数组，每一个元素代表一个HP
    UV = []
    hpu = []
    hpv = []
    #HPs是一个数组，有多个元素，每一个元素代表一个PH,假设HPs中有不超过5个HP
    U = []
    V = []
    eta = Theta-1
    qq = 0
    for qj in Q:
        qq += qj
    lamda = float(1 / qq)
    n = len(HPs)
    fn1 = 1 / factorial(n)
    for i in range(len(HPs)): #将UV分离出来，分别进行处理
        ui=[]
        vi=[]
        for j in range(len(HPs[i])):
            ui.append(HPs[i][j][0])
            vi.append(HPs[i][j][1])
        U.append(ui)
        V.append(vi)
    for hp0ui in U[0]: #
        for hp1ui in U[1]:  #
            if (len(HPs)>2):
                for hp2ui in U[2]:  #
                    if(len(HPs)>3):
                        for hp3ui in U[3]:  #
                            if (len(HPs) > 4):
                                for hp4ui in U[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui), str(hp4ui)]
                                    n = len(HPs)
                                    UI = []
                                    # 排列数总数为n!个
                                    S = 1
                                    for i in itertools.permutations(Ott, n):
                                        # 每一个i代表一个排列数，每一个排列数进行U计算
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M = 1
                                        for k in range(0, len(j)):
                                            m = pow(j[k],2*Q[k])
                                            M *= m
                                        s = pow(1-M,fn1)
                                        S *= s
                                    u = math.sqrt(pow(1-S,lamda))  # log(b,a):表示以a底b的对数
                                    hpu.append(u)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui)]
                                n = len(HPs)
                                UI = []
                                # 排列数总数为n!个
                                S = 1
                                for i in itertools.permutations(Ott, n):
                                    # 每一个i代表一个排列数，每一个排列数进行U计算
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M = 1
                                    for k in range(0, len(j)):
                                        m = pow(j[k], 2 * Q[k])
                                        M *= m
                                    s = pow(1 - M, fn1)
                                    S *= s
                                u = math.sqrt(pow(1 - S, lamda))  # log(b,a):表示以a底b的对数
                                hpu.append(u)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0ui), str(hp1ui), str(hp2ui)]
                        n = len(HPs)
                        UI = []
                        # 排列数总数为n!个
                        S = 1
                        for i in itertools.permutations(Ott, n):
                            # 每一个i代表一个排列数，每一个排列数进行U计算
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M = 1
                            for k in range(0, len(j)):
                                m = pow(j[k], 2 * Q[k])
                                M *= m
                            s = pow(1 - M, fn1)
                            S *= s
                        u = math.sqrt(pow(1 - S, lamda))  # log(b,a):表示以a底b的对数
                        hpu.append(u)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0ui), str(hp1ui)]
                n = len(HPs)
                UI = []
                # 排列数总数为n!个
                S = 1
                for i in itertools.permutations(Ott, n):
                    # 每一个i代表一个排列数，每一个排列数进行U计算
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M = 1
                    for k in range(0, len(j)):
                        m = pow(j[k], 2 * Q[k])
                        M *= m
                    s = pow(1 - M, fn1)
                    S *= s
                u = math.sqrt(pow(1 - S, lamda))  # log(b,a):表示以a底b的对数
                hpu.append(u)
    # #计算V
    for hp0vi in V[0]:  #
        for hp1vi in V[1]:  #
            if (len(HPs) > 2):
                for hp2vi in V[2]:  #
                    if(len(HPs) > 3):
                        for hp3vi in V[3]:  #
                            if (len(HPs) > 4):
                                for hp4vi in V[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi), str(hp4vi)]
                                    n = len(HPs)
                                    VI = []
                                    S = 1
                                    for i in itertools.permutations(Ott, len(HPs)):
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M = 1
                                        for k in range(0, len(j)):
                                            m = pow(1-pow(j[k],2),Q[k])
                                            M *= m
                                        s = pow(1-M,fn1)
                                        S *=s
                                    ro = pow(1-S,lamda)
                                    v = math.sqrt(1-ro)
                                    hpv.append(v)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi)]
                                n = len(HPs)
                                VI = []
                                S = 1
                                for i in itertools.permutations(Ott, len(HPs)):
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M = 1
                                    for k in range(0, len(j)):
                                        m = pow(1 - pow(j[k], 2), Q[k])
                                        M *= m
                                    s = pow(1 - M, fn1)
                                    S *= s
                                ro = pow(1 - S, lamda)
                                v = math.sqrt(1 - ro)
                                hpv.append(v)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0vi), str(hp1vi), str(hp2vi)]
                        n = len(HPs)
                        VI = []
                        S = 1
                        for i in itertools.permutations(Ott, len(HPs)):
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M = 1
                            for k in range(0, len(j)):
                                m = pow(1 - pow(j[k], 2), Q[k])
                                M *= m
                            s = pow(1 - M, fn1)
                            S *= s
                        ro = pow(1 - S, lamda)
                        v = math.sqrt(1 - ro)
                        hpv.append(v)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0vi), str(hp1vi)]
                n = len(HPs)
                VI = []
                S = 1
                for i in itertools.permutations(Ott, len(HPs)):
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M = 1
                    for k in range(0, len(j)):
                        m = pow(1 - pow(j[k], 2), Q[k])
                        M *= m
                    s = pow(1 - M, fn1)
                    S *= s
                ro = pow(1 - S, lamda)
                v = math.sqrt(1 - ro)
                hpv.append(v)
    UV.append([hpu,hpv])
    return UV # 返回这个计算后的结果

def PHFAMM_Einstein_5(HPs,w,Q,Theta):#Mλp,#HPs表示多个PHFNs，是一个多维数组，每一个元素代表一个HP
    UV = []
    hpu = []
    hpv = []
    #HPs是一个数组，有多个元素，每一个元素代表一个PH,假设HPs中有不超过5个HP
    U = []
    V = []
    eta = Theta-1
    qq = 0
    for qj in Q:
        qq += qj
    lamda = float(1 / qq)
    n = len(HPs)
    fn1 = 1 / factorial(n)
    for i in range(len(HPs)): #将UV分离出来，分别进行处理
        ui=[]
        vi=[]
        for j in range(len(HPs[i])):
            ui.append(HPs[i][j][0])
            vi.append(HPs[i][j][1])
        U.append(ui)
        V.append(vi)
    for hp0ui in U[0]: #
        for hp1ui in U[1]:  #
            if (len(HPs)>2):
                for hp2ui in U[2]:  #
                    if(len(HPs)>3):
                        for hp3ui in U[3]:  #
                            if (len(HPs) > 4):
                                for hp4ui in U[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui), str(hp4ui)]
                                    n = len(HPs)
                                    Sa = 1
                                    Sb = 1
                                    # 排列数总数为n!个
                                    for i in itertools.permutations(Ott, n):
                                        # 每一个i代表一个排列数，每一个排列数进行U计算
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            M1 *= pow(2-pow(j[k],2),Q[k])
                                            M2 *= pow(j[k],2*Q[k])
                                        Sa *= pow(M1+3*M2,fn1)
                                        Sb *= pow(M1-M2,fn1)
                                    mol = 2*pow(Sa-Sb,lamda)
                                    del1 = pow(Sa+3*Sb,lamda)
                                    del2 = pow(Sa-Sb,lamda)
                                    u = math.sqrt(mol/(del1+del2))  # log(b,a):表示以a底b的对数
                                    hpu.append(u)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui)]
                                n = len(HPs)
                                Sa = 1
                                Sb = 1
                                # 排列数总数为n!个
                                for i in itertools.permutations(Ott, n):
                                    # 每一个i代表一个排列数，每一个排列数进行U计算
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        M1 *= pow(2 - pow(j[k], 2), Q[k])
                                        M2 *= pow(j[k], 2 * Q[k])
                                    Sa *= pow(M1 + 3 * M2, fn1)
                                    Sb *= pow(M1 - M2, fn1)
                                mol = 2 * pow(Sa - Sb, lamda)
                                del1 = pow(Sa + 3 * Sb, lamda)
                                del2 = pow(Sa - Sb, lamda)
                                u = math.sqrt(mol / (del1 + del2))  # log(b,a):表示以a底b的对数
                                hpu.append(u)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0ui), str(hp1ui), str(hp2ui)]
                        n = len(HPs)
                        Sa = 1
                        Sb = 1
                        # 排列数总数为n!个
                        for i in itertools.permutations(Ott, n):
                            # 每一个i代表一个排列数，每一个排列数进行U计算
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                M1 *= pow(2 - pow(j[k], 2), Q[k])
                                M2 *= pow(j[k], 2 * Q[k])
                            Sa *= pow(M1 + 3 * M2, fn1)
                            Sb *= pow(M1 - M2, fn1)
                        mol = 2 * pow(Sa - Sb, lamda)
                        del1 = pow(Sa + 3 * Sb, lamda)
                        del2 = pow(Sa - Sb, lamda)
                        u = math.sqrt(mol / (del1 + del2))  # log(b,a):表示以a底b的对数
                        hpu.append(u)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0ui), str(hp1ui)]
                n = len(HPs)
                Sa = 1
                Sb = 1
                # 排列数总数为n!个
                for i in itertools.permutations(Ott, n):
                    # 每一个i代表一个排列数，每一个排列数进行U计算
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        M1 *= pow(2 - pow(j[k], 2), Q[k])
                        M2 *= pow(j[k], 2 * Q[k])
                    Sa *= pow(M1 + 3 * M2, fn1)
                    Sb *= pow(M1 - M2, fn1)
                mol = 2 * pow(Sa - Sb, lamda)
                del1 = pow(Sa + 3 * Sb, lamda)
                del2 = pow(Sa - Sb, lamda)
                u = math.sqrt(mol / (del1 + del2))  # log(b,a):表示以a底b的对数
                hpu.append(u)
    # #计算V
    for hp0vi in V[0]:  #
        for hp1vi in V[1]:  #
            if (len(HPs) > 2):
                for hp2vi in V[2]:  #
                    if(len(HPs) > 3):
                        for hp3vi in V[3]:  #
                            if (len(HPs) > 4):
                                for hp4vi in V[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi), str(hp4vi)]
                                    Sa = 1
                                    Sb = 1
                                    for i in itertools.permutations(Ott, len(HPs)):
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            M1 *= pow(1+pow(j[k],2),Q[k])
                                            M2 *= pow(1-pow(j[k],2),Q[k])
                                        Sa *= pow(M1+3*M2,fn1)
                                        Sb *= pow(M1-M2,fn1)
                                    A = pow(Sa+3*Sb,lamda)
                                    B = pow(Sa-Sb,lamda)
                                    v = math.sqrt((A-B)/(A+B))
                                    hpv.append(v)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi)]
                                Sa = 1
                                Sb = 1
                                for i in itertools.permutations(Ott, len(HPs)):
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        M1 *= pow(1 + pow(j[k], 2), Q[k])
                                        M2 *= pow(1 - pow(j[k], 2), Q[k])
                                    Sa *= pow(M1 + 3 * M2, fn1)
                                    Sb *= pow(M1 - M2, fn1)
                                A = pow(Sa + 3 * Sb, lamda)
                                B = pow(Sa - Sb, lamda)
                                v = math.sqrt((A - B) / (A + B))
                                hpv.append(v)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0vi), str(hp1vi), str(hp2vi)]
                        Sa = 1
                        Sb = 1
                        for i in itertools.permutations(Ott, len(HPs)):
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                M1 *= pow(1 + pow(j[k], 2), Q[k])
                                M2 *= pow(1 - pow(j[k], 2), Q[k])
                            Sa *= pow(M1 + 3 * M2, fn1)
                            Sb *= pow(M1 - M2, fn1)
                        A = pow(Sa + 3 * Sb, lamda)
                        B = pow(Sa - Sb, lamda)
                        v = math.sqrt((A - B) / (A + B))
                        hpv.append(v)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0vi), str(hp1vi)]
                Sa = 1
                Sb = 1
                for i in itertools.permutations(Ott, len(HPs)):
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        M1 *= pow(1 + pow(j[k], 2), Q[k])
                        M2 *= pow(1 - pow(j[k], 2), Q[k])
                    Sa *= pow(M1 + 3 * M2, fn1)
                    Sb *= pow(M1 - M2, fn1)
                A = pow(Sa + 3 * Sb, lamda)
                B = pow(Sa - Sb, lamda)
                v = math.sqrt((A - B) / (A + B))
                hpv.append(v)
    UV.append([hpu,hpv])
    return UV # 返回这个计算后的结果

def WPHFAMM_Einstein_5(HPs,w,Q,Theta):#Mλp,#HPs表示多个PHFNs，是一个多维数组，每一个元素代表一个HP
    UV = []
    hpu = []
    hpv = []
    #HPs是一个数组，有多个元素，每一个元素代表一个PH,假设HPs中有不超过5个HP
    U = []
    V = []
    eta = Theta-1
    qq = 0
    for qj in Q:
        qq += qj
    lamda = float(1 / qq)
    n = len(HPs)
    fn1 = 1 / factorial(n)
    for i in range(len(HPs)): #将UV分离出来，分别进行处理
        ui=[]
        vi=[]
        for j in range(len(HPs[i])):
            ui.append(HPs[i][j][0])
            vi.append(HPs[i][j][1])
        U.append(ui)
        V.append(vi)
    for hp0ui in U[0]: #
        for hp1ui in U[1]:  #
            if (len(HPs)>2):
                for hp2ui in U[2]:  #
                    if(len(HPs)>3):
                        for hp3ui in U[3]:  #
                            if (len(HPs) > 4):
                                for hp4ui in U[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui), str(hp4ui)]
                                    n = len(HPs)
                                    # 排列数总数为n!个
                                    S1 = 1
                                    S2 = 1
                                    for i in itertools.permutations(Ott, n):
                                        # 每一个i代表一个排列数，每一个排列数进行U计算
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            m10 = pow(1+pow(j[k],2),n*w[k])
                                            m11 = pow(1-pow(j[k],2),n*w[k])
                                            m1 = pow(m10+3*m11,Q[k])
                                            m2 = pow(m10-m11,Q[k])
                                            M1 *= m1
                                            M2 *= m2
                                        s1 = pow(M1+3*M2,fn1)
                                        s2 = pow(M1-M2,fn1)
                                        S1 *= s1
                                        S2 *= s2
                                    mols = 2*pow(S1-S2,lamda)
                                    dels = pow(S1+3*S2,lamda)+pow(S1-S2,lamda)
                                    u = math.sqrt(mols/dels)
                                    hpu.append(u)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui)]
                                n = len(HPs)
                                # 排列数总数为n!个
                                S1 = 1
                                S2 = 1
                                for i in itertools.permutations(Ott, n):
                                    # 每一个i代表一个排列数，每一个排列数进行U计算
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        m10 = pow(1 + pow(j[k], 2), n * w[k])
                                        m11 = pow(1 - pow(j[k], 2), n * w[k])
                                        m1 = pow(m10 + 3 * m11, Q[k])
                                        m2 = pow(m10 - m11, Q[k])
                                        M1 *= m1
                                        M2 *= m2
                                    s1 = pow(M1 + 3 * M2, fn1)
                                    s2 = pow(M1 - M2, fn1)
                                    S1 *= s1
                                    S2 *= s2
                                mols = 2 * pow(S1 - S2, lamda)
                                dels = pow(S1 + 3 * S2, lamda) + pow(S1 - S2, lamda)
                                u = math.sqrt(mols / dels)
                                hpu.append(u)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0ui), str(hp1ui), str(hp2ui)]
                        n = len(HPs)
                        # 排列数总数为n!个
                        S1 = 1
                        S2 = 1
                        for i in itertools.permutations(Ott, n):
                            # 每一个i代表一个排列数，每一个排列数进行U计算
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                m10 = pow(1 + pow(j[k], 2), n * w[k])
                                m11 = pow(1 - pow(j[k], 2), n * w[k])
                                m1 = pow(m10 + 3 * m11, Q[k])
                                m2 = pow(m10 - m11, Q[k])
                                M1 *= m1
                                M2 *= m2
                            s1 = pow(M1 + 3 * M2, fn1)
                            s2 = pow(M1 - M2, fn1)
                            S1 *= s1
                            S2 *= s2
                        mols = 2 * pow(S1 - S2, lamda)
                        dels = pow(S1 + 3 * S2, lamda) + pow(S1 - S2, lamda)
                        u = math.sqrt(mols / dels)
                        hpu.append(u)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0ui), str(hp1ui)]
                n = len(HPs)
                # 排列数总数为n!个
                S1 = 1
                S2 = 1
                for i in itertools.permutations(Ott, n):
                    # 每一个i代表一个排列数，每一个排列数进行U计算
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        m10 = pow(1 + pow(j[k], 2), n * w[k])
                        m11 = pow(1 - pow(j[k], 2), n * w[k])
                        m1 = pow(m10 + 3 * m11, Q[k])
                        m2 = pow(m10 - m11, Q[k])
                        M1 *= m1
                        M2 *= m2
                    s1 = pow(M1 + 3 * M2, fn1)
                    s2 = pow(M1 - M2, fn1)
                    S1 *= s1
                    S2 *= s2
                mols = 2 * pow(S1 - S2, lamda)
                dels = pow(S1 + 3 * S2, lamda) + pow(S1 - S2, lamda)
                u = math.sqrt(mols / dels)
                hpu.append(u)
    # #计算V
    for hp0vi in V[0]:  #
        for hp1vi in V[1]:  #
            if (len(HPs) > 2):
                for hp2vi in V[2]:  #
                    if(len(HPs) > 3):
                        for hp3vi in V[3]:  #
                            if (len(HPs) > 4):
                                for hp4vi in V[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi), str(hp4vi)]
                                    n = len(HPs)
                                    S = 1
                                    for i in itertools.permutations(Ott, len(HPs)):
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            m10 = pow(2- pow(j[k], 2), n * w[k])
                                            m11 = pow(j[k], 2*n * w[k])
                                            m1 = pow(m10+3*m11,Q[k])
                                            m2 = pow(m10 - m11, Q[k])
                                            M1 *= m1
                                            M2 *= m2
                                        s1 = pow(M1+3*M2,fn1)
                                        s2 = pow(M1-M2,fn1)
                                        S1 *=s1
                                        S2 *=s2
                                    mols = pow(S1+3*S2,lamda)-pow(S1-S2,lamda)
                                    dels = pow(S1+3*S2,lamda)+pow(S1-S2,lamda)
                                    v = math.sqrt(mols/dels)
                                    hpv.append(v)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi)]
                                n = len(HPs)
                                S = 1
                                for i in itertools.permutations(Ott, len(HPs)):
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        m10 = pow(2 - pow(j[k], 2), n * w[k])
                                        m11 = pow(j[k], 2 * n * w[k])
                                        m1 = pow(m10 + 3 * m11, Q[k])
                                        m2 = pow(m10 - m11, Q[k])
                                        M1 *= m1
                                        M2 *= m2
                                    s1 = pow(M1 + 3 * M2, fn1)
                                    s2 = pow(M1 - M2, fn1)
                                    S1 *= s1
                                    S2 *= s2
                                mols = pow(S1 + 3 * S2, lamda) - pow(S1 - S2, lamda)
                                dels = pow(S1 + 3 * S2, lamda) + pow(S1 - S2, lamda)
                                v = math.sqrt(mols / dels)
                                hpv.append(v)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0vi), str(hp1vi), str(hp2vi)]
                        n = len(HPs)
                        S = 1
                        for i in itertools.permutations(Ott, len(HPs)):
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                m10 = pow(2 - pow(j[k], 2), n * w[k])
                                m11 = pow(j[k], 2 * n * w[k])
                                m1 = pow(m10 + 3 * m11, Q[k])
                                m2 = pow(m10 - m11, Q[k])
                                M1 *= m1
                                M2 *= m2
                            s1 = pow(M1 + 3 * M2, fn1)
                            s2 = pow(M1 - M2, fn1)
                            S1 *= s1
                            S2 *= s2
                        mols = pow(S1 + 3 * S2, lamda) - pow(S1 - S2, lamda)
                        dels = pow(S1 + 3 * S2, lamda) + pow(S1 - S2, lamda)
                        v = math.sqrt(mols / dels)
                        hpv.append(v)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0vi), str(hp1vi)]
                n = len(HPs)
                S = 1
                for i in itertools.permutations(Ott, len(HPs)):
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        m10 = pow(2 - pow(j[k], 2), n * w[k])
                        m11 = pow(j[k], 2 * n * w[k])
                        m1 = pow(m10 + 3 * m11, Q[k])
                        m2 = pow(m10 - m11, Q[k])
                        M1 *= m1
                        M2 *= m2
                    s1 = pow(M1 + 3 * M2, fn1)
                    s2 = pow(M1 - M2, fn1)
                    S1 *= s1
                    S2 *= s2
                mols = pow(S1 + 3 * S2, lamda) - pow(S1 - S2, lamda)
                dels = pow(S1 + 3 * S2, lamda) + pow(S1 - S2, lamda)
                v = math.sqrt(mols / dels)
                hpv.append(v)
    UV.append([hpu,hpv])
    return UV # 返回这个计算后的结果

def WPHFAMM_Frank_5(HPs,w,Q,Theta):#Mλp,#HPs表示多个PHFNs，是一个多维数组，每一个元素代表一个HP
    UV = []
    hpu = []
    hpv = []
    #HPs是一个数组，有多个元素，每一个元素代表一个PH,假设HPs中有不超过5个HP
    U = []
    V = []
    eta = Theta-1
    qq = 0
    for qj in Q:
        qq += qj
    lamda = float(1 / qq)
    n = len(HPs)
    fn1 = 1 / factorial(n)
    for i in range(len(HPs)): #将UV分离出来，分别进行处理
        ui=[]
        vi=[]
        for j in range(len(HPs[i])):
            ui.append(HPs[i][j][0])
            vi.append(HPs[i][j][1])
        U.append(ui)
        V.append(vi)
    for hp0ui in U[0]: #
        for hp1ui in U[1]:  #
            if (len(HPs)>2):
                for hp2ui in U[2]:  #
                    if(len(HPs)>3):
                        for hp3ui in U[3]:  #
                            if (len(HPs) > 4):
                                for hp4ui in U[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui), str(hp4ui)]
                                    n = len(HPs)
                                    # 排列数总数为n!个
                                    S1 = 1
                                    S2 = 1
                                    for i in itertools.permutations(Ott, n):
                                        # 每一个i代表一个排列数，每一个排列数进行U计算
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            m11 = pow(pow(Theta,1-pow(j[k],2))-1,n*w[k])
                                            m1 = pow(pow(eta,2)*m11+pow(eta,n*w[k]+1),Q[k])
                                            m2 = pow(pow(eta,n*w[k]+1)-eta*m11,Q[k])
                                            M1 *= m1
                                            M2 *= m2
                                        s1 = pow(eta*(M1-M2),fn1)
                                        s2 = pow(pow(eta,2)*M2+eta*M1,fn1)
                                        S1 *= s1
                                        S2 *= s2
                                    Au = pow(pow(eta,2)*S1+eta*S2,lamda)
                                    Bu = pow(eta*(S2-S1),lamda)
                                    hh = (eta*Bu+Au)/ Au
                                    u = math.sqrt(math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                                    hpu.append(u)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui)]
                                n = len(HPs)
                                # 排列数总数为n!个
                                S1 = 1
                                S2 = 1
                                for i in itertools.permutations(Ott, n):
                                    # 每一个i代表一个排列数，每一个排列数进行U计算
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        m11 = pow(pow(Theta, 1 - pow(j[k], 2)) - 1, n * w[k])
                                        m1 = pow(pow(eta, 2) * m11 + pow(eta, n * w[k] + 1), Q[k])
                                        m2 = pow(pow(eta, n * w[k] + 1) - eta * m11, Q[k])
                                        M1 *= m1
                                        M2 *= m2
                                    s1 = pow(eta * (M1 - M2), fn1)
                                    s2 = pow(pow(eta, 2) * M2 + eta * M1, fn1)
                                    S1 *= s1
                                    S2 *= s2
                                Au = pow(pow(eta, 2) * S1 + eta * S2, lamda)
                                Bu = pow(eta * (S2 - S1), lamda)
                                hh = (eta * Bu + Au) / Au
                                u = math.sqrt(math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                                hpu.append(u)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0ui), str(hp1ui), str(hp2ui)]
                        n = len(HPs)
                        # 排列数总数为n!个
                        S1 = 1
                        S2 = 1
                        for i in itertools.permutations(Ott, n):
                            # 每一个i代表一个排列数，每一个排列数进行U计算
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                m11 = pow(pow(Theta, 1 - pow(j[k], 2)) - 1, n * w[k])
                                m1 = pow(pow(eta, 2) * m11 + pow(eta, n * w[k] + 1), Q[k])
                                m2 = pow(pow(eta, n * w[k] + 1) - eta * m11, Q[k])
                                M1 *= m1
                                M2 *= m2
                            s1 = pow(eta * (M1 - M2), fn1)
                            s2 = pow(pow(eta, 2) * M2 + eta * M1, fn1)
                            S1 *= s1
                            S2 *= s2
                        Au = pow(pow(eta, 2) * S1 + eta * S2, lamda)
                        Bu = pow(eta * (S2 - S1), lamda)
                        hh = (eta * Bu + Au) / Au
                        u = math.sqrt(math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                        hpu.append(u)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0ui), str(hp1ui)]
                n = len(HPs)
                # 排列数总数为n!个
                S1 = 1
                S2 = 1
                for i in itertools.permutations(Ott, n):
                    # 每一个i代表一个排列数，每一个排列数进行U计算
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        m11 = pow(pow(Theta, 1 - pow(j[k], 2)) - 1, n * w[k])
                        m1 = pow(pow(eta, 2) * m11 + pow(eta, n * w[k] + 1), Q[k])
                        m2 = pow(pow(eta, n * w[k] + 1) - eta * m11, Q[k])
                        M1 *= m1
                        M2 *= m2
                    s1 = pow(eta * (M1 - M2), fn1)
                    s2 = pow(pow(eta, 2) * M2 + eta * M1, fn1)
                    S1 *= s1
                    S2 *= s2
                Au = pow(pow(eta, 2) * S1 + eta * S2, lamda)
                Bu = pow(eta * (S2 - S1), lamda)
                hh = (eta * Bu + Au) / Au
                u = math.sqrt(math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                hpu.append(u)
    # #计算V
    for hp0vi in V[0]:  #
        for hp1vi in V[1]:  #
            if (len(HPs) > 2):
                for hp2vi in V[2]:  #
                    if(len(HPs) > 3):
                        for hp3vi in V[3]:  #
                            if (len(HPs) > 4):
                                for hp4vi in V[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi), str(hp4vi)]
                                    n = len(HPs)
                                    # 排列数总数为n!个
                                    S1 = 1
                                    S2 = 1
                                    for i in itertools.permutations(Ott, n):
                                        # 每一个i代表一个排列数，每一个排列数进行U计算
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            m11 = pow(pow(Theta, pow(j[k], 2)) - 1, n * w[k])
                                            m1 = pow(pow(eta, 2) * m11 + pow(eta, n * w[k] + 1), Q[k])
                                            m2 = pow(pow(eta, n * w[k] + 1) - eta * m11, Q[k])
                                            M1 *= m1
                                            M2 *= m2
                                        s1 = pow(eta * (M1 - M2), fn1)
                                        s2 = pow(pow(eta, 2) * M2 + eta * M1, fn1)
                                        S1 *= s1
                                        S2 *= s2
                                    Av = pow(pow(eta, 2) * S1 + eta * S2, lamda)
                                    Bv = pow(eta * (S2 - S1), lamda)
                                    hh = (eta * Bv + Av) / Av
                                    v = math.sqrt(1-math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                                    hpv.append(v)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi)]
                                n = len(HPs)
                                # 排列数总数为n!个
                                S1 = 1
                                S2 = 1
                                for i in itertools.permutations(Ott, n):
                                    # 每一个i代表一个排列数，每一个排列数进行U计算
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        m11 = pow(pow(Theta, pow(j[k], 2)) - 1, n * w[k])
                                        m1 = pow(pow(eta, 2) * m11 + pow(eta, n * w[k] + 1), Q[k])
                                        m2 = pow(pow(eta, n * w[k] + 1) - eta * m11, Q[k])
                                        M1 *= m1
                                        M2 *= m2
                                    s1 = pow(eta * (M1 - M2), fn1)
                                    s2 = pow(pow(eta, 2) * M2 + eta * M1, fn1)
                                    S1 *= s1
                                    S2 *= s2
                                Av = pow(pow(eta, 2) * S1 + eta * S2, lamda)
                                Bv = pow(eta * (S2 - S1), lamda)
                                hh = (eta * Bv + Av) / Av
                                v = math.sqrt(1 - math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                                hpv.append(v)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0vi), str(hp1vi), str(hp2vi)]
                        n = len(HPs)
                        # 排列数总数为n!个
                        S1 = 1
                        S2 = 1
                        for i in itertools.permutations(Ott, n):
                            # 每一个i代表一个排列数，每一个排列数进行U计算
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                m11 = pow(pow(Theta, pow(j[k], 2)) - 1, n * w[k])
                                m1 = pow(pow(eta, 2) * m11 + pow(eta, n * w[k] + 1), Q[k])
                                m2 = pow(pow(eta, n * w[k] + 1) - eta * m11, Q[k])
                                M1 *= m1
                                M2 *= m2
                            s1 = pow(eta * (M1 - M2), fn1)
                            s2 = pow(pow(eta, 2) * M2 + eta * M1, fn1)
                            S1 *= s1
                            S2 *= s2
                        Av = pow(pow(eta, 2) * S1 + eta * S2, lamda)
                        Bv = pow(eta * (S2 - S1), lamda)
                        hh = (eta * Bv + Av) / Av
                        v = math.sqrt(1 - math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                        hpv.append(v)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0vi), str(hp1vi)]
                n = len(HPs)
                # 排列数总数为n!个
                S1 = 1
                S2 = 1
                for i in itertools.permutations(Ott, n):
                    # 每一个i代表一个排列数，每一个排列数进行U计算
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        m11 = pow(pow(Theta, pow(j[k], 2)) - 1, n * w[k])
                        m1 = pow(pow(eta, 2) * m11 + pow(eta, n * w[k] + 1), Q[k])
                        m2 = pow(pow(eta, n * w[k] + 1) - eta * m11, Q[k])
                        M1 *= m1
                        M2 *= m2
                    s1 = pow(eta * (M1 - M2), fn1)
                    s2 = pow(pow(eta, 2) * M2 + eta * M1, fn1)
                    S1 *= s1
                    S2 *= s2
                Av = pow(pow(eta, 2) * S1 + eta * S2, lamda)
                Bv = pow(eta * (S2 - S1), lamda)
                hh = (eta * Bv + Av) / Av
                v = math.sqrt(1 - math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                hpv.append(v)
    UV.append([hpu,hpv])
    return UV # 返回这个计算后的结果

def PHFAMM_Frank_5(HPs,w,Q,Theta):#Mλp,#HPs表示多个PHFNs，是一个多维数组，每一个元素代表一个HP
    UV = []
    hpu = []
    hpv = []
    #HPs是一个数组，有多个元素，每一个元素代表一个PH,假设HPs中有不超过5个HP
    U = []
    V = []
    eta = Theta-1
    qq = 0
    for qj in Q:
        qq += qj
    lamda = float(1 / qq)
    n = len(HPs)
    fn1 = 1 / factorial(n)
    for i in range(len(HPs)): #将UV分离出来，分别进行处理
        ui=[]
        vi=[]
        for j in range(len(HPs[i])):
            ui.append(HPs[i][j][0])
            vi.append(HPs[i][j][1])
        U.append(ui)
        V.append(vi)
    for hp0ui in U[0]: #
        for hp1ui in U[1]:  #
            if (len(HPs)>2):
                for hp2ui in U[2]:  #
                    if(len(HPs)>3):
                        for hp3ui in U[3]:  #
                            if (len(HPs) > 4):
                                for hp4ui in U[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui), str(hp4ui)]
                                    n = len(HPs)
                                    # 排列数总数为n!个
                                    S1 = 1
                                    S2 = 1
                                    for i in itertools.permutations(Ott, n):
                                        # 每一个i代表一个排列数，每一个排列数进行U计算
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            M1 *= pow(pow(Theta,pow(j[k],2))-1,Q[k])
                                            M2 *= pow(eta,Q[k])
                                        S1 *= pow(pow(eta,2)*M1+eta*M2,fn1)
                                        S2 *= pow(eta*(M2-M1),fn1)
                                    ro1 = pow(eta*(S1-S2),lamda)
                                    ro2 = pow(pow(eta,2)*S2+eta*S1,lamda)
                                    hh = (eta*ro1+ro2)/ro2
                                    u = math.sqrt(math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                                    hpu.append(u)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0ui), str(hp1ui), str(hp2ui), str(hp3ui)]
                                n = len(HPs)
                                # 排列数总数为n!个
                                S1 = 1
                                S2 = 1
                                for i in itertools.permutations(Ott, n):
                                    # 每一个i代表一个排列数，每一个排列数进行U计算
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        M1 *= pow(pow(Theta, pow(j[k], 2)) - 1, Q[k])
                                        M2 *= pow(eta, Q[k])
                                    S1 *= pow(pow(eta, 2) * M1 + eta * M2, fn1)
                                    S2 *= pow(eta * (M2 - M1), fn1)
                                ro1 = pow(eta * (S1 - S2), lamda)
                                ro2 = pow(pow(eta, 2) * S2 + eta * S1, lamda)
                                hh = (eta * ro1 + ro2) / ro2
                                u = math.sqrt(math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                                hpu.append(u)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0ui), str(hp1ui), str(hp2ui)]
                        n = len(HPs)
                        # 排列数总数为n!个
                        S1 = 1
                        S2 = 1
                        for i in itertools.permutations(Ott, n):
                            # 每一个i代表一个排列数，每一个排列数进行U计算
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                M1 *= pow(pow(Theta, pow(j[k], 2)) - 1, Q[k])
                                M2 *= pow(eta, Q[k])
                            S1 *= pow(pow(eta, 2) * M1 + eta * M2, fn1)
                            S2 *= pow(eta * (M2 - M1), fn1)
                        ro1 = pow(eta * (S1 - S2), lamda)
                        ro2 = pow(pow(eta, 2) * S2 + eta * S1, lamda)
                        hh = (eta * ro1 + ro2) / ro2
                        u = math.sqrt(math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                        hpu.append(u)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0ui), str(hp1ui)]
                n = len(HPs)
                # 排列数总数为n!个
                S1 = 1
                S2 = 1
                for i in itertools.permutations(Ott, n):
                    # 每一个i代表一个排列数，每一个排列数进行U计算
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        M1 *= pow(pow(Theta, pow(j[k], 2)) - 1, Q[k])
                        M2 *= pow(eta, Q[k])
                    S1 *= pow(pow(eta, 2) * M1 + eta * M2, fn1)
                    S2 *= pow(eta * (M2 - M1), fn1)
                ro1 = pow(eta * (S1 - S2), lamda)
                ro2 = pow(pow(eta, 2) * S2 + eta * S1, lamda)
                hh = (eta * ro1 + ro2) / ro2
                u = math.sqrt(math.log(hh, Theta))  # log(b,a):表示以a底b的对数
                hpu.append(u)
    # #计算V
    for hp0vi in V[0]:  #
        for hp1vi in V[1]:  #
            if (len(HPs) > 2):
                for hp2vi in V[2]:  #
                    if(len(HPs) > 3):
                        for hp3vi in V[3]:  #
                            if (len(HPs) > 4):
                                for hp4vi in V[4]:  #
                                    # 在将每一个排列做运算后再综合运算，得出结果
                                    Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi), str(hp4vi)]
                                    S1 = 1
                                    S2 = 1
                                    for i in itertools.permutations(Ott, len(HPs)):
                                        j = []
                                        for ii in i:
                                            j.append(float(ii))
                                        M1 = 1
                                        M2 = 1
                                        for k in range(0, len(j)):
                                            M1 *= pow(pow(Theta,1-pow(j[k],2))-1,Q[k])
                                            M2 *= pow(eta,Q[k])
                                        S1 *= pow(pow(eta,2)*M1+eta*M2,fn1)
                                        S2 *= pow(eta*(M2-M1),fn1)
                                    ro1 = pow(eta*(S1-S2),lamda)
                                    ro2 = pow(pow(eta,2)*S2+eta*S1,lamda)
                                    hh = (eta*ro1+ro2)/ro2
                                    v = math.sqrt(1-math.log(hh,Theta))
                                    hpv.append(v)
                            else:
                                # 在将每一个排列做运算后再综合运算，得出结果
                                Ott = [str(hp0vi), str(hp1vi), str(hp2vi), str(hp3vi)]
                                S1 = 1
                                S2 = 1
                                for i in itertools.permutations(Ott, len(HPs)):
                                    j = []
                                    for ii in i:
                                        j.append(float(ii))
                                    M1 = 1
                                    M2 = 1
                                    for k in range(0, len(j)):
                                        M1 *= pow(pow(Theta, 1 - pow(j[k], 2)) - 1, Q[k])
                                        M2 *= pow(eta, Q[k])
                                    S1 *= pow(pow(eta, 2) * M1 + eta * M2, fn1)
                                    S2 *= pow(eta * (M2 - M1), fn1)
                                ro1 = pow(eta * (S1 - S2), lamda)
                                ro2 = pow(pow(eta, 2) * S2 + eta * S1, lamda)
                                hh = (eta * ro1 + ro2) / ro2
                                v = math.sqrt(1 - math.log(hh, Theta))
                                hpv.append(v)
                    else:
                        # 在将每一个排列做运算后再综合运算，得出结果
                        Ott = [str(hp0vi), str(hp1vi), str(hp2vi)]
                        S1 = 1
                        S2 = 1
                        for i in itertools.permutations(Ott, len(HPs)):
                            j = []
                            for ii in i:
                                j.append(float(ii))
                            M1 = 1
                            M2 = 1
                            for k in range(0, len(j)):
                                M1 *= pow(pow(Theta, 1 - pow(j[k], 2)) - 1, Q[k])
                                M2 *= pow(eta, Q[k])
                            S1 *= pow(pow(eta, 2) * M1 + eta * M2, fn1)
                            S2 *= pow(eta * (M2 - M1), fn1)
                        ro1 = pow(eta * (S1 - S2), lamda)
                        ro2 = pow(pow(eta, 2) * S2 + eta * S1, lamda)
                        hh = (eta * ro1 + ro2) / ro2
                        v = math.sqrt(1 - math.log(hh, Theta))
                        hpv.append(v)
            else:
                # 在将每一个排列做运算后再综合运算，得出结果
                Ott = [str(hp0vi), str(hp1vi)]
                S1 = 1
                S2 = 1
                for i in itertools.permutations(Ott, len(HPs)):
                    j = []
                    for ii in i:
                        j.append(float(ii))
                    M1 = 1
                    M2 = 1
                    for k in range(0, len(j)):
                        M1 *= pow(pow(Theta, 1 - pow(j[k], 2)) - 1, Q[k])
                        M2 *= pow(eta, Q[k])
                    S1 *= pow(pow(eta, 2) * M1 + eta * M2, fn1)
                    S2 *= pow(eta * (M2 - M1), fn1)
                ro1 = pow(eta * (S1 - S2), lamda)
                ro2 = pow(pow(eta, 2) * S2 + eta * S1, lamda)
                hh = (eta * ro1 + ro2) / ro2
                v = math.sqrt(1 - math.log(hh, Theta))
                hpv.append(v)
    UV.append([hpu,hpv])
    return UV # 返回这个计算后的结果

def testPFMM():
    #使用此算法来测试论文《PFMM》中的数据，结果是一致的：A2 > A5 > A3 > A4 > A1
    scores = []
    PHF1s = [
        [[0.4160,0.8421]],
        [[0.7319,0.6000]],
        [[0.5313,0.7718]],
        [[0.3557,0.6656]]
    ]
    w = [0.2, 0.1, 0.3, 0.4]
    Q = [1, 1, 1, 1]
    #Q = [0.5, 0.5, 0.5, 0.5]
    result1 = WPHFAMM_Hamacher_5(PHF1s, w, Q,2)

    score1 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result1))
    scores.append(score1)
    PHF2s = [
        [[0.6649,0.4702]],
        [[0.8653,0.2860]],
        [[0.8000,0.3123]],
        [[0.3000,0.5372]]
    ]
    result2 = WPHFAMM_Hamacher_5(PHF2s, w, Q, 2)
    score2 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result2))
    scores.append(score2)
    PHF3s = [
        [[0.3302,0.5919]],
        [[0.5278,0.6481]],
        [[0.6649,0.4354]],
        [[0.3634,0.6987]]
    ]
    result3 = WPHFAMM_Hamacher_5(PHF3s, w, Q, 2)
    score3 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result3))
    scores.append(score3)
    PHF4s = [
        [[0.3634,0.8047]],
        [[0.6257,0.5372]],
        [[0.6316,0.2000]],
        [[0.2520,0.6377]]
    ]
    result4 = WPHFAMM_Hamacher_5(PHF4s, w, Q, 2)
    score4 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result4))
    scores.append(score4)
    PHF5s = [
        [[0.2154,0.6707]],
        [[0.7268,0.2860]],
        [[0.8653,0.3131]],
        [[0.6000,0.6119]]
    ]
    result5 = WPHFAMM_Hamacher_5(PHF5s, w, Q, 2)
    score5 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result5))
    scores.append(score5)

    for i in range(0, len(scores)):
        print("The score of A%d is: %f" % (i + 1, scores[i]))

def testArc_in_PHF():
    #测试论文《Arc in PHF》中的数据，结果几乎一致：Z3 > Z4 > Z5 > Z2 > Z1
    w = [0.4, 0.1, 0.2, 0.3]
    Q = [1, 1, 1, 1]
    Theta = 3
    scores = []
    PHF1s = [[[0.2, 0.3], [0.3, 0.4]],#第一个PHFNs
             [[0.4, 0.6], [0.5, 0.4], [0.7, 0.2]],#第2个PHFNs
             [[0.4, 0.5], [0.6, 0.3]],#第3个PHFNs
             [[0.6, 0.3], [0.7, 0.4]]]#第4个PHFNs

    result1 = WPHFAMM_Hamacher_5(PHF1s, w, Q,Theta)
    score1 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result1))
    scores.append(score1)
    PHF2s = [[[0.4, 0.3], [0.6, 0.4]],
             [[0.5, 0.6], [0.6, 0.4]],
             [[0.5, 0.3], [0.5, 0.6]],
             [[0.5, 0.4], [0.7, 0.6]]]
    result2 = WPHFAMM_Hamacher_5(PHF2s, w, Q, Theta)
    score2 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result2))
    scores.append(score2)
    PHF3s = [[[0.6, 0.2], [0.7, 0.3]],
             [[0.5, 0.3], [0.5, 0.4]],
             [[0.5, 0.2], [0.8, 0.6], [0.8, 0.2]],
             [[0.4, 0.5], [0.6, 0.4]]]
    result3 = WPHFAMM_Hamacher_5(PHF3s, w, Q, Theta)
    score3 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result3))
    scores.append(score3)
    PHF4s = [[[0.6, 0.5], [0.7, 0.4]],
             [[0.5, 0.2], [0.6, 0.5]],
             [[0.4, 0.3], [0.5, 0.4]],
             [[0.6, 0.2], [0.6, 0.3], [0.8, 0.4]]]
    result4 = WPHFAMM_Hamacher_5(PHF4s, w, Q, Theta)
    score4 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result4))
    scores.append(score4)
    PHF5s = [[[0.3, 0.3], [0.6, 0.4]],
             [[0.5, 0.4], [0.7, 0.4]],
             [[0.6, 0.4], [0.7, 0.4]],
             [[0.4, 0.6], [0.5, 0.3]]]
    result5 = WPHFAMM_Hamacher_5(PHF5s, w, Q, Theta)
    score5 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result5))
    scores.append(score5)

    for i in range(0, len(scores)):
        print("The score of Z%d is: %f" % (i + 1, scores[i]))

def testHamacher_in_PHF():
    #测试论文《Hamacher in PHF》中的数据(与论文《Arc in PHF》中的数据一样)，结果几乎一致：Z3 > Z4 > Z5 > Z2 > Z1
    scores = []
    w = [0.4, 0.1, 0.2, 0.3]
    Q = [1, 1, 1, 1]
    Theta = 1 #(无论如何变化，不影响等级排名)
    PHF1s = [[[0.2, 0.3], [0.3, 0.4]],#第一个PHFNs
             [[0.4, 0.6], [0.5, 0.4], [0.7, 0.2]],#第2个PHFNs
             [[0.4, 0.5], [0.6, 0.3]],#第3个PHFNs
             [[0.6, 0.3], [0.7, 0.4]]]#第4个PHFNs
    result1 = WPHFAMM_Aglebraic_5(PHF1s, w, Q, Theta)
    score1 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result1))
    scores.append(score1)
    PHF2s = [[[0.4, 0.3], [0.6, 0.4]],
             [[0.5, 0.6], [0.6, 0.4]],
             [[0.5, 0.3], [0.5, 0.6]],
             [[0.5, 0.4], [0.7, 0.6]]]
    result2 = WPHFAMM_Aglebraic_5(PHF2s, w, Q, Theta)
    score2 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result2))
    scores.append(score2)
    PHF3s = [[[0.6, 0.2], [0.7, 0.3]],
             [[0.5, 0.3], [0.5, 0.4]],
             [[0.5, 0.2], [0.8, 0.6], [0.8, 0.2]],
             [[0.4, 0.5], [0.6, 0.4]]]
    result3 = WPHFAMM_Aglebraic_5(PHF3s, w, Q, Theta)
    score3 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result3))
    scores.append(score3)
    PHF4s = [[[0.6, 0.5], [0.7, 0.4]],
             [[0.5, 0.2], [0.6, 0.5]],
             [[0.4, 0.3], [0.5, 0.4]],
             [[0.6, 0.2], [0.6, 0.3], [0.8, 0.4]]]
    result4 = WPHFAMM_Aglebraic_5(PHF4s, w, Q, Theta)
    score4 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result4))
    scores.append(score4)
    PHF5s = [[[0.3, 0.3], [0.6, 0.4]],
             [[0.5, 0.4], [0.7, 0.4]],
             [[0.6, 0.4], [0.7, 0.4]],
             [[0.4, 0.6], [0.5, 0.3]]]
    result5 = WPHFAMM_Aglebraic_5(PHF5s, w, Q, Theta)
    score5 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result5))
    scores.append(score5)

    for i in range(0, len(scores)):
        print("The score of Z%d is: %f" % (i + 1, scores[i]))

def DealWith_PHFAAWMM():
    print("MMM")
    #testArc_in_PHF()
    testPFMM()
    #testHamacher_in_PHF()

if __name__ == "__main__":
    DealWith_PHFAAWMM()