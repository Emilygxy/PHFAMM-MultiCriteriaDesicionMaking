import itertools
import math
import numpy as np
import matplotlib.pyplot as plt

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

def draw(Thetas,scores):
    x = np.array(Thetas)  # 将a转化为numpy的array
    y = np.array(scores)  # 将a转化为numpy的array
    # 创建一幅线图，x轴是年份，y轴是gdb
    plt.plot(x, y, color='green', marker='o',linestyle='solid')
    # 添加一个标题
    plt.title('Scores for alternatives(WPHFAFMM)')
    # 给y轴加标记
    plt.ylabel('y = Scores ')
    plt.xlabel('x = Theta ')
    # 输出线图
    plt.show()

def draw_Multi_F(Thetas,Mscores):
    x = np.array(Thetas)  # 将a转化为numpy的array

    # 创建一幅线图，x轴是年份，y轴是gdb
    colors = ['blue','red','green','purple','black']
    lables = ['Z1', 'Z2', 'Z3', 'Z4', 'Z5']
    n = 0
    for yi in Mscores:
        y = np.array(yi)  # 将a转化为numpy的array
        plt.plot(x,y , color=colors[n], linewidth=1.0, marker='.',linestyle='-')
        n = n + 1
    # 添加一个标题
    plt.title('Scores for alternatives(WPHFAFMM)')
    # 给y轴加标记
    plt.ylabel('y = Scores ')
    plt.xlabel('x = Theta ')
    plt.legend(lables,loc='upper right')#loc(设置图例显示的位置)
    # 输出线图
    plt.show()

def draw_Multi_H(Thetas,Mscores):
    x = np.array(Thetas)  # 将a转化为numpy的array

    # 创建一幅线图，x轴是年份，y轴是gdb
    colors = ['blue','red','green','purple','black']
    lables = ['Z1', 'Z2', 'Z3', 'Z4', 'Z5']
    n = 0
    for yi in Mscores:
        y = np.array(yi)  # 将a转化为numpy的array
        plt.plot(x,y , color=colors[n], linewidth=1.0, marker='.',linestyle='-')
        n = n + 1
    # 添加一个标题
    plt.title('Scores for alternatives(WPHFAHMM)')
    # 给y轴加标记
    plt.ylabel('y = Scores ')
    plt.xlabel('x = Theta ')
    plt.legend(lables,loc='upper right')#loc(设置图例显示的位置)
    # 输出线图
    plt.show()

def testArc_in_PHF_MultiTheta_Hamacher():
    w = [0.4, 0.1, 0.2, 0.3]
    Q = [1, 1, 1, 1]
    Sscores = []
    PHF1s = [[[0.2, 0.3], [0.3, 0.4]],#第一个PHFNs
             [[0.4, 0.6], [0.5, 0.4], [0.7, 0.2]],#第2个PHFNs
             [[0.4, 0.5], [0.6, 0.3]],#第3个PHFNs
             [[0.6, 0.3], [0.7, 0.4]]]#第4个PHFNs
    scores1 = []
    Thetas = []
    for k in range(1,32,1):
        result1 = WPHFAMM_Hamacher_5(PHF1s, w, Q,k)
        score1 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result1))
        scores1.append(score1)
        Thetas.append(k)
    Sscores.append(scores1)
    PHF2s = [[[0.4, 0.3], [0.6, 0.4]],
             [[0.5, 0.6], [0.6, 0.4]],
             [[0.5, 0.3], [0.5, 0.6]],
             [[0.5, 0.4], [0.7, 0.6]]]
    scores2 = []
    for k in range(1,32,1):
        result2 = WPHFAMM_Hamacher_5(PHF2s, w, Q,k)
        score2 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result2))
        scores2.append(score2)
    Sscores.append(scores2)
    PHF3s = [[[0.6, 0.2], [0.7, 0.3]],
             [[0.5, 0.3], [0.5, 0.4]],
             [[0.5, 0.2], [0.8, 0.6], [0.8, 0.2]],
             [[0.4, 0.5], [0.6, 0.4]]]
    scores3 = []
    for k in range(1, 32, 1):
        result3 = WPHFAMM_Hamacher_5(PHF3s, w, Q, k)
        score3 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result3))
        scores3.append(score3)
    Sscores.append(scores3)
    PHF4s = [[[0.6, 0.5], [0.7, 0.4]],
             [[0.5, 0.2], [0.6, 0.5]],
             [[0.4, 0.3], [0.5, 0.4]],
             [[0.6, 0.2], [0.6, 0.3], [0.8, 0.4]]]
    scores4 = []
    for k in range(1, 32, 1):
        result4 = WPHFAMM_Hamacher_5(PHF4s, w, Q, k)
        score4 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result4))
        scores4.append(score4)
    Sscores.append(scores4)
    PHF5s = [[[0.3, 0.3], [0.6, 0.4]],
             [[0.5, 0.4], [0.7, 0.4]],
             [[0.6, 0.4], [0.7, 0.4]],
             [[0.4, 0.6], [0.5, 0.3]]]
    scores5 = []
    for k in range(1, 32, 1):
        result5 = WPHFAMM_Hamacher_5(PHF5s, w, Q, k)
        score5 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result5))
        scores5.append(score5)
    Sscores.append(scores5)
    draw_Multi_H(Thetas,Sscores)
    # 测试论文《Arc in PHF》中的数据，结果几乎一致：Z3 > Z4 > Z5 > Z2 > Z1

def testArc_in_PHF_MultiTheta_Frank():
    w = [0.4, 0.1, 0.2, 0.3]
    Q = [1, 1, 1, 1]
    Sscores = []
    PHF1s = [[[0.2, 0.3], [0.3, 0.4]],#第一个PHFNs
             [[0.4, 0.6], [0.5, 0.4], [0.7, 0.2]],#第2个PHFNs
             [[0.4, 0.5], [0.6, 0.3]],#第3个PHFNs
             [[0.6, 0.3], [0.7, 0.4]]]#第4个PHFNs
    scores1 = []
    Thetas = []
    for k in range(2,32,1):
        result1 = WPHFAMM_Frank_5(PHF1s, w, Q,k)
        score1 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result1))
        scores1.append(score1)
        Thetas.append(k)
    Sscores.append(scores1)
    PHF2s = [[[0.4, 0.3], [0.6, 0.4]],
             [[0.5, 0.6], [0.6, 0.4]],
             [[0.5, 0.3], [0.5, 0.6]],
             [[0.5, 0.4], [0.7, 0.6]]]
    scores2 = []
    for k in range(2,32,1):
        result2 = WPHFAMM_Frank_5(PHF2s, w, Q,k)
        score2 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result2))
        scores2.append(score2)
    Sscores.append(scores2)
    PHF3s = [[[0.6, 0.2], [0.7, 0.3]],
             [[0.5, 0.3], [0.5, 0.4]],
             [[0.5, 0.2], [0.8, 0.6], [0.8, 0.2]],
             [[0.4, 0.5], [0.6, 0.4]]]
    scores3 = []
    for k in range(2, 32, 1):
        result3 = WPHFAMM_Frank_5(PHF3s, w, Q, k)
        score3 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result3))
        scores3.append(score3)
    Sscores.append(scores3)
    PHF4s = [[[0.6, 0.5], [0.7, 0.4]],
             [[0.5, 0.2], [0.6, 0.5]],
             [[0.4, 0.3], [0.5, 0.4]],
             [[0.6, 0.2], [0.6, 0.3], [0.8, 0.4]]]
    scores4 = []
    for k in range(2, 32, 1):
        result4 = WPHFAMM_Frank_5(PHF4s, w, Q, k)
        score4 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result4))
        scores4.append(score4)
    Sscores.append(scores4)
    PHF5s = [[[0.3, 0.3], [0.6, 0.4]],
             [[0.5, 0.4], [0.7, 0.4]],
             [[0.6, 0.4], [0.7, 0.4]],
             [[0.4, 0.6], [0.5, 0.3]]]
    scores5 = []
    for k in range(2, 32, 1):
        result5 = WPHFAMM_Frank_5(PHF5s, w, Q, k)
        score5 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result5))
        scores5.append(score5)
    Sscores.append(scores5)
    draw_Multi_F(Thetas,Sscores)
    # 测试论文《Arc in PHF》中的数据，结果几乎一致：Z3 > Z4 > Z5 > Z2 > Z1

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

    result1 = WPHFAMM_Frank_5(PHF1s, w, Q,Theta)
    score1 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result1))
    scores.append(score1)
    PHF2s = [[[0.4, 0.3], [0.6, 0.4]],
             [[0.5, 0.6], [0.6, 0.4]],
             [[0.5, 0.3], [0.5, 0.6]],
             [[0.5, 0.4], [0.7, 0.6]]]
    result2 = WPHFAMM_Frank_5(PHF2s, w, Q, Theta)
    score2 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result2))
    scores.append(score2)
    PHF3s = [[[0.6, 0.2], [0.7, 0.3]],
             [[0.5, 0.3], [0.5, 0.4]],
             [[0.5, 0.2], [0.8, 0.6], [0.8, 0.2]],
             [[0.4, 0.5], [0.6, 0.4]]]
    result3 = WPHFAMM_Frank_5(PHF3s, w, Q, Theta)
    score3 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result3))
    scores.append(score3)
    PHF4s = [[[0.6, 0.5], [0.7, 0.4]],
             [[0.5, 0.2], [0.6, 0.5]],
             [[0.4, 0.3], [0.5, 0.4]],
             [[0.6, 0.2], [0.6, 0.3], [0.8, 0.4]]]
    result4 = WPHFAMM_Frank_5(PHF4s, w, Q, Theta)
    score4 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result4))
    scores.append(score4)
    PHF5s = [[[0.3, 0.3], [0.6, 0.4]],
             [[0.5, 0.4], [0.7, 0.4]],
             [[0.6, 0.4], [0.7, 0.4]],
             [[0.4, 0.6], [0.5, 0.3]]]
    result5 = WPHFAMM_Frank_5(PHF5s, w, Q, Theta)
    score5 = getScore_PHFAWMM(getAccuracy_PHFAWMM(result5))
    scores.append(score5)

    for i in range(0, len(scores)):
        print("The score of Z%d is: %f" % (i + 1, scores[i]))

def DealWith_WPHFAWMM():
    #testArc_in_PHF_MultiTheta_Frank()
    testArc_in_PHF_MultiTheta_Hamacher()

if __name__ == "__main__":
    DealWith_WPHFAWMM()