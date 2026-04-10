# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.linalg import cholesky
import time

def MoM_old_with_additive(Z, y, nmc=40):
    n = Z.shape[0]
    m = Z.shape[1]
    p = m*(m-1)/2
    u = np.random.randn(n,nmc)
    ZZ = Z*Z
    
    ### three important elSements
    Ku = Z @ (Z.T @ u) / m
    Du = ZZ @ (ZZ.T @ u) 
    KKu = np.zeros((n, nmc))
    for i in range(nmc):
        Z_tilde = u[:, i].reshape(-1, 1) * Z   
        M = Z.T @ Z_tilde
        ZM = Z @ M
        KKu[:, i] = np.sum(ZM * Z, axis=1)

    
    ############Trace part
    ###trK_app
    trK_app = np.mean(np.sum(u * Ku, axis=0))
    ###trW_app
    trW_app = 0.5*(np.mean(np.sum(u * KKu, axis=0)) - np.mean(np.sum(u * Du, axis=0)))/p
    ###trK2_app
    trK2_app = np.mean(np.sum(Ku * Ku, axis=0))
    ###trW2_app
    trKK2 =np.mean(np.sum(KKu *KKu, axis=0))
    trDK =np.mean(np.sum(Du *KKu, axis=0))
    trDD = np.mean(np.sum(Du * Du, axis=0))
    trW2_app = (trKK2 - 2*trDK +trDD)/(4*p**2)
    ###trWK_app
    trWK_app = 0.5*(np.mean(np.sum(Ku * KKu, axis=0)) - np.mean(np.sum(Ku * Du, axis=0)))/p

    
    ############# Quadratic forms part
    yTy = y @ y
    yTWy = 0.5*(y.T@(np.sum((Z @(Z.T@(y.reshape(-1, 1) * Z)))*Z,axis = 1)) - (y @ (ZZ @ (ZZ.T @ y))))/p
    yTKy =(y @ (Z @ (Z.T @ y))) / m

 

    A = np.array([
        [trK2_app, trWK_app, trK_app],
        [trWK_app, trW2_app, trW_app],
        [trK_app,  trW_app,  n]
    ])
    
    
    b = np.array([
        yTKy,
        yTWy,
        yTy
    ]).reshape(-1, 1)
     
    x1 = np.linalg.solve(A, b)
    return x1[0],x1[1],x1[2],A


def SE_old_with_additive(Z,y, s2a_hat,s2gxg_hat,s2e_hat,T):
    n = Z.shape[0]
    m = Z.shape[1]
    ZZ = Z*Z
    p = m*(m-1)/2
    
    def k_middle(left,right):
        return float(np.sum((left@Z)* (right@Z))/m)
        
    def w_middle(left,right):
        return float((0.5*(left@(np.sum((Z @(Z.T@(right.reshape(-1, 1) * Z)))*Z,axis = 1)) - (left @ (ZZ @ (ZZ.T @ right))))/p))
        
    def no_middle(left,right):
        return float(left@right)

    Ky = Z @ (Z.T @ y) / m
    Wy = (0.5*((np.sum((Z @(Z.T@(y.reshape(-1, 1) * Z)))*Z,axis = 1)) - ((ZZ @ (ZZ.T @ y))))/p)
    side = [Ky,Wy,y]
    middle = [k_middle,w_middle,no_middle]
    var_component= [s2a_hat,s2gxg_hat,s2e_hat]
    cov_q = np.zeros((3,3))
    for i in range(3):
         for j in range(i, 3):
            val = 0
            for l in range(3):
                val += 2*(var_component[l].item() * middle[l](side[i], side[j]))
            cov_q[i, j] = val
            cov_q[j,i] = val
             
    X = np.linalg.solve(T, cov_q)
    Cov_theta = np.linalg.solve(T, X.T).T  
    return Cov_theta[0,0],Cov_theta[1,1],Cov_theta[2,2]


def SE_gxg(Z, y, s2gxg_hat, s2e_hat, T):
    n = Z.shape[0]
    m = Z.shape[1]
    ZZ = Z * Z
    p = m * (m - 1) / 2
    
    def w_middle(left, right):
        return float((0.5 * (left @ (np.sum((Z @ (Z.T @ (right.reshape(-1, 1) * Z))) * Z, axis=1)) 
                     - (left @ (ZZ @ (ZZ.T @ right)))) / p))
    
    def no_middle(left, right):
        return float(left @ right)
    
    Wy = (0.5 * ((np.sum((Z @ (Z.T @ (y.reshape(-1, 1) * Z))) * Z, axis=1)) 
          - ((ZZ @ (ZZ.T @ y)))) / p)
    
    side = [Wy, y]
    middle = [w_middle, no_middle]
    var_component = [s2gxg_hat, s2e_hat]
    
    cov_q = np.zeros((2, 2))
    for i in range(2):
        for j in range(i, 2):
            val = 0
            for l in range(2):
                val += 2 * (var_component[l].item() * middle[l](side[i], side[j]))
            cov_q[i, j] = val
            cov_q[j, i] = val
    
    X = np.linalg.solve(T, cov_q)
    Cov_theta = np.linalg.solve(T, X.T).T
    
    return Cov_theta[0, 0], Cov_theta[1, 1]


############################################# New method (std mom)

def compute_weight_matrices(Z):
    n, m = Z.shape
    
    ZTZ = Z.T @ Z
    E_ZaZb = ZTZ / n
    
    D = Z * Z
    DTD = D.T @ D
    E_ZaZb_sq = DTD / n
    
    V = E_ZaZb_sq - E_ZaZb ** 2
    
    S = 1.0 / V
    R = E_ZaZb / V
    T = E_ZaZb ** 2 / V
    
    np.fill_diagonal(S, 0)
    np.fill_diagonal(R, 0)
    np.fill_diagonal(T, 0)
    
    return S, R, T


def compute_WU(Z, U, S, R, T):

    n, m = Z.shape
    nmc = U.shape[1]
    p = m * (m - 1) // 2
    
    # Precompute shared terms
    ZR = Z @ R
    term2_base = 0.5 * np.sum(Z * ZR, axis=1, keepdims=True)
    sum_T = 0.5 * np.sum(T)
    sum_U = np.sum(U, axis=0, keepdims=True)
    
    # Term 2 & 4 
    term2 = term2_base @ sum_U
    term4 = sum_T * np.ones((n, 1)) @ sum_U
    
    # Term 1 & 3 
    term1 = np.zeros((n, nmc))
    term3 = np.zeros((n, nmc))
    
    for k in range(nmc):
        u = U[:, k]
        M = Z.T @ (u[:, np.newaxis] * Z)
        SM = S * M
        term1[:, k] = 0.5 * np.sum(Z * (Z @ SM), axis=1)
        term3[:, k] = 0.5 * np.sum(R * M)
    
    return (term1 - term2 - term3 + term4) / p


def MoM_only_M_std(Z, y, nmc=40):
    
    n = Z.shape[0]
    m = Z.shape[1]
    p = m * (m - 1) / 2 

    S, R, T = compute_weight_matrices(Z)
    
    # Random vectors for stochastic trace estimation
    u = np.random.choice([-1, 1], size=(n, nmc))
    
    # Compute Wu using standardized W
    Wu = compute_WU(Z, u, S, R, T)  # This should return W @ u
    
    # tr(W) estimation: tr(W) = mean(u.T @ W @ u)
    trW_app = np.mean(np.sum(u * Wu, axis=0))
    
    # tr(W^2) estimation: tr(W^2) = mean(||W @ u||^2)
    trW2_app = np.mean(np.sum(Wu * Wu, axis=0))
    
    # Quadratic forms
    yTy = y @ y
    
    # y.T @ W @ y - need to compute W @ y first
    y_col = y.reshape(-1, 1)
    Wy = compute_WU(Z, y_col, S, R, T).flatten()  # W @ y
    yTWy = y @ Wy
    
    # Solve moment equations
    # E[y.T @ W @ y] = s2gxg * tr(W^2) + s2e * tr(W)
    # E[y.T @ y] = s2gxg * tr(W) + s2e * n
    A = np.array([
        [trW2_app, trW_app],
        [trW_app, n]
    ])
    
    b = np.array([yTWy, yTy])
    
    # Solve A @ [s2gxg, s2e].T = b
    x = np.linalg.solve(A, b)

    return x[0], x[1], A


def MoM_exact(Z, y, Lgxg):

    n = len(y)
    
    W = Lgxg @ Lgxg.T
    
    q1 = y.T @ W @ y  # y'Wy
    q2 = y.T @ y      # y'y
        
    trW = np.trace(W)
    trWW = np.trace(W @ W)
       
    T = np.array([
        [trWW, trW],
        [trW,  n]
    ])
    
    q = np.array([q1, q2])
    
    theta_hat = np.linalg.solve(T, q)
    
    s2gxg_hat = theta_hat[0]
    s2e_hat = theta_hat[1]
    
    return s2gxg_hat, s2e_hat,T
