import numpy as np
import scipy.linalg

debuge = True 

def logm(A):
    
    if debuge==True :
        print("In the logm function")
        print(f'The dimension of dim is :',A.ndim)
        print(f'The type of A is :',A.dtype)
        
    if not np.issubdtype(A.dtype, float) or not A.ndim == 2:
        raise ValueError("Input must be a 2D array of floats.")

    if A.shape[0] != A.shape[1]:
        raise ValueError("Input matrix must be square.")

    maxroots = 100
    exitflag = 0
    
    if debuge==True :
        print("In the logm function 2")
    if not np.all(np.isfinite(A)):
        L = np.full_like(A, np.nan, dtype=A.dtype)
        return L, exitflag
    
    if debuge==True :
        print("In the logm function 2.2")

    # Check for triangularity.
    schur_input, T = scipy.linalg.schur(A)
    if debuge==True :
        print("In the logm function 2.3")
        print(f'schur_input is : \n{schur_input}')
        print(f'T is : \n {T}')
    print("In the logm function 2.3..0")    
    if not schur_input:
        print("In the logm function 2.3..1")
        Q, T = scipy.linalg.schur(A, output='full')
        if debuge==True :
            print("In the logm function 2.3..2")
            print(f'Q is : \n{Q}')
            print(f'T is : \n {T}')
    else:
        print("In the logm function 2.3..1")
        Q = np.eye(A.shape[0], dtype=A.dtype)
        if debuge==True :
            print("In the logm function 2.3..2")
            print(f'Q is : \n{Q}')
            
        
    if debuge==True :
        print("In the logm function 2.4")
    stay_real = np.isreal(A).all()

    if debuge==True :
        print("In the logm function 3")
        
    # Compute the logarithm.
    if np.all(np.diag(T)):
        d = np.diag(T)
        if np.any((np.real(d) <= 0) & (np.imag(d) == 0)):
            print("Warning: Non-positive real eigenvalues.")
        if schur_input:
            L = np.diag(np.log(d))
        else:
            logd = np.log(d)
            L = np.dot(Q * logd, Q.T)
            if np.isreal(logd).all():
                L = (L + L.T) / 2
    else:
        n = T.shape[0]
        ei = np.linalg.eigvals(T)
        warns = np.any(ei == 0)
        if np.any((np.real(ei) < 0) & (np.imag(ei) == 0)):
            warns = True
            if stay_real:
                if schur_input:
                    Q = np.eye(n, dtype=A.dtype)
                    schur_input = False  # Need to undo rsf2csf at the end.
                Q, T = scipy.linalg.rsf2csf(Q, T)
        if warns:
            print("Warning: Non-positive real eigenvalues.")
        
        if debuge==True :
            print("In the logm function 4")
        # Get block structure of Schur factor.
        blockformat = qtri_struct(T)

        # Get parameters.
        s, m, Troot, exitflag = logm_params(T, maxroots)

        # Compute Troot - I = T(1/2^s) - I more accurately.
        Troot = recompute_diag_blocks_sqrt(Troot, T, blockformat, s)

        # Compute Pade approximant.
        L = pade_approx(Troot, m)

        # Scale back up.
        L = 2**s * L
        
        if debuge==True :
            print("In the logm function 5")
        # Recompute diagonal blocks.
        L = recompute_diag_blocks_log(L, T, blockformat)

        # Combine if needed
        if not schur_input:
            L = np.dot(np.dot(Q, L), Q.T)

    return L, exitflag

def logm_params(T, maxroots):
    if debuge==True :
        print("In the logm_params function")
    exitflag = 0
    n = T.shape[0]
    I = np.eye(n, dtype=T.dtype)

    xvals = np.array([
        1.586970738772063e-005,
        2.313807884242979e-003,
        1.938179313533253e-002,
        6.209171588994762e-002,
        1.276404810806775e-001,
        2.060962623452836e-001,
        2.879093714241194e-001
    ])

    mmax = 7
    foundm = False

    # Get initial s0 so that T^(1/2^s0) < xvals(mmax).
    s = 0
    d = np.linalg.eigvals(T)
    while np.linalg.norm(d - 1, np.inf) > xvals[mmax - 1] and s < maxroots:
        d = np.sqrt(d)
        s += 1

    s0 = s
    if s == maxroots:
        print("Warning: Too many matrix square roots.")
        exitflag = 1

    Troot = T.copy()
    for k in range(1, min(s, maxroots) + 1):
        Troot = sqrtm_tri(Troot)

    # Compute value of s and m needed.
    TrootmI = Troot - I
    d2 = np.linalg.norm(normAm(TrootmI, 2))**(1/2)
    d3 = np.linalg.norm(normAm(TrootmI, 3))**(1/3)
    a2 = max(d2, d3)

    if a2 <= xvals[1]:
        m = np.where(a2 <= xvals[:2])[0][0] + 1
        foundm = True

    p = 0
    while not foundm:
        more = False  # More norm checks needed.
        if s > s0:
            d3 = np.linalg.norm(normAm(TrootmI, 3))**(1/3)

        d4 = np.linalg.norm(normAm(TrootmI, 4))**(1/4)
        a3 = max(d3, d4)

        if a3 <= xvals[mmax - 1]:
            j = np.where(a3 <= xvals[2:mmax])[0][0] + 2
            if j <= 6:
                m = j
                break
            else:
                if a3/2 <= xvals[4] and p < 2:
                    more = True
                    p += 1

        if not more:
            d5 = np.linalg.norm(normAm(TrootmI, 5))**(1/5)
            a4 = max(d4, d5)
            eta = min(a3, a4)

            if eta <= xvals[mmax - 1]:
                m = np.where(eta <= xvals[5:mmax])[0][0] + 5
                break

        if s == maxroots:
            if exitflag == 0:
                print("Warning: Too many matrix square roots.")
            exitflag = 1
            m = mmax  # No good value found so take the largest.
            break

        Troot = sqrtm_tri(Troot)
        TrootmI = Troot - I
        s += 1

    return s, m, Troot, exitflag

def recompute_diag_blocks_log(L, T, blockStruct):
    if debuge==True :
        print("In the recompute_diag_blocks_log function")
    n = len(T)
    last_block = 0

    for j in range(n - 1):
        if blockStruct[j] == 0:
            if last_block != 0:
                last_block = 0
                continue
            else:
                last_block = 0
                L[j, j] = np.log(T[j, j])
        elif blockStruct[j] == 1:
            last_block = 1
            a1 = T[j, j]
            a2 = T[j + 1, j + 1]
            loga1 = np.log(a1)
            loga2 = np.log(a2)
            L[j, j] = loga1
            L[j + 1, j + 1] = loga2

            if (a1 < 0 and np.imag(a1) == 0) or (a2 < 0 and np.imag(a1) == 0):
                continue

            if a1 == a2:
                a12 = T[j, j + 1] / a1
            else:
                z = (a2 - a1) / (a2 + a1)
                if check_condition(z):
                    a12 = T[j, j + 1] * (loga2 - loga1) / (a2 - a1)
                else:
                    dd = (2 * np.arctanh(z) + 2j * np.pi * unwinding(loga2 - loga1)) / (a2 - a1)
                    a12 = T[j, j + 1] * dd

            L[j, j + 1] = a12
        elif blockStruct[j] == 2:
            last_block = 2
            f = 0.5 * np.log(T[j, j]**2 - T[j, j + 1] * T[j + 1, j])
            t = np.arctan2(np.sqrt(-T[j, j + 1] * T[j + 1, j]), T[j, j]) / np.sqrt(-T[j, j + 1] * T[j + 1, j])
            L[j, j] = f
            L[j + 1, j] = t * T[j + 1, j]
            L[j, j + 1] = t * T[j, j + 1]
            L[j + 1, j + 1] = f

    if blockStruct[-1] == 0:
        L[n - 1, n - 1] = np.log(T[n - 1, n - 1])

    return L

def sqrt_obo(a, s):
    if debuge==True :
        print("In the sqrt_obo function")
    if s == 0:
        val = a - 1
    else:
        n0 = s
        if np.angle(a) >= np.pi / 2:
            a = np.sqrt(a)
            n0 = s - 1
        z0 = a - 1
        a = np.sqrt(a)
        r = 1 + a
        for i in range(1, n0):
            a = np.sqrt(a)
            r = r * (1 + a)
        val = z0 / r

    return val

def recompute_diag_blocks_sqrt(Troot, T, blockStruct, s):
    if debuge==True :
        print("In the recompute_diag_blocks_sqrt function")
    n = len(T)
    last_block = 0

    for j in range(n - 1):
        if blockStruct[j] == 0:
            if last_block != 0:
                last_block = 0
                continue
            else:
                last_block = 0
                Troot[j, j] = sqrt_obo(T[j, j], s)
        else:
            last_block = blockStruct[j]
            I = np.eye(2, dtype=T.dtype)
            if s == 0:
                Troot[j:j + 2, j:j + 2] = T[j:j + 2, j:j + 2] - I
                continue
            A = sqrtm_tbt(T[j:j + 2, j:j + 2])
            Z0 = A - I
            if s == 1:
                Troot[j:j + 2, j:j + 2] = Z0
                continue
            A = sqrtm_tbt(A)
            P = A + I
            for i in range(s - 2):
                A = sqrtm_tbt(A)
                P = P * (I + A)

            Troot[j:j + 2, j:j + 2] = np.linalg.solve(P, Z0)

            if T[j + 1, j] == 0 and T[j, j] >= 0 and T[j + 1, j + 1] >= 0:
                Troot[j, j + 1] = powerm2by2(T[j:j + 2, j:j + 2], 1. / (2.**s))

    if blockStruct[-1] == 0:
        Troot[n - 1, n - 1] = sqrt_obo(T[n - 1, n - 1], s)

    return Troot

def powerm2by2(A, p):
    if debuge==True :
        print("In the powerm2by2 function")
    a1 = A[0, 0]
    a2 = A[1, 1]

    if a1 == a2:
        x12 = p * A[0, 1] * a1**(p - 1)
    else:
        z = (a2 - a1) / (a2 + a1)
        if check_condition(z):
            x12 = A[0, 1] * (a2**p - a1**p) / (a2 - a1)
        else:
            loga1 = np.log(a1)
            loga2 = np.log(a2)
            w = np.arctanh(z) + 1j * np.pi * unwinding(loga2 - loga1)
            dd = 2 * np.exp(p * (loga1 + loga2) / 2) * np.sinh(p * w) / (a2 - a1)
            x12 = A[0, 1] * dd

    return x12

def gauss_legendre(n):
    if debuge==True :
        print("In the gauss_legendre function")
    k = np.arange(1, n)
    v = k / np.sqrt((2 * k)**2 - 1)
    V, x = np.linalg.eig(np.diag(v, -1) + np.diag(v, 1))
    w = 2 * V[0]**2

    return x, w

def pade_approx(T, m):
    if debuge==True :
        print("In the pade_approx function")
    nodes, wts = gauss_legendre(m)
    # Convert from [-1,1] to [0,1].
    nodes = (nodes + 1) / 2
    wts = wts / 2
    n = T.shape[0]
    L = np.zeros_like(T, dtype=T.dtype)

    for j in range(m):
        K = nodes[j] * T
        K[np.diag_indices(n)] += 1
        L = L + wts[j] * np.linalg.solve(K, T)

    return L
