from sage.crypto.lwe import LWE
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler

def lwegen(n,m,q,sig):
    lwe = LWE(n=n, q=q, D=DiscreteGaussianDistributionIntegerSampler(sig))
    samps = [lwe() for _ in range(m)]

    A = matrix([ai for (ai,bi) in samps])
    b = vector(ZZ, [bi for (ai,bi) in samps])

    return A, b, lwe._LWE__s

def solvelwe(A,b,sig,beta):
    q = A.base_ring().order()
    m = A.nrows()

    Aext  = (A.T).change_ring(ZZ).stack(q*identity_matrix(m))
    basis = Aext.hermite_form()[:m,:]

    L = block_matrix([[matrix(ZZ,1,1,[sig]), matrix(b)], \
                    [matrix(ZZ,m,1,[0]*m), basis    ]])

    Lred = L.BKZ(block_size=beta)

    ecand = Lred[0][1:]

    try:
        scand = A \ (b - ecand)
        return scand
    except:
        return None
    

def testsolvelwe(n=20,m=40,q=1001,sig=3.0,beta=25):
    A,b,s = lwegen(n,m,q,sig)
    scand = solvelwe(A,b,sig,beta)

    if scand == s:
        print "Success!"
    else:
        print "Failure."

# vim: ft=python
