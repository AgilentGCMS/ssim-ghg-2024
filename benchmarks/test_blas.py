import numpy as np
import time, tqdm, math

def str_with_err(value, error):
    digits = -int(math.floor(math.log10(error)))
    return "{0:.{2}f} ± {1:.{2}f}".format(value, error, digits)
    # return "{0:.{2}f}({1:.0f})".format(value, error*10**digits, digits)

n = 3000
p = 1000
q = 5000

num_repeats = 20
matmul_times = np.zeros(num_repeats, dtype=np.float64)
svd_times = np.zeros(num_repeats, dtype=np.float64)
create_times = np.zeros(num_repeats, dtype=np.float64)
cholesky_times = np.zeros(num_repeats, dtype=np.float64)

for rep in tqdm.tqdm(range(num_repeats), desc='Benchmarking matrix operations'):
    t1 = time.perf_counter()
    A = np.random.random((n,p))
    B = np.random.random((p,q))
    H = np.random.random((n,n))
    H = 0.5*(H+H.T)
    t2 = time.perf_counter()
    create_times[rep] = t2-t1

    t1 = time.perf_counter()
    C = np.matmul(A,B)
    t2 = time.perf_counter()
    matmul_times[rep] = t2-t1

    t1 = time.perf_counter()
    D = np.linalg.svd(C)
    t2 = time.perf_counter()
    svd_times[rep] = t2-t1

    t1 = time.perf_counter()
    D = np.linalg.cholesky(H)
    t2 = time.perf_counter()
    cholesky_times[rep] = t2-t1

print("Time for creating matrices = %s s"%str_with_err(create_times.mean(), create_times.std()))
print("Time for multiplying matrices = %s s"%str_with_err(matmul_times.mean(), matmul_times.std()))
print("Time for SVD of %i×%i matrix = %s s"%(C.shape[0], C.shape[1], str_with_err(svd_times.mean(), svd_times.std())))
print("Time for Cholesky decomposition of %i×%i matrix = %s s"%(H.shape[0], H.shape[1], str_with_err(cholesky_times.mean(), cholesky_times.std())))
