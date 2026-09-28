import sys; sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
from methods import *; from b1up import *; from ctf import *
telescope_env(); H = Harness(); E = np.logspace(-1, 2, 501)
p = DCP_Parameters(boost=-1.0)
H.call(E, p); s = H.state()
ref = np.load("exp2.npz")["b1up256x8_U8x8"]
for cfg in [(64,4,64,4,1),(64,4,64,8,1),(64,4,128,8,1),(128,5,128,8,1)]:
    t=time.time(); K, nr, nd = ctf_kernel(H, s, *cfg); k0,k1=K.arrays(); out=H.call(E,p,kernel=(k0,k1))
    print(">> ctf", cfg, "rays", nr, "deposits", nd, "max %.2e rms %.2e"%(np.max(np.abs(out/ref-1)), np.sqrt(np.mean((out/ref-1)**2))), "%.1fs"%(time.time()-t))
