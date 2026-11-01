# Sanity test of libdfm.so: 60 synthetic days with a 1.5 cm rain on day 30; prints 14:00 moisture (%) of 1-h and 10-h sticks.
import ctypes, numpy as np, sys
L = ctypes.CDLL(sys.argv[1]); f = np.ctypeslib.ndpointer(np.float32, flags='C_CONTIGUOUS'); L.run_dfm.argtypes = [ctypes.c_int, ctypes.c_int, f, f, f, f, f]
n = 24 * 60; h = np.arange(n) % 24
T = (20 + 8 * np.sin((h - 9) / 24 * 2 * np.pi)).astype(np.float32); RH = (0.5 - 0.25 * np.sin((h - 9) / 24 * 2 * np.pi)).astype(np.float32)
SW = np.maximum(0, 900 * np.sin((h - 6) / 12 * np.pi)).astype(np.float32); R = np.zeros(n, np.float32); R[24 * 30 + 14:24 * 30 + 20] = 0.25
for lag in (1, 10):
    o = np.zeros(n, np.float32); L.run_dfm(lag, n, T, RH, SW, R, o); print(lag, np.round(o[24 * 28 + 14::24][:12] * 100, 1))
