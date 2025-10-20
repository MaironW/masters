import numpy as np

def qprod(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    q = np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])
    return q


def qtrans(q):
    q0, q1, q2, q3 = q
    q_conj = np.array([q0, -q1, -q2, -q3])
    return q_conj

def qvecrot(v, q):
    q = q / np.linalg.norm(q)
    vq = np.array([0.0, *v])
    v_new = qprod(qprod(q, vq), qtrans(q))
    return v_new[1:]
