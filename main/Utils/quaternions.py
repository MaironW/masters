import numpy as np

# Quaternion product
def qprod(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    w1, x1, y1, z1 = np.moveaxis(a, -1, 0)
    w2, x2, y2, z2 = np.moveaxis(b, -1, 0)
    q = np.stack([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ], axis=-1)
    return q

# Quaternion transpose (conjugate)
def qtrans(q):
    q = np.asarray(q)
    q_conj = q.copy()
    q_conj[..., 1:] *= -1
    return q_conj

# Rotate vector by quaternion
def qvecrot(v, q):
    v = np.asarray(v)
    q = np.asanyarray(q)
    q = q / np.linalg.norm(q, axis=-1, keepdims=True)
    q0 = q[..., 0]
    qv = q[..., 1:]
    tmp = 2*np.cross(qv, v)
    return v + q0[..., None] * tmp + np.cross(qv, tmp)
