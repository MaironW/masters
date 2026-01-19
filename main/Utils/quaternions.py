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

def qprod_ref(Aq_B, Bq_C):
    Aq_B = np.asarray(Aq_B)
    Bq_C = np.asarray(Bq_C)
    w1, x1, y1, z1 = np.moveaxis(Aq_B, -1, 0)
    w2, x2, y2, z2 = np.moveaxis(Bq_C, -1, 0)

    s_AB = w1
    s_BC = w2
    v_AB = np.stack([x1,y1,z1], axis=-1)
    v_BC = np.stack([x2,y2,z2], axis=-1)

    s_AC = s_AB * s_BC - np.sum(v_AB * v_BC, axis=-1)
    v_AC = s_AB[..., None] * v_BC + s_BC[..., None] * v_AB - np.cross(v_AB, v_BC)

    Aq_C = np.concatenate([s_AC[..., None], v_AC], axis=-1)
    return Aq_C

# Difference between two quaternions
def qerr(Aq_B, Cq_B):
    Aq_B = np.asarray(Aq_B)
    Cq_B = np.asarray(Cq_B)
    Bq_A = qtrans(Aq_B)
    Cq_A = qchksign(qprod(Cq_B, Bq_A))
    return Cq_A

# Quaternion transpose (conjugate)
def qtrans(q):
    q = np.asarray(q)
    q_conj = q.copy()
    q_conj[..., 1:] *= -1
    return q_conj

# Rotate vector by quaternion
def qvecprod(Bq_A, v_A):
    v_A = np.asarray(v_A)
    Bq_A = np.asanyarray(Bq_A)
    vq = np.zeros((*v_A.shape[:-1], 4))
    vq[..., 1:] = v_A
    Aq_B = qtrans(Bq_A)
    vq_B = qprod(Aq_B, qprod(vq, Bq_A))
    v_B = vq_B[..., 1:]
    return v_B

# Convert small angle vector to quaternion
def rotvec2q(v):
    theta = np.linalg.norm(v)
    if theta < 1e-12:
        return np.array([1, 0, 0, 0])
    axis = v/theta
    half = 0.5 * theta
    q = np.array([
        np.cos(half),
        axis[0] * np.sin(half),
        axis[1] * np.sin(half),
        axis[2] * np.sin(half),
    ])
    q = qnorm(q)
    q = qchksign(q)
    return q

# Returns the quaternion which gives the minimum rotation from vector a to vector b
def vecqvec(a, b):
    a = qnorm(a)
    b = qnorm(b)

    v = np.cross(a, b)
    c = np.dot(a, b)

    if c < -0.999999:
        raise ValueError("180 deg rotation: choose secondary constraint")

    q = np.zeros(4)
    q[0] = 1.0 + c
    q[1:] = v
    q = qnorm(q)
    q = qchksign(q)
    return q

# Normalize a quaternion
def qnorm(q):
    q = q / np.linalg.norm(q, axis=-1, keepdims=True)
    return q

# Make a quaternion aways scalar positive
def qchksign(q):
    q = np.asarray(q)
    sign = np.where(q[..., 0] < 0, -1, 1)
    return q*sign[..., None]

# Return 3 axis vectors of the frame B relative to frame A
def q2axis(Bq_A):
    Bq_A = np.asarray(Bq_A)

    xA = np.array([1,0,0])
    yA = np.array([0,1,0])
    zA = np.array([0,0,1])

    Aq_B = qtrans(Bq_A)

    xB = qvecprod(Aq_B, xA)
    yB = qvecprod(Aq_B, yA)
    zB = qvecprod(Aq_B, zA)

    return xB, yB, zB
