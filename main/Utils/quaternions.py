import numpy as np

# Quaternion product
def qprod(Aq_B, Bq_C):
    Aq_B = np.asarray(Aq_B)
    Bq_C = np.asarray(Bq_C)
    s_AB, x_AB, y_AB, z_AB = np.moveaxis(Aq_B, -1, 0)
    s_BC, x_BC, y_BC, z_BC = np.moveaxis(Bq_C, -1, 0)

    v_AB = np.stack([x_AB,y_AB,z_AB], axis=-1)
    v_BC = np.stack([x_BC,y_BC,z_BC], axis=-1)

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

# Express in frame A a vector initially expressed in frame B
def qvecprod(Aq_B, v_B):
    Aq_B  = np.asarray(Aq_B)
    v_B   = np.asarray(v_B)
    Aq_B  = np.broadcast_to(Aq_B, v_B.shape[:-1] + (4,))
    zeros = np.zeros(v_B.shape[:-1] + (1,))
    q_vec = np.concatenate((zeros, v_B), axis=-1)
    q_new = qprod(qprod(Aq_B, q_vec), qtrans(Aq_B))
    v_A   = q_new[..., 1:4]
    return v_A

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

# Convert a quaternion to a cosine matrix
def q2cosmat(Aq_B):
    Aq_B = np.asarray(Aq_B)

    w = Aq_B[..., 0]
    x = Aq_B[..., 1]
    y = Aq_B[..., 2]
    z = Aq_B[..., 3]

    C_AB = np.zeros(Aq_B.shape[:-1] + (3, 3))

    C_AB[..., 0, 0] = 2 * (w*w + x*x) - 1
    C_AB[..., 0, 1] = 2 * (x*y + w*z)
    C_AB[..., 0, 2] = 2 * (x*z - w*y)

    C_AB[..., 1, 0] = 2 * (x*y - w*z)
    C_AB[..., 1, 1] = 2 * (w*w + y*y) - 1
    C_AB[..., 1, 2] = 2 * (y*z + w*x)

    C_AB[..., 2, 0] = 2 * (x*z + w*y)
    C_AB[..., 2, 1] = 2 * (y*z - w*x)
    C_AB[..., 2, 2] = 2 * (w*w + z*z) - 1

    return C_AB

# Convert a cosine matrix to quaternion
def cosmat2q(C_AB):
    C_AB = np.asarray(C_AB)

    b0 = 0.25 * (1.0 + np.trace(C_AB))
    b1 = 0.25 * (1.0 + C_AB[0, 0] - C_AB[1, 1] - C_AB[2, 2])
    b2 = 0.25 * (1.0 + C_AB[1, 1] - C_AB[0, 0] - C_AB[2, 2])
    b3 = 0.25 * (1.0 + C_AB[2, 2] - C_AB[0, 0] - C_AB[1, 1])

    b = np.array([b0, b1, b2, b3])
    ind = np.argmax(b)

    qvec = np.zeros(3)
    qsca = 0.0

    if ind == 0:
        qsca = np.sqrt(b0)
        qvec[0] = 0.25 * (C_AB[1, 2] - C_AB[2, 1]) / qsca
        qvec[1] = 0.25 * (C_AB[2, 0] - C_AB[0, 2]) / qsca
        qvec[2] = 0.25 * (C_AB[0, 1] - C_AB[1, 0]) / qsca
    elif ind == 1:
        qvec[0] = np.sqrt(b1)
        qsca = 0.25 * (C_AB[1, 2] - C_AB[2, 1]) / qvec[0]
        qvec[1] = 0.25 * (C_AB[0, 1] + C_AB[1, 0]) / qvec[0]
        qvec[2] = 0.25 * (C_AB[2, 0] + C_AB[0, 2]) / qvec[0]
    elif ind == 2:
        qvec[1] = np.sqrt(b2)
        qsca = 0.25 * (C_AB[2, 0] - C_AB[0, 2]) / qvec[1]
        qvec[0] = 0.25 * (C_AB[0, 1] + C_AB[1, 0]) / qvec[1]
        qvec[2] = 0.25 * (C_AB[1, 2] + C_AB[2, 1]) / qvec[1]
    elif ind == 3:
        qvec[2] = np.sqrt(b3)
        qsca = 0.25 * (C_AB[0, 1] - C_AB[1, 0]) / qvec[2]
        qvec[0] = 0.25 * (C_AB[2, 0] + C_AB[0, 2]) / qvec[2]
        qvec[1] = 0.25 * (C_AB[1, 2] + C_AB[2, 1]) / qvec[2]

    # Assemble quaternion: [scalar, vector]
    Aq_B = np.hstack((qsca, qvec))

    # Normalize and enforce sign convention
    Aq_B = qnorm(Aq_B)
    Aq_B = qchksign(Aq_B)

    return Aq_B

# Return 3 axis defining the A frame relative to B frame
def q2axis(Aq_B):
    C_AB = q2cosmat(Aq_B)
    Ax_B = C_AB[..., 0, :]
    Ay_B = C_AB[..., 1, :]
    Az_B = C_AB[..., 2, :]
    return Ax_B, Ay_B, Az_B

# Convert roll, pitch and yaw to cosine matrix
def rpy2cosmat(rpy):
    rpy = np.asarray(rpy)

    # Compute cos and sin
    c1 = np.cos(rpy[0])
    s1 = np.sin(rpy[0])
    c2 = np.cos(rpy[1])
    s2 = np.sin(rpy[1])
    c3 = np.cos(rpy[2])
    s3 = np.sin(rpy[2])

    # Compute transformation matrix columns
    x = np.array([
        c2 * c3,
        c2 * s3,
        -s2
    ])
    y = np.array([
        -c1 * s3 + s1 * s2 * c3,
         c1 * c3 + s1 * s2 * s3,
         s1 * c2
    ])
    z = np.array([
         s1 * s3 + c1 * s2 * c3,
        -s1 * c3 + c1 * s2 * s3,
         c1 * c2
    ])

    # Direction cosine matrix
    C = np.vstack((x, y, z))

    return C

# Convert roll, pitch and yaw to quaternion
def rpy2q(rpy):
    C = rpy2cosmat(rpy)
    q = cosmat2q(C)
    return q
