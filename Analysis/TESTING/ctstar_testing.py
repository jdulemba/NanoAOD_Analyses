from pdb import set_trace
import uproot_methods
import numpy as np
import Utilities.make_variables as mv


## (px, py, pz, E)
thad_vecs = np.array([
    [-129.437, -7.20629, -393.285, 419.209],
    [-131.242, -7.27465, -398.786, 425.07],
    [9.14759, 10.0404, -271.21, 309.339],
    [-1.06654, -49.0759, -270.41, 312.948],
    [11.2448, 9.414, -270.852, 309.472],
    [-117.639, -1.00322, -278.847, 339.129],
    [-118.723, -0.910418, -283.005, 344.077],
    [-116.555, -1.09601, -274.688, 334.181],
    [-117.54, -1.42974, -277.23, 337.407],
    [-117.739, -0.576686, -280.464, 340.851],
])

tlep_vecs = np.array([
    [82.3734, -25.6928, -39.7761, 196.936],
    [81.6899, -26.6738, -42.2783, 197.303],
    [-82.952, -52.38, -221.113, 297.106],
    [-74.5939, 7.59137, -239.597, 304.606],
    [-85.3283, -53.4918, -230.396, 304.927],
    [118.805, 13.6745, -194.384, 286.082],
    [118.586, 15.3518, -199.306, 289.444],
    [119.163, 11.98, -189.509, 282.866],
    [118.852, 14.2936, -195.328, 286.774],
    [118.764, 13.0471, -193.447, 285.4],
])

thad_cth = np.array([
    0.862137,
    0.862534,
    0.29129,
    0.314984,
    0.238442,
    0.269041,
    0.260589,
    0.27686,
    0.264161,
    0.273875
])

def compare_outputs(expected, calculated):
    for val in range(len(calculated)):
        if calculated[val].round(4) != expected[val].round(4):
            print('Values not the same. Calculated=', calculated[val], ', expected=', expected[val])
#set_trace()

tops = uproot_methods.TLorentzVectorArray.from_cartesian(thad_vecs[:,0], thad_vecs[:,1], thad_vecs[:,2], thad_vecs[:,3])
tbars = uproot_methods.TLorentzVectorArray.from_cartesian(tlep_vecs[:,0], tlep_vecs[:,1], tlep_vecs[:,2], tlep_vecs[:,3])

top_ctstar, tbar_ctstar = mv.ctstar(tops, tbars, debug=True)

print('Calculated topctstar: ', top_ctstar)
print('Expected topctstar: ', thad_cth)

#set_trace()
compare_outputs(thad_cth, top_ctstar)

