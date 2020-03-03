import numpy as np


## (px, py, pz, E)
test_jets = np.array([
        [109.975, 2.8439, -34.4725, 117.936],
        [19.4617, -7.17918, -13.2483, 24.8724],
        [-135.97, 18.7476, -175.697, 223.586],
        [111.105, 2.87311, -34.8266, 119.147],
        [19.8691, -7.32946, -13.5256, 25.3931],
        [-137.405, 18.9455, -177.552, 225.946],
        [108.846, 2.81469, -34.1185, 116.724],
        [19.0543, -7.0289, -12.971, 24.3517],
        [-134.535, 18.5497, -173.842, 221.226],
        [110.136, 2.84805, -34.5228, 118.107],
        [19.3452, -7.13619, -13.169, 24.7235],
        [-136.194, 18.7785, -175.986, 223.955],
        [109.815, 2.83976, -34.4223, 117.764],
        [19.5783, -7.22217, -13.3276, 25.0213],
        [-135.746, 18.7167, -175.407, 223.218],
])

## lep (px, py, pz, E)
test_leps = np.array([
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
        [20.6456, 17.215, -1.57129, 26.927],
])

## met (px, py, pz, E)
test_met = np.array([
        [-61.1258, 6.08048, 0, 61.4275],
        [-61.1258, 6.08048, 0, 61.4275],
        [-61.1258, 6.08048, 0, 61.4275],
        [-58.4479, 5.60997, 0, 58.7165],
        [-58.4479, 5.60997, 0, 58.7165],
        [-58.4479, 5.60997, 0, 58.7165],
        [-63.8038, 6.55098, 0, 64.1392],
        [-63.8038, 6.55098, 0, 64.1392],
        [-63.8038, 6.55098, 0, 64.1392],
        [-60.5248, 5.9798, 0, 60.8195],
        [-60.5248, 5.9798, 0, 60.8195],
        [-60.5248, 5.9798, 0, 60.8195],
        [-61.7269, 6.18116, 0, 62.0356],
        [-61.7269, 6.18116, 0, 62.0356],
        [-61.7269, 6.18116, 0, 62.0356],
])

## nu (px, py, pz, E, nuChi2)
test_nu = np.array([
        [-24.2269, -24.5548, -83.8295, 90.6491, 2300.05],
        [69.4059, 264.345, 187.287, 331.318, 83739],
        [-68.4777, 3.72165, -18.2811, 70.9736, 59.6144],
        [-22.3191, -25.0685, -85.5329, 91.8828, 2246.46],
        [65.2455, 254.588, 185.985, 321.966, 77290],
        [-67.9439, 2.44385, -19.0715, 70.6121, 100.198],
        [-26.1209, -24.0112, -82.1058, 89.4439, 2354.05],
        [73.9188, 274.51, 188.667, 341.196, 90769.6],
        [-68.9743, 4.951, -17.5083, 71.3337, 29.2937],
        [-23.9098, -24.5728, -84.188, 90.9017, 2274.12],
        [71.1807, 266.946, 188.206, 334.288, 85449.8],
        [-68.3698, 3.44135, -18.3904, 70.8835, 67.9869],
        [-24.5427, -24.5356, -83.4717, 90.3984, 2326.18],
        [67.6535, 261.779, 186.368, 328.388, 82069.6],
        [-68.584, 3.99945, -18.1728, 71.0634, 51.7799],
])