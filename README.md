# HEFDiM-preview-version
Highly Efficient Feasible Direction Method (HEFDiM) for Structural Topology Optimization (Preview Version)

# HEFDiM
Highly Efficient Feasible Direction Method (HEFDiM) for Structural Topology Optimization

Feasible Direction Method (FDM) is a concise yet rigorous mathematical method for structural topology optimization, which can be easily applied to different types of problems with less modification. In addition, FDM always converges to a near optimum rapidly. We advance the state-of-the-art by proposing a highly efficient feasible direction method (HEFDiM), which substantially improves the efficiency of the FDM with negligible loss of accuracy. The proposed method can benefit us in at least four aspects: 
- analytical gradient projection 
- few heuristics and clear physical meaning 
- negligible memory and time-cost for updating 
- directly applied to different problems without extra efforts. 

The MBB problem [(CODE)](./MBB_300_100_OUR.m):

![MBB](./MBB_300_100_OUR.gif)

Force inverter mechanism [(CODE)](./MSE_SE_KO.m):

![Force inverter](./MSE_SE_KO.gif)

3D cantilever beam [(CODE)](./TO_3D_OUR.m):

![3D cantilever beam](./TO_3D_60_20_10_OUR.gif)
