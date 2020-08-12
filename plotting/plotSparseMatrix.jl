# Plot sparse matrices

using TrajOptSOCPs
using Plots, SparseArrays


function plotHessianAL(al::augLag, traj)
    h = evalHessAl(al, traj)
    plth = spy(sparse(h))
    return plth
end
