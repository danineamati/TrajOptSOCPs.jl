# Parse Primals and Duals

struct primal_dual
    primals
    duals
end

function parsePrimalDualVec(pdVec, primalSize::Int64)
    return primal_dual(pdVec[1:primalSize], pdVec[primalSize + 1:end])
end
