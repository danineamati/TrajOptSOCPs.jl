# Other Utils

"""
    safeNorm(arr, floor = 10^(-20))

Acts as a safe 2-norm that operates on an array that would approach zero.
Usually this is used with semilog or log-log plots of the residuals

See also: [`calcNormGradResiduals`](@ref)
"""
function safeNorm(arr, floor = 10^(-20))
    # Take the norm to get non-negative.
    # Take max with floor to get positive.
    return max.(norm.(arr), floor)
end
