# Histogram of Violation

using Plots

function plotDynConstVio(dynVio, constVio, useLog10 = true)

    if useLog10
        dV = log10.(dynVio)
        cV = log10.(constVio)
        header = "Log10 of "
    else
        dV = dynVio
        cV = constVio
        header = ""
    end

    plt = scatter(dV, cV)
    xlabel!(header * "Dynamics Constraint Violation")
    ylabel!(header * "Inequality Constraint Violation")
    title!("Clustering of Constraint Violation")

    display(plt)

    return plt
end
