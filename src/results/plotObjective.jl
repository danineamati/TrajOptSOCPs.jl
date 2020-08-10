# Plot the Objective over the course of the solve

include("..\\objective\\QP_Linear_objectives.jl")


function plotObjective(obj::objectiveFunc, trajList)
    fList = [fObjQP(obj, traj)[1] for traj in trajList]

    pltObj = plot(fList, markershape = :square)
    title!("Cost Function")
    ylabel!("Cost Function")
    xlabel!("Newton Step")

    return pltObj
end
