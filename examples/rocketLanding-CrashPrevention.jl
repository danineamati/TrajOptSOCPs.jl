#=

In this file we simulate a simple rocket landing

Units in kg, m, s

=#
include("..\\src\\rocket\\rocket-setup.jl")
include("..\\src\\rocket\\ground.jl")

include("..\\src\\dynamics\\trajectory-setup.jl")
include("..\\src\\objective\\LQR_objective.jl")
include("..\\src\\constraints\\constraintManager.jl")
include("..\\src\\constraints\\constraintManagerDynamics.jl")
include("..\\src\\auglag\\auglag-core.jl")

include("..\\src\\solver\\AL-Primal-Main-Solver.jl")
include("..\\src\\solver\\QP-AffineEquality-Solver.jl")
include("..\\src\\other_utils\\parsePrimalDual.jl")

include("..\\src\\results\\trajectoryParsing.jl")
include("..\\src\\results\\plotTrajectory.jl")
include("..\\src\\results\\plotConstraintViolation.jl")
include("..\\src\\results\\plotObjective.jl")
include("..\\src\\results\\batchPlots.jl")


# Based on the Falcon 9
# 549,054 kg (Mass)
# 282 s (Specific Impulse)
# select "y" as the vertical
mass = 1
isp = 1
grav = [0; -9.81]
deltaTime = 1
rocket = rocket_simple(mass, isp, grav, deltaTime)

# in m
# The Karman Line (100 km)
rocketStart = [2.0; 20.0; 0.0; -45.0]
rocketEnd = [0.0; 0.0; 0.0; 0.0]#[-5.0; 0.0; 0.0; 0.0]

uHover = mass * grav

# Number of time steps to discretize the trajectory
NSteps = 60
# Initialize the trajectory with a line
initTraj = initializeTraj(rocketStart, rocketEnd, uHover, uHover, NSteps)

# Use a Linear Quadratic Regulator as the cost function
lqrQMat = 0.0001 * Diagonal(I, size(rocketStart, 1))
lqrRMat = 0.0025 * Diagonal(I, Int64(size(rocketStart, 1) / 2))
costFun = makeLQR_TrajReferenced(lqrQMat, lqrRMat, NSteps, initTraj)

# Create the Dynamics Constraints
ADyn, BDyn = rocketDynamicsFull(rocket, rocketStart, rocketEnd, NSteps)
dynConstraint = AL_AffineEquality(ADyn, BDyn)
lambdaInit = -1 * ones(size(BDyn))

# Create the Ground Constraints
groundConstraint = makeGroundConstraint(NSteps, size(grav, 1), size(grav, 1))
groundLambda = zeros(size(groundConstraint.b, 1))

# Create the constraint manager
# cMRocket = constraintManager_Base([dynConstraint], [lambdaInit])
# cMRocket = constraintManager_Dynamics([], [], dynConstraint, lambdaInit)
cMRocket = constraintManager_Dynamics([groundConstraint], [groundLambda],
                                        dynConstraint, lambdaInit)

# Initialize the primal-dual vector
initTrajPD = [initTraj; lambdaInit]

# Initialize the Augmented Lagrangian penalty on the constraints
penaltyStart = 1.0

# Test that the evaluations work
println("\n--------------------------------------------")
println(" Testing constraints are inputted correctly ")
println("--------------------------------------------")
println("Starting constraint violation: ")
println([evalConstraints(cMRocket, initTraj, penaltyStart)])
println("Starting gradient of constraints: ")
println(evalGradConstraints(cMRocket, initTraj, penaltyStart))
println("Starting hessian of constraints: ")
println(evalHessConstraints(cMRocket, initTraj))


# Equiped with the constraint term and the objective term, I now build the
# Augmented Lagrangian
alRocket = augLag(costFun, cMRocket, penaltyStart)

# Test that the evaluations work
println("--------------------------------------------")
println("           Testing AL is Functional         ")
println("--------------------------------------------")
println("Evaluating augmented lagrangian: ")
println([evalAL(alRocket, initTrajPD)])
println("Evaluating gradient of augmented lagrangian: ")
println(evalGradAL(alRocket, initTrajPD))
println("Evaluating hessian of augmented lagrangian: ")
println(size(evalHessAl(alRocket, initTrajPD)))

# Next we select resonable solver parameters
currSolveParams = solverParams(0.1, 0.5,
                                18, 4,
                                10^-4,
                                10, 10^6,
                                2.5, 2, 0.2, 0.2, 0.4)
solParamPrint(currSolveParams)
println()
println()

println("--------------------------------------------")
println("             Beginning Solve                ")
println("--------------------------------------------")

# Solve the Trajectory Optimization problem
useALMethod = true

if useALMethod
    # Use an augmented lagrangian method
    trajLambdaSolved, resArr = ALPrimalNewtonMain(initTrajPD, alRocket,
                                            currSolveParams)
    trajStatesAllPD = [parsePrimalDualVec(trajL, size(initTraj, 1))
                                                for trajL in trajLambdaSolved]
    trajStatesAll = [trajPDThis.primals for trajPDThis in trajStatesAllPD]
else
    # Solve only the objective and linear dynamics
    costQ = costFun.lqr_qp.QR_Full
    costP = - (initTraj' * costQ)'
    trajLambdaSolvedTrue = solveQP_AffineEq(costQ, costP, ADyn, BDyn)
    trajStateLastTrue = parsePrimalDualVec(trajLambdaSolvedTrue,
                                                        size(initTraj, 1))
end

# Blocked so that it can be run independently after the fact
if true
    # Get the parsed list of trajectories
    nDim = size(grav, 1)
    pltTraj, pltCV, pltCV2, pltObj, plts, pltv, pltu =
                                    batchPlot(trajStatesAll, cMRocket, nDim)
end

# Blocked so that it can be run independently after the fact
if true
    header = "freefallingAllPlots18_4" * string(Int64(rocketStart[4])) * "_"
    saveBulk(pltTraj, pltCV, pltCV2, pltObj, plts, pltv, pltu, header)
end
