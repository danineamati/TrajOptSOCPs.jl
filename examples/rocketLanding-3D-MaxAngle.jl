#=

In this file we simulate a simple rocket landing

Units in kg, m, s

=#

# RUN:
# ] activate TrajOptSOCPs

using TrajOptSOCPs
using LinearAlgebra, SparseArrays

# If you don't want plotting comment the include and pyplot statements.
# Then, sent runplots = false
# include("../plotting/plotTrajectory.jl")
# include("../plotting/plotConstraintViolation.jl")
# include("../plotting/plotObjective.jl")
include("../plotting/batchPlots.jl")
pyplot()
runplots = true
saveplots = true


println("\n--------------------------------------------")
println("          Setting Up Rocket Landing ")
println("--------------------------------------------")

# Based on the Falcon 9
# 549,054 kg (Mass)
# 282 s (Specific Impulse)
# select "y" as the vertical
mass = 1
isp = 1
grav = [0; 0; -9.81]
deltaTime = 1
rocket = rocket_simple(mass, isp, grav, deltaTime)

# in m
# The Karman Line (100 km)
const rocketStart = [2.0; 5.0; 20.0; 4.0; -1.0; -15.0]
const rocketEnd = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0]#[-5.0; 0.0; 0.0; 0.0]

uHover = mass * grav

# Number of time steps to discretize the trajectory
const NSteps = 60
# Initialize the trajectory with a line
initTraj = initializeTraj(rocketStart, rocketEnd, uHover, uHover, NSteps)

# Use a Linear Quadratic Regulator as the cost function
const lqrQMat = 0.0 * Diagonal(I, size(rocketStart, 1))
const lqrRMat = 0.25 * Diagonal(I, Int64(size(rocketStart, 1) / 2))
costFun = makeLQR_TrajReferenced(lqrQMat, lqrRMat, NSteps, initTraj)

# Create the Dynamics Constraints
const ADyn, BDyn = rocketDynamicsFull(rocket, rocketStart, rocketEnd, NSteps)
dynConstraint = AL_AffineEquality(ADyn, BDyn)
lambdaInit = -1 * ones(size(BDyn))

# Create the Ground Constraints
groundConstraint = makeGroundConstraint(NSteps, size(grav, 1), size(grav, 1))
groundLambda = zeros(size(groundConstraint.b, 1))

# Create the Max Thrust Constraints
thrustMax = 20.0
maxThrustConstraint = makeMaxThrustConstraint(NSteps, size(grav, 1), thrustMax)
maxThrustLambda = zeros(size(maxThrustConstraint.indicatorList))

# Create the Max Thrust Angle Constraints
thrustAngleMax = 10.0 # Deg
maxAngleConstraint = makeMaxAngleConstraint(NSteps, size(grav, 1),
                                                thrustAngleMax, true)
maxAngleLambda = zeros(size(maxAngleConstraint.indicatorList))

# Create the constraint manager
cMRocket = constraintManager_Dynamics(
            [groundConstraint, maxThrustConstraint, maxAngleConstraint],
            [groundLambda, maxThrustLambda, maxAngleLambda],
            dynConstraint, lambdaInit
            )
# cMRocket = constraintManager_Dynamics(
#             [groundConstraint, maxThrustConstraint],
#             [groundLambda, maxThrustLambda],
#             dynConstraint, lambdaInit
#             )

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
println(size(evalGradConstraints(cMRocket, initTraj, penaltyStart)))
println("Starting hessian of constraints: ")
println(size(evalHessConstraints(cMRocket, initTraj)))


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
println(size(evalGradAL(alRocket, initTrajPD)))
println("Evaluating hessian of augmented lagrangian: ")
println(size(evalHessAl(alRocket, initTrajPD)))

# Next we select resonable solver parameters
currSolveParams = solverParams(0.1, 0.5,
                                8, 10,
                                10^-4,
                                10, 10^6,
                                2.5, 2, 0.2, 0.2, 0.4,
                                TrajOptSOCPs.NONE)
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
if runplots

    hDual = heatmap(hcat(trajStatesAllPD[end].duals))
    yflip!()
    ylabel!("Column")
    title!("Dual Vector for Dynamics Constraints")
    display(hDual)

    # Get the parsed list of trajectories
    nDim = size(grav, 1)
    pltDict = batchPlot(trajStatesAll, cMRocket, nDim, penaltyStart,
                                    thrustMax, thrustAngleMax)
end

# Blocked so that it can be run independently after the fact
if runplots && saveplots
    header = "3DTest_SlantAngle_" * string(currSolveParams.maxOuterIters) *
             "Outer_" * string(currSolveParams.maxNewtonSteps) * "Newton" *
             string(Int64(rocketStart[2 * nDim])) * "Vel" * "_"
    saveBulk(pltDict, header)
end
