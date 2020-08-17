# Integrate the Rocket Landing example with MonteCarlo

# RUN:
# ] activate TrajOptSOCPs

using TrajOptSOCPs
using LinearAlgebra, SparseArrays
using Dates

include("../plotting/batchPlots.jl")
include("../plotting/batchPlots_MonteCarlo.jl")
pyplot()
runplots = true
saveplots = true


println("\n--------------------------------------------")
println("          Setting Up Rocket Landing ")
println("--------------------------------------------")

# We will keep the same rocket in all runs
mass = 1
isp = 1
grav = [0; 0; -9.81]
deltaTime = 1
rocket = rocket_simple(mass, isp, grav, deltaTime)

# Then ending point is always a soft landing at the origin.
# But, we want to vary the starting point.
minStateStart = [-10.0; -5.0; 20.0; -4.0; -4.0; -8.0]
maxStateStart = [  0.0;  5.0; 30.0;  4.0;  4.0;  0.0]
const numRuns = 50

# Stores each randomly generated starting case
startList = TrajOptSOCPs.montecarlo(minStateStart, maxStateStart, numRuns)

# Maintain the same end point
const rocketEnd = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0]

# Initialize the hover control, which is the base control to prevent freefall.
uHover = mass * grav

# Number of time steps to discretize the trajectory
const NSteps = 60

# The following trajectories do not change with starting condition:
# ----------------------------
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

# Set up the LQR matrices
const lqrQMat = 0.0 * Diagonal(I, size(rocketEnd, 1))
const lqrRMat = 0.025 * Diagonal(I, Int64(size(rocketEnd, 1) / 2))

# Initialize the Augmented Lagrangian penalty on the constraints
penaltyStart = 1.0

# Next we select resonable solver parameters
currSolveParams = solverParams(0.1, 0.5,
                                8, 10,
                                10^-4,
                                10, 10^6,
                                2.5, 2, 0.2, 0.2, 0.4,
                                TrajOptSOCPs.NONE)
# ----------------------------

solParamPrint(currSolveParams)
println()
println()


# Arrays to store the results
trajPDFinals = []
trajPFinals = []

dynVio = []
constVio = []


# Now we loop through the different starting conditions and solve the
# optimization problem
for (ind, rocketStart) in enumerate(startList)
    # Initialize the trajectory with a line
    initTraj = initializeTraj(rocketStart, rocketEnd, uHover, uHover, NSteps)

    # Use a Linear Quadratic Regulator as the cost function
    costFun = makeLQR_TrajReferenced(lqrQMat, lqrRMat, NSteps, initTraj)

    # Create the Dynamics Constraints
    ADyn, BDyn = rocketDynamicsFull(rocket, rocketStart, rocketEnd, NSteps)
    dynConstraint = AL_AffineEquality(ADyn, BDyn)
    lambdaInit = -1 * ones(size(BDyn))

    # Create the constraint manager
    cMRocket = constraintManager_Dynamics(
                [groundConstraint, maxThrustConstraint, maxAngleConstraint],
                [groundLambda, maxThrustLambda, maxAngleLambda],
                dynConstraint, lambdaInit
                )

    # Initialize the primal-dual vector
    initTrajPD = [initTraj; lambdaInit]

    # Equiped with the constraint term and the objective term, I now build the
    # Augmented Lagrangian
    alRocket = augLag(costFun, cMRocket, penaltyStart)

    println("--------------------------------------------")
    println("             Beginning Solve $ind           ")
    println("--------------------------------------------")

    trajLambdaSolved, resArr = ALPrimalNewtonMain(initTrajPD, alRocket,
                                            currSolveParams)
    trajPD = parsePrimalDualVec(trajLambdaSolved[end], size(initTraj, 1))
    trajP = trajPD.primals

    push!(trajPDFinals, trajPD)
    push!(trajPFinals, trajP)

    # Next get the final constraint violations
    # (1) Dynamics Vioaltion
    push!(dynVio, safeNorm(norm(evalAffineEq(cMRocket, trajP))))
    # (2) Other Constraint Violation
    push!(constVio, safeNorm(evalConstraints(cMRocket, trajP, penaltyStart)))

end

if runplots
    println("--------------------------------------------")
    println("             Beginning Plotting             ")
    println("--------------------------------------------")

    allPlots = batchPlotMonteCarlo(trajPFinals, dynVio, constVio,
                            size(grav, 1), penaltyStart,
                            thrustMax, thrustAngleMax)
end

if runplots && saveplots
    curr = Dates.format(now(), "yyyymmdd_HH-MM-SS")
    header = "MonteCarlo_$(numRuns)Runs_" * curr * "_"
    saveBulk(allPlots, header)
end
