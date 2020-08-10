using LinearAlgebra

function backtrackLineSearch(xInit, dirΔ, f, dfdx, paramA, paramB, verbose = false)
    # Output an updated x

    if verbose
        println("\nBeginning Line Search")
        println("xInit = $xInit")
        println("dirΔ = $dirΔ")
        println("objective function (current x) = $(f(xInit))")
        println("Parameters: $paramA and $paramB")
        println("Verbose = $verbose \n")
    end

    # First, initializae α = 1. This is the "learning rate"

    α = 1
    if verbose
        print("Left side: ")
        println(f(xInit + α * dirΔ))
        print("Right side: ")
        print(f(xInit))
        print(" + ")
        print(paramA * α * dirΔ'dfdx(xInit))
        print(" = ")
        println(f(xInit) + paramA * α * dirΔ'dfdx(xInit))
    end

    # Check the condition:
    # f(x + alpha * dir) > f(x) + a * grad f of (alpha * delta)
    while f(xInit + α * dirΔ) > (f(xInit) + paramA * α * (dirΔ'dfdx(xInit))[1])
        # If the condition passes, set α to b * α
        if verbose
            print("Left side: ")
            println(f(xInit + α * dirΔ))
            print("Right side: ")
            print(f(xInit))
            print(" + ")
            print(paramA * α * dirΔ'dfdx(xInit))
            print(" = ")
            println(f(xInit) + paramA * α * dirΔ'dfdx(xInit))
        end

        α = paramB * α

        if verbose
            println(α)
        end
    end

    # the updated x is xInit - alpha * gradf(xInit)
    xUpdated = xInit + α * dirΔ
    # Return updated x
    return xUpdated, α
end


runTest = false

if runTest
    # Test Script

    # One Dimensional Function
    xInit1 = 600
    fFun(x) = x^2
    dfdx(x) = 2 * x
    dir = -dfdx(xInit1)
    aTest = 0.3
    bTest = 0.707

    xNew = backtrackLineSearch(xInit1, dir, fFun, dfdx, aTest, bTest, true)
    println("New x = $xNew")
    println("Result should by 0.181")


    # Two Dimensional Function
    xInit2 = [600, 400]
    fFun(x) = x'x
    dfdx(x) = 2 * x
    dir2 = -dfdx(xInit2)
    aTest = 0.3
    bTest = 0.707

    xNew = backtrackLineSearch(xInit2, dir2, fFun, dfdx, aTest, bTest, true)
    println("New x = $xNew")
    println("Result should by 0.181")
end
