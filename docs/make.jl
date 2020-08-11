using Documenter, TrajOptSOCPs

# push!(LOAD_PATH,"../../")
# include("..\\src\\TrajOptSOCPs.jl")


# using SolversWithProjections

makedocs(modules = [TrajOptSOCPs],
         sitename = "SOCP Trajectory Optimization Documentation",
         pages = [
                  "Home" => "index.md",
                  "Using the Module" => Any[
                        "Overview" => "UI/ui.md",
                        "(1) Rockets" => "UI/rocket_ui.md",
                        "(2) Objective" => "UI/objective_ui.md",
                        "(3) Constraints" => "UI/constraints_ui.md",
                        "(4) Augmented Lagrangian" => "UI/augLag_ui.md",
                        "(5) Solve" => "UI/solver_ui.md"
                  ]
                 ],
         format = Documenter.HTML(prettyurls = false))

deploydocs(
 repo = "github.com/danineamati/TrajOptSOCPs.git",
 target = "build",
 push_preview = true,
)
