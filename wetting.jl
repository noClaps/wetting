using Roots
using DelimitedFiles
using ColorSchemes
import Contour
using Plots
pythonplot()

# Equation whose root defines the homogeneous equilibrium solutions:
function func(beta, mu)
    return (rho) -> rho - (1.0 - rho) * exp(beta * (mu + 5.0 * rho))
end

function problem14(beta_epsilon_wall)
    Lx = 100 # Number of sites along x-axis
    Ly = 40 # Number of sites along y-axis

    beta = 1.2 # Inverse temperature beta*epsilon
    epsilon_wall = beta_epsilon_wall / beta

    mu_coex = -5 / 2
    tol = 1e-12 # Convergence tolerance
    count = 300000 # Upper limit for iterations
    alpha = 0.01 # Mixing parameter

    # Solve equations iteratively:
    conv = 1
    cnt = 1

    minsol = minimum(fzeros(func(beta, mu_coex), 0, 1))
    maxsol = maximum(fzeros(func(beta, mu_coex), 0, 1))

    rho = ones(Lx, Ly) * minsol

    squarewidth = 30
    squareheight = 20
    shape = ones(squarewidth, squareheight) * maxsol
    rho[Int(Lx / 2)-Int(squarewidth / 2)+1:Int(Lx / 2)+Int(squarewidth / 2), 2:squareheight+1] = shape

    rho[:, begin] .= 0

    rho_initial = transpose(rho)
    rho_new = zeros(Lx, Ly)

    while conv >= tol && cnt < count
        cnt += 1
        for i = 1:Lx
            for j = 2:Ly-1
                # Handle the periodic boundaries for x and y:
                left = i == 1 ? Lx : i - 1
                right = i == Lx ? 1 : i + 1
                down = j == 1 ? Ly : j - 1
                up = j == Ly ? 1 : j + 1

                v_j = -epsilon_wall * (j - 1)^(-3)
                rho_new[i, j] = (1 - rho[i, j]) * exp(beta * (rho[i, down] + rho[i, up] + rho[left, j] + rho[right, j] + (1 / 4) * (rho[left, down] + rho[right, down] + rho[left, up] + rho[right, up]) + mu_coex - v_j))
            end
        end

        conv = sum((rho - rho_new) .^ 2) # Compute the convergence parameter.
        rho = alpha * rho_new + (1 - alpha) * rho # Mix the new and old solutions.

        # Normalization step
        N = sum(rho_initial)
        N_current = sum(rho)
        rho *= N / N_current

        rho[:, begin] .= 0
        rho[:, end] .= minsol
    end

    rho = transpose(rho)

    writedlm("problem14-$beta_epsilon_wall.csv", rho, ",")
end

# for beta = 0.1:0.1:2.0
#     problem14(beta)
# end

rho = readdlm("problem14-0.5.csv", ',')
x = 1:size(rho, 1)
y = 1:size(rho, 2)

colorscheme = cgrad(ColorScheme(get(ColorSchemes.jet, range(0.2, 0.7, 256))))
# heatmap(rho, color=colorscheme)
contourf(rho, color=colorscheme, title="Contour plot", xlabel="x positions", ylabel="y positions", colorbar_title="Density gradient", size=(1000, 400), levels=256)

level = Contour.levels(Contour.contours(x, y, rho, 1))[1]
line = first(Contour.lines(level))
ys, xs = Contour.coordinates(line)

indices = []
grads = []

for i = 1:length(xs)-1
    gradient = (ys[i+1] - ys[i]) / (xs[i+1] - xs[i])
    angle = 180 - rad2deg(atan(gradient))

    if floor(angle) == 115
        push!(indices, i)
        push!(grads, gradient)
    end
end

for i = eachindex(indices)
    index = indices[i]
    grad = grads[i]

    local x = xs[index]:(xs[index]+5)
    local y = (grad .* (x .- xs[index])) .+ ys[index]
    plot!(x, y, color=:black, legend=false)
end

gui()
readline()
