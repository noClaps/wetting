using ColorSchemes
using DelimitedFiles
using Plots
pythonplot()

colorscheme = cgrad(ColorScheme(get(ColorSchemes.jet, range(0.2, 0.7, 256))))

for bew = 0.1:0.1:2.0
    rho = readdlm("problem14-$bew.csv", ',')
    p = contourf(rho, color=colorscheme, title="Droplet shape on 100x40 lattice, βε_w = $bew", xlabel="x positions", ylabel="y positions", colorbar_title="Density gradient", size=(100, 40) .* 12.5, levels=256)

    # change the path in the second argument to wherever you want to save the images to
    png(p, "/Users/user/Downloads/bew$bew.png")
end
