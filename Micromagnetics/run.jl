#=
    Calculates the magnetostatic field (demagnetizing field)
    for a given magnetization, using FEM-BEM.

    All the FEM-BEM calculations are done in C++
=#

include("../src/gmsh_wrapper.jl")

# For plots | Uncomment the plot section of "main()"
using GLMakie

function main(meshSize=0,showGmsh=true,saveMesh=false)
    #=
        This creates the mesh. C++ then handles the simulation
    =#
    
    mu0 = pi*4e-7                       # vacuum magnetic permeability
    Hext::Vector{Float64} = [1,0,0]     # T

    # Create a geometry
    gmsh.initialize()

    # >> Model
    L::Vector{Float64} = [100, 100, 5] # Dimensions of the object

    # Add the body
    addCuboid([0,0,0],L)

    # Generate Mesh
    mesh = Mesh([],meshSize,0.0,saveMesh)
    
    if showGmsh
        gmsh.fltk.run()
    end
    gmsh.finalize()

    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))
    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))

    # Adjust to 0 indexing
    t::Matrix{Int32} = mesh.t .- 1
    surfaceT::Matrix{Int32} = mesh.surfaceT .- 1

    # Magnetization
    m::Matrix{Float64} = zeros(3,mesh.nv)
    m[1,:] .= 1

    maxAtt::Int32 = 15000
    M_avg::Matrix{Float64} = zeros(3,maxAtt)

    @time begin
    # Calculate the magnetostatic field
    @ccall "julia_wrapper.so".LL(
        m::Ptr{Float64},
        mesh.p::Ptr{Float64},
        t::Ptr{Int32},
        surfaceT::Ptr{Int32},
        mesh.normal::Ptr{Float64},
        mesh.AE::Ptr{Float64},
        mesh.VE::Ptr{Float64},
        mesh.nv::Int32,
        mesh.nt::Int32,
        mesh.ne::Int32,
        M_avg::Ptr{Float64},
        maxAtt::Int32
    )::Cvoid
    end

    time::Vector{Int32} = 1:maxAtt

    fig = Figure()
    ax = Axis( fig[1,1] )

    scatter!(ax,time,M_avg[1,:], label = "M_x")
    scatter!(ax,time,M_avg[2,:], label = "M_y")
    scatter!(ax,time,M_avg[3,:], label = "M_z")
    axislegend()

    # save("M_time_Sphere.png",fig)
    wait(display(fig))

    # save("H.png",fig)

end # end of main

meshSize = 0
showGmsh = false
saveMesh = false

main(meshSize,showGmsh,saveMesh)

