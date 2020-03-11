#!/bin/bash
#=
exec julia --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#

using Chemfiles

function cif2lammps(cif)
    @assert isfile(cif)
    dir, extension = splitext(cif)
    @assert extension == ".cif"
    name = basename(dir)

    cd(dirname(dir))
    newdir = "./generated_input"
    io = IOBuffer()
    if !ispath(newdir)
        error = nothing
        for cutoff in Float64[12.5, 11, 9, 7, 5, 3]
            io = IOBuffer()
            try
                run(pipeline(`lammps-interface -ff UFF4MOF --cutoff $cutoff $name.cif`;
                             stdout=io, stderr=io))
                error = nothing
                break
            catch e
                println("Build for cutoff = $cutoff failed")
                error = e
            end
        end

        println(String(take!(io)))
        if !isnothing(error)
            throw(error)
        end

        mkpath(newdir)
        mv("in.$name", joinpath(newdir, "in.$name"))
        mv("data.$name", joinpath(newdir, "data.$name"))
        #mv("lammpstrj_to_element.txt", joinpath(newdir, "lammpstrj_to_element.txt"))
        cd(newdir)
        open("in.test.$name", "w") do f
            flag = false
            each = eachline("in.$name", keep=true)
            iterate(each) # Remove the log line
            for l in each
                write(f, l)
            end
        end
        cp("in.test.$name", "in.$name"; force=true)
        open("in.test.$name", "a") do f
            println(f, "fix\t1 all nvt temp 298.00 298.00 100.0\\run\t 1")
        end

        open("in.$name", "a") do f
            println(f, """
            timestep 0.2
            thermo_style custom step temp epair emol etotal press ke pe
            thermo\t1000

            min_style\tfire
            minimize\t1.0e-12 1.0e-12 10000 100000



            # compute 1 all temp
            # variable t equal c_1

            # dump\t$(name)_lammpstrj all atom 100 $name.lammpstrj

            dump\t$(name)_langevin_lammpstrj all atom 100 $(name)_langevin.lammpstrj
            fix\tpre1 fram langevin 298.00 298.00 100.0 2325813 zero yes
            fix\tpre2 fram nve
            thermo\t100
            run\t10000
            unfix\tpre1
            unfix\tpre2
            undump $(name)_langevin_lammpstrj

            dump\t$(name)_lammpstrj all atom 100 $name.lammpstrj
            fix\t1 all nvt temp 298.00 298.00 100.0
            # fix\t2 all recenter INIT INIT INIT
            thermo\t10
            run\t100000
            undump\t$(name)_lammpstrj
            """)
        end
        cd("..")
    end
    cd(newdir)
    try
        run(pipeline(`lammps -in in.test.$name`, devnull), wait=true)
        run(`lammps -in in.$name`)
    catch e
        if e isa ProcessFailedException
            println("Not working with GPU")
            run(`lammps -sf opt -in in.$name`)
        else
            rethrow()
        end
    end
end

cif2lammps(ARGS[1])
include("analyze_log.jl")

display(plot_energy("./log.lammps"))
