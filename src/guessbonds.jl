# Bond guessing algorithm derived from from Chemfiles, itself from VMD

# Data from Blue Obelisk's data repository, corrected for H, C, O, N, S and F
# (see https://github.com/chemfiles/chemfiles/issues/301#issuecomment-574100048)
const vdwradii = Float32[1.0, 1.4, 2.2, 1.9, 1.8, 1.5, 1.4, 1.3, 1.2, 1.54, 2.4,
                         2.2, 2.1, 2.1, 1.95, 1.9, 1.8, 1.88, 2.8, 2.4, 2.3,
                         2.15, 2.05, 2.05, 2.05, 2.05, 2.0, 2.0, 2.0, 2.1, 2.1,
                         2.1, 2.05, 1.9, 1.9, 2.02, 2.9, 2.55, 2.4, 2.3, 2.15,
                         2.1, 2.05, 2.05, 2.0, 2.05, 2.1, 2.2, 2.2, 2.25, 2.2,
                         2.1, 2.1, 2.16, 3.0, 2.7, 2.5, 2.48, 2.47, 2.45, 2.43,
                         2.42, 2.4, 2.38, 2.37, 2.35, 2.33, 2.32, 2.3, 2.28,
                         2.27, 2.25, 2.2, 2.1, 2.05, 2.0, 2.0, 2.05, 2.1, 2.05,
                         2.2, 2.3, 2.3, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.4, 2.0,
                         2.3, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
                         2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
                         2.0, 2.0, 2.0, 2.0, 2.0, 2.0]

# Data from PeriodicTable.jl
const ismetal = Bool[0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1]

const ismetalloid = Bool[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

"""
    guess_bonds(pos, types, mat, options)

Return the bonds guessed from the positions, types and cell matrix, given as a
`Vector{Vector{Tuple{Int,Float32}}}`.

The `i`-th entry of the list is a list, whose entries are of the form `(j, dist)` which
indicates that the representatives of vertices `i` and `j` distant of at most `dist` are
bonded together.
"""
function guess_bonds(pos, types, mat, options)
    # Algorithm from chemfiles, itself from VMD
    @ifwarn if options.bonding != Bonding.Guess
        @warn "Guessing bonds with custom algorithm (from Chemfiles and VMD). This may take a while for big structures and may be inexact."
        @info "To avoid guessing bonds, use a file format that contains the bonds."
    end
    n = length(pos)
    bonds = [Tuple{Int,Float32}[] for _ in 1:n]
    @toggleassert n == length(types)
    radii = Vector{Float32}(undef, n)
    for (i, typ) in enumerate(types)
        t = get(atomic_numbers, typ, nothing)
        if t isa Int
            radii[i] = vdwradii[t]*(1 + options.wider_metallic_bonds*(ismetal[t]|ismetalloid[t])*0.5)
        else
            radii[i] = 0.0
            @ifwarn @warn lazy"Unrecognized atom type \"$t\" will be considered a dummy atom."
        end
    end
    cutoff = 3*(options.cutoff_coeff^3.1) * max(maximum(radii), 0.833)
    cutoff2 = 13*options.cutoff_coeff/15
    buffer, ortho, safemin = prepare_periodic_distance_computations(mat)
    for i in 1:n
        radius_i = radii[i]
        iszero(radius_i) && continue
        posi = pos[i]
        typi = types[i]
        #ignoreifmetallic = typi === :C
        #ignoreifC = ismetal[atomic_numbers[typi]]
        skiphomoatomic = typi ∈ options.ignore_homoatomic_bonds ||
                         (options.ignore_homometallic_bonds && ismetal[atomic_numbers[typi]])
        acceptonlyO = options.structure === StructureType.Zeolite && typi !== :O
        acceptallbutO = options.structure === StructureType.Zeolite && typi === :O
        for j in (i+1):n
            typj = types[j]
            skiphomoatomic && typj === typi && continue
            acceptonlyO && typj !== :O && continue
            acceptallbutO && typj === :O && continue
            #(ignoreifC & (typj === :C)) && continue
            #ignoreifmetallic && ismetal[atomic_numbers[typj]] && continue
            radius_j = radii[j]
            iszero(radius_j) && continue
            posj = pos[j]
            d1 = periodic_distance!(buffer, posi .- posj, mat, ortho, safemin)
            maxdist = cutoff2*(radius_i + radius_j)
            if d1 < cutoff && 0.5 < d1 < maxdist
                push!(bonds[i], (j, maxdist))
                push!(bonds[j], (i, maxdist))
            end
        end
    end
    return bonds
end
