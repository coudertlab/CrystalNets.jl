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

# Si is considered not to be a metalloid below to allow Si-Si bonds
const ismetalormetalloid = Bool[0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
                                0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 0, 0, 1]


function guess_bonds(pos, types, mat, options)
    # Algorithm from chemfiles, itself from VMD
    @ifwarn begin
        @warn "Guessing bonds with Chemfiles algorithm (from VMD). This may take a while for big structures and may be inexact."
        @info "To avoid guessing bonds, use a file format that contains the bonds."
    end
    bonds = Tuple{Int,Int,Float64}[]
    typs = [atomic_numbers[t] for t in types]
    mof = options.clustering_mode == ClusteringMode.MOF ||
          options.clustering_mode == ClusteringMode.Guess
    radii = [vdwradii[t]*(1 + mof*ismetalormetalloid[t]*0.5) for t in typs]
    cutoff = 3*(options.cutoff_coeff^3.1) * max(maximum(radii), 0.833)
    cutoff2 = 13*options.cutoff_coeff/15
    n = length(pos)
    for i in 1:n
        radius_i = radii[i]
        posi = pos[i]
        typi = types[i]
        skiphomoatomic = typi ∈ options.ignore_homoatomic_bonds ||
                         (options.ignore_homometallic_bonds &&
                          ismetalormetalloid[atomic_numbers[typi]])
        for j in (i+1):n
            skiphomoatomic && types[j] === typi && continue
            radius_j = radii[j]
            posj = pos[j]
            d1 = periodic_distance(posi, posj, mat)
            maxdist = cutoff2*(radius_i + radius_j)
            if d1 < cutoff && 0.5 < d1 < maxdist
                push!(bonds, (i, j, maxdist))
            end
        end
    end
    return bonds
end
