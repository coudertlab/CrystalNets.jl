# Data from Blue Obelisk's data repository, corrected for H, C, O, N, S and F
# (see https://github.com/chemfiles/chemfiles/issues/301#issuecomment-574100048)

const atomic_number = Dict{Symbol, Int}(
    :H => 1,
    :He => 2,
    :Li => 3,
    :Be => 4,
    :B => 5,
    :C => 6,
    :N => 7,
    :O => 8,
    :F => 9,
    :Ne => 10,
    :Na => 11,
    :Mg => 12,
    :Al => 13,
    :Si => 14,
    :P => 15,
    :S => 16,
    :Cl => 17,
    :Ar => 18,
    :K => 19,
    :Ca => 20,
    :Sc => 21,
    :Ti => 22,
    :V => 23,
    :Cr => 24,
    :Mn => 25,
    :Fe => 26,
    :Co => 27,
    :Ni => 28,
    :Cu => 29,
    :Zn => 30,
    :Ga => 31,
    :Ge => 32,
    :As => 33,
    :Se => 34,
    :Br => 35,
    :Kr => 36,
    :Rb => 37,
    :Sr => 38,
    :Y => 39,
    :Zr => 40,
    :Nb => 41,
    :Mo => 42,
    :Tc => 43,
    :Ru => 44,
    :Rh => 45,
    :Pd => 46,
    :Ag => 47,
    :Cd => 48,
    :In => 49,
    :Sn => 50,
    :Sb => 51,
    :Te => 52,
    :I => 53,
    :Xe => 54,
    :Cs => 55,
    :Ba => 56,
    :La => 57,
    :Ce => 58,
    :Pr => 59,
    :Nd => 60,
    :Pm => 61,
    :Sm => 62,
    :Eu => 63,
    :Gd => 64,
    :Tb => 65,
    :Dy => 66,
    :Ho => 67,
    :Er => 68,
    :Tm => 69,
    :Yb => 70,
    :Lu => 71,
    :Hf => 72,
    :Ta => 73,
    :W => 74,
    :Re => 75,
    :Os => 76,
    :Ir => 77,
    :Pt => 78,
    :Au => 79,
    :Hg => 80,
    :Tl => 81,
    :Pb => 82,
    :Bi => 83,
    :Po => 84,
    :At => 85,
    :Rn => 86,
    :Fr => 87,
    :Ra => 88,
    :Ac => 89,
    :Th => 90,
    :Pa => 91,
    :U => 92,
    :Np => 93,
    :Pu => 94,
    :Am => 95,
    :Cm => 96,
    :Bk => 97,
    :Cf => 98,
    :Es => 99,
    :Fm => 100,
    :Md => 101,
    :No => 102,
    :Lr => 103,
    :Rf => 104,
    :Db => 105,
    :Sg => 106,
    :Bh => 107,
    :Hs => 108,
    :Mt => 109,
    :Ds => 110,
    :Rg => 111,
    :Cn => 112,
    :Nh => 113,
    :Fl => 114,
    :Mc => 115,
    :Lv => 116,
    :Ts => 117,
    :Og => 118,
    :Uue => 119
)

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

const ismetal = Bool[0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1]



function guess_bonds(pos, types, mat, options)
    # Algorithm from chemfiles, itself from VMD
    @ifwarn begin
        @warn "Guessing bonds with Chemfiles algorithm (from VMD). This may take a while for big structures and may be inexact."
        @info "To avoid guessing bonds, use a file format that contains the bonds."
    end
    bonds = Tuple{Int,Int}[]
    typs = [atomic_number[t] for t in types]
    mof = options.clustering == MOFClustering
    radii = [vdwradii[t]*(1 + mof*ismetal[t]/2) for t in typs]
    cutoff = 3*(options.cutoff_coeff^3.1) * max(maximum(radii), 0.833)
    cutoff2 = 13*options.cutoff_coeff/15
    n = length(pos)
    for i in 1:n
        radius_i = radii[i]
        posi = pos[i]
        typi = types[i]
        skiphomoatomic = typi ∈ options.ignore_homoatomic_bonds
        for j in (i+1):n
            skiphomoatomic && types[j] === typi && continue
            radius_j = radii[j]
            posj = pos[j]
            d1 = periodic_distance(posi, posj, mat)
            d2 = radius_i + radius_j
            if d1 < cutoff && 0.5 < d1 < cutoff2*d2
                push!(bonds, (i,j))
            end
        end
    end
    return bonds
end
