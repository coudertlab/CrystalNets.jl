# Data from Blue Obelisk's data repository, corrected for H, C, O, N, S and F
# (see https://github.com/chemfiles/chemfiles/issues/301#issuecomment-574100048)
const vdwradii = Dict{Symbol, Float64}(
    :H => 1.0,
    :He => 1.4,
    :Li => 2.2,
    :Be => 1.9,
    :B => 1.8,
    :C => 1.5,
    :N => 1.4,
    :O => 1.3,
    :F => 1.2,
    :Ne => 1.54,
    :Na => 2.4,
    :Mg => 2.2,
    :Al => 2.1,
    :Si => 2.1,
    :P => 1.95,
    :S => 1.9,
    :Cl => 1.8,
    :Ar => 1.88,
    :K => 2.8,
    :Ca => 2.4,
    :Sc => 2.3,
    :Ti => 2.15,
    :V => 2.05,
    :Cr => 2.05,
    :Mn => 2.05,
    :Fe => 2.05,
    :Co => 2.0,
    :Ni => 2.0,
    :Cu => 2.0,
    :Zn => 2.1,
    :Ga => 2.1,
    :Ge => 2.1,
    :As => 2.05,
    :Se => 1.9,
    :Br => 1.9,
    :Kr => 2.02,
    :Rb => 2.9,
    :Sr => 2.55,
    :Y => 2.4,
    :Zr => 2.3,
    :Nb => 2.15,
    :Mo => 2.1,
    :Tc => 2.05,
    :Ru => 2.05,
    :Rh => 2.0,
    :Pd => 2.05,
    :Ag => 2.1,
    :Cd => 2.2,
    :In => 2.2,
    :Sn => 2.25,
    :Sb => 2.2,
    :Te => 2.1,
    :I => 2.1,
    :Xe => 2.16,
    :Cs => 3.0,
    :Ba => 2.7,
    :La => 2.5,
    :Ce => 2.48,
    :Pr => 2.47,
    :Nd => 2.45,
    :Pm => 2.43,
    :Sm => 2.42,
    :Eu => 2.4,
    :Gd => 2.38,
    :Tb => 2.37,
    :Dy => 2.35,
    :Ho => 2.33,
    :Er => 2.32,
    :Tm => 2.3,
    :Yb => 2.28,
    :Lu => 2.27,
    :Hf => 2.25,
    :Ta => 2.2,
    :W => 2.1,
    :Re => 2.05,
    :Os => 2.0,
    :Ir => 2.0,
    :Pt => 2.05,
    :Au => 2.1,
    :Hg => 2.05,
    :Tl => 2.2,
    :Pb => 2.3,
    :Bi => 2.3,
    :Po => 2.0,
    :At => 2.0,
    :Rn => 2.0,
    :Fr => 2.0,
    :Ra => 2.0,
    :Ac => 2.0,
    :Th => 2.4,
    :Pa => 2.0,
    :U => 2.3,
    :Np => 2.0,
    :Pu => 2.0,
    :Am => 2.0,
    :Cm => 2.0,
    :Bk => 2.0,
    :Cf => 2.0,
    :Es => 2.0,
    :Fm => 2.0,
    :Md => 2.0,
    :No => 2.0,
    :Lr => 2.0,
    :Rf => 2.0,
    :Db => 2.0,
    :Sg => 2.0,
    :Bh => 2.0,
    :Hs => 2.0,
    :Mt => 2.0,
    :Ds => 2.0,
    :Rg => 2.0,
    :Cn => 2.0,
    :Uut => 2.0,
    :Fl => 2.0,
    :Uup => 2.0,
    :Lv => 2.0,
    :Uus => 2.0,
    :Uuo => 2.0
)

function guess_bonds(pos, types, mat)
    # Algorithm from chemfiles, itself from VMD
    @ifwarn begin
        @warn "Guessing bonds with Chemfiles algorithm (from VMD). This may take a while for big structures and may be inexact."
        @info "To avoid guessing bonds, use a file format that contains the bonds."
    end
    bonds = Tuple{Int,Int}[]
    radii = [vdwradii[t] for t in types]
    cutoff = 1.2 * max(maximum(radii), 0.833)
    n = length(pos)
    for i in 1:n
        radius_i = radii[i]
        posi = pos[i]
        for j in (i+1):n
            radius_j = radii[j]
            posj = pos[j]
            d1 = periodic_distance(posi, posj, mat)
            d2 = radius_i + radius_j
            if d1 < cutoff && 0.03 < d1 < 0.6*d2
                push!(bonds, (i,j))
            end
        end
    end
    return bonds
end
