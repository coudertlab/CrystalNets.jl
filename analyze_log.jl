using Plots

function find_thermo!(each)
    x = iterate(each)
    isnothing(x) && return ""
    l,_ = x
    while length(l) < 3 || l[1:3] != "run"
        x = iterate(each)
        isnothing(x) && return ""
        l,_ = x
    end
    while length(l) < 4 || l[1:4] != "Step"
        x = iterate(each)
        isnothing(x) && return ""
        l,_ = x
    end
    return l
end

function parse_line(l)
    ret = Float64[]
    next = true
    start = 0
    for i in 1:length(l)
        if l[i] == ' '
            if !next
                next = true
                push!(ret, parse(Float64, l[start:i]))
            end
        else
            if next
                start = i
                next = false
            end
        end
    end
    if !next
        push!(ret, parse(Float64, l[start:end]))
    end
    return ret
end


function plot_energy(log)
    step = Int[]
    eng = Float64[]
    pot = Float64[]
    kin = Float64[]
    open(log, "r") do f
        each = eachline(f)
        _thermo = find_thermo!(each)
        while _thermo != ""
            thermo = split(_thermo, ' ')
            i_step = findfirst(==("Step"), thermo)
            i_eng = findfirst(==("TotEng"), thermo)
            i_kin = findfirst(==("KinEng"), thermo)
            i_pot = findfirst(==("PotEng"), thermo)
            i_max = max(i_step, i_eng, i_kin, i_pot)
            l = readline(f)
            while !isempty(l) && !isletter(l[1])
                parsed = parse_line(l)
                length(parsed) >= i_max || break
                @inbounds push!(step, Int(parsed[i_step]))
                @inbounds push!(eng, parsed[i_eng])
                @inbounds push!(kin, parsed[i_kin])
                @inbounds push!(pot, parsed[i_pot])
                l = readline(f)
            end
            _thermo = find_thermo!(each)
        end
    end
    p = scatter(step, eng, marker=:+, markersize=1.2, label="Total energy")
    scatter!(p, step, kin, marker=:+, markersize=1.2, label="Kinetic energy")
    scatter!(p, step, pot, marker=:+, markersize=1.2, label="Potential energy")
    return p
end
