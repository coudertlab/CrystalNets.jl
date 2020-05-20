const general_positions = deserialize("general_positions")::Vector{Vector{EquivalentPosition}}
# @test all(parse.(EquivalentPosition, string.(x)) == x for x in CIFTypes.general_positions)

