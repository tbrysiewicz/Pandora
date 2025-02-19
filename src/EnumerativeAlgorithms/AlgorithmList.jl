
monodromy_group_algorithm_1 = EnumerativeAlgorithm(
    "monodromy from monodromy dict",
    [monodromy_dictionary],
    Vector{Symbol}(),
    function get_group_generated(md)
        G = subgroup(unique(keys(md)))
    end,
    :monodromy_group,
    NullCitation,
    :uncertified
)

global MainAlgorithms = [monodromy_group_algorithm_1]