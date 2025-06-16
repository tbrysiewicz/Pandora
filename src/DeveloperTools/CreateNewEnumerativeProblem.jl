#Copy and paste this code into a new julia terminal and answer the questions to create a boilerplate for a new enumerative problem
function prompt(msg)
    print(msg, " ")
    chomp(readline())
end

function to_all_caps(name::AbstractString)
    uppercase(replace(String(name), r"\s+" => "_"))
end

function to_lower_snake(name::AbstractString)
    lowercase(replace(String(name), r"\s+" => "_"))
end

function generate_algorithm_block()
    println("=== Enumerative Algorithm Code Generator ===")

    # Prompt for EnumerativeProperty name and transform to styles
    raw_prop_name = prompt("Enter the EnumerativeProperty name (e.g. Blah):")
    prop_name = to_all_caps(raw_prop_name)
    prop_string = to_lower_snake(raw_prop_name)    
    prop_type = prompt("Enter the output type of the EnumerativeProperty (e.g. Int, Float64, Group):")

    # Public-facing function name is lowercase snake of raw name
    public_fn = to_lower_snake(raw_prop_name)

    # Prompt for the input property names (comma-separated, e.g. SYSTEM, ANOTHER_PROP)
    input_props_str = prompt("Enter input EnumerativeProperties required (comma-separated, e.g. SYSTEM):")
    input_props = [to_all_caps(strip(p)) for p in split(input_props_str, ',') if !isempty(strip(p))]

    # Prompt for algorithm metadata
    alg_name = prop_string
    alg_desc = "#TODO: Enter a short description of the algorithm:"
    reliability = "#TODO: Enter reliability (:certified, :heuristic, :null, etc.):"

    # Generate a list of argument strings like "#1::System", "#2::Group", etc.
    input_arg_strings = [
        begin
            EProp = getfield(Pandora, Symbol(p))
            argtype = Pandora.get_type(EProp)
            "#$i$EProp::$(string(argtype))"
        end
        for (i, p) in enumerate(input_props)
    ]
    fn_arg_list = join(input_arg_strings, ", ")
    

    # Combine them into the function argument list
    fn_arg_list = join(input_arg_strings, ", ")

    # Compose the code block
    println("\n########## Implementation of Enumerative Property ",prop_name,"#########\n")

    println("const $prop_name = EnumerativeProperty{$(Symbol(prop_type))}(\"$prop_string\")\n")

    println("""\"\"\"
    $public_fn(EP::EnumerativeProblem; kwargs...)

Return the $public_fn of the enumerative problem.
\"\"\"
function $public_fn(EP::EnumerativeProblem; kwargs...)
    $prop_name(EP; kwargs...)
end
""")

    println("""function compute_$public_fn($fn_arg_list)::$prop_type
    # TODO: Implement core algorithm here
end
""")

    input_props_joined = join(input_props, ", ")
    println("""compute_$(public_fn)_datum = AlgorithmDatum(
    name = \"$alg_name\",
    description = \"$alg_desc\",
    input_properties = [$input_props_joined],
    output_property = $prop_name,
    reliability = $reliability
)

ALGORITHM_DATA[compute_$public_fn] = compute_$(public_fn)_datum
""")
end

generate_algorithm_block()
