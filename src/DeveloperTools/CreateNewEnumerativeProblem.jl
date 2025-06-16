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

    # Public-facing function name is lowercase snake of raw name
    public_fn = to_lower_snake(raw_prop_name)

    # Prompt for the input property names (comma-separated, e.g. SYSTEM, ANOTHER_PROP)
    input_props_str = prompt("Enter input EnumerativeProperties required (comma-separated, e.g. SYSTEM):")
    input_props = [strip(p) for p in split(input_props_str, ',') if !isempty(strip(p))]

    # Prompt for algorithm metadata
    alg_name = prompt("Enter the name of the algorithm (for documentation):")
    alg_desc = prompt("Enter a short description of the algorithm:")
    reliability = prompt("Enter reliability (:certified, :heuristic, :null, etc.):")

    # Prompt for function input type (usually the first input property)
    fn_input_type = ""
    if !isempty(input_props)
        fn_input_type = input_props[1] == "SYSTEM" ? "System" : "TODO_Type_for_" * input_props[1]
    else
        fn_input_type = "TODO_InputType"
    end

    # Compose the code block
    println("\n########## Implementation of Enumerative Property ",prop_name,"#########\n")

    println("const $prop_name = EnumerativeProperty{$(Symbol(fn_input_type))}(\"$prop_string\")\n")

    println("""\"\"\"
    $public_fn(EP::EnumerativeProblem; kwargs...)

Return the $alg_name of the enumerative problem.
\"\"\"
function $public_fn(EP::EnumerativeProblem; kwargs...)
    $prop_name(EP; kwargs...)
end
""")

    println("""function compute_$public_fn(#::$fn_input_type)::$fn_input_type
    # TODO: Implement core algorithm here
end
""")

    input_props_joined = join(input_props, ", ")
    println("""$(public_fn)_datum = AlgorithmDatum(
    name = \"$alg_name\",
    description = \"$alg_desc\",
    input_properties = [$input_props_joined],
    output_property = $prop_name,
    reliability = $reliability
)

ALGORITHM_DATA[compute_$public_fn] = $(public_fn)_datum
""")
end

generate_algorithm_block()
