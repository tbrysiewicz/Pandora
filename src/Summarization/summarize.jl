export 
    summarize

struct LatexDocumentMetadata
    title::String
    author::String
    date::String
end

function LatexDocumentMetadata(; 
        title="Pandora.jl summary of enumerative problem", 
        author="", 
        date=string(today()))
    full_author = isempty(author) ? "Pandora.jl" : "$author and Pandora.jl"
    return LatexDocumentMetadata(title, full_author, date)
end

function latex_metadata_block(meta::LatexDocumentMetadata)::String
    return join([
        "\\title{" * meta.title * "}",
        "\\author[" * meta.author * "]{" * meta.author * "}",
        "\\date{" * meta.date * "}"
    ], "\n") * "\n"
end

function generate_latex_header(meta::LatexDocumentMetadata)::String
    preamble = read("Latex/Preamble.txt", String)
    return preamble * "\n" * latex_metadata_block(meta) * "\\begin{document}\n\\maketitle\n"
end

function sanitize_latex_text(s::String)::String
    replacement_dictionary = Dict(
        "₋" => "",  "₁" => "1", "₂" => "2", "₃" => "3", "₄" => "4",
        "₅" => "5", "₆" => "6", "₇" => "7", "₈" => "8", "₉" => "9", "₀" => "0"
    )
    for (k, v) in replacement_dictionary
        s = replace(s, k => v)
    end
    return s
end

include("abstract_summary.jl")
include("background_summary.jl")
include("degree_summary.jl")
include("reality_summary.jl")
include("monodromy_summary.jl")
include("bibliography.jl")

function compile_latex_to_pdf(tex_path::String; output_dir::String = dirname(tex_path))
    run(`pdflatex -output-directory=$output_dir $tex_path`)
    run(`pdflatex -output-directory=$output_dir $tex_path`)
end


function summarize(EP::EnumerativeProblem;
    metadata::LatexDocumentMetadata = LatexDocumentMetadata(),
    output_path::String = "OutputFiles/latex_summary.tex",
    section_list::Vector{Function} = [abstract_summary, background_summary, degree_summary, reality_summary, monodromy_summary, bibliography_summary],
    compile_pdf::Bool = true)

    # Generate the LaTeX content
    latex_string = generate_latex_header(metadata)

    for section_fn in section_list
       latex_string *= "\n\n" * section_fn(EP)
    end

    latex_string *= "\n\\end{document}\n"

    # Sanitize LaTeX for common unicode issues
    latex_string = sanitize_latex_text(latex_string)

    # Write to file
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        write(io, latex_string)
    end
    if compile_pdf
        compile_latex_to_pdf(output_path)
    end
    return output_path
end
