function initialize_summary(EP::EnumerativeProblem; kwargs...)
    preamble = latex_preamble()
    meta_data = latex_metadata(;kwargs...)
    begin_document = raw"""
    \begin{document}

    \maketitle"""
    return(preamble*meta_data*begin_document)
end

function latex_preamble()
    preamble = open("Latex/Preamble.txt","r")
    s = read(preamble,String)
    return(s)
end


function latex_metadata(;kwargs...)
    title = get(kwargs,:title,"Pandora.jl summary of enumerative problem")
    date = string(today())
    author = get(kwargs,:author,"")
    author = author*" and Pandora.jl"

    meta_strings = []

	title_string = raw"""\title{"""
    title_string = title_string*title*raw"""}"""
    push!(meta_strings,title_string)

	author_string = raw"""\author["""
    author_string = author_string*author*raw"""]{"""
    author_string = author_string*author*raw"""}"""
    push!(meta_strings,author_string)

    date_string = raw"""\date{"""
    date_string = date_string*date*raw"""}"""
    push!(meta_strings,date_string)

    meta_string = ""
    for ms in meta_strings
        meta_string = meta_string*ms*"\n"
    end

    return(meta_string)
end

