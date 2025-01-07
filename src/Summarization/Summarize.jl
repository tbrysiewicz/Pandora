

function summarize(EP::EnumerativeProblem;kwargs...)
    file_name = "latex_summary"
    rm("OutputFiles/"*file_name*".tex", force=true)
    create_document(file_name*".tex")
    io = open("OutputFiles/"*file_name*".tex", "a");
    write(io, initialize_summary(EP;kwargs...))
    write(io,"\n\n\n\n\n %Abstract\n\n")
    write(io, abstract_summary(EP;kwargs...))
    write(io,"\n\n\n\n\n %Background\n\n")
    write(io, background_summary(EP;kwargs...))
    write(io, degree_summary(EP;kwargs...))
    write(io, end_summary(EP;kwargs...))
    close(io)

    replacement_dictionary = Dict()
    replacement_dictionary["₋"] = ""
    replacement_dictionary["₁"] = "1"
    replacement_dictionary["₂"] = "2"
    replacement_dictionary["₃"] = "3"
    replacement_dictionary["₄"] = "4"
    replacement_dictionary["₅"] = "5"
    replacement_dictionary["₆"] = "6"
    replacement_dictionary["₇"] = "7"
    replacement_dictionary["₈"] = "8"
    replacement_dictionary["₉"] = "9"
    replacement_dictionary["₀"] = "0"

    doc = readlines("OutputFiles/"*file_name*".tex")
    for k in keys(replacement_dictionary)
        doc = map(x->replace(x,k=>replacement_dictionary[k]),doc)
    end

    rm("OutputFiles/"*file_name*".tex", force=true)
    create_document(file_name*".tex")
    io = open("OutputFiles/"*file_name*".tex", "a");
    for r in doc
        println(io,r)
    end
    close(io)
end