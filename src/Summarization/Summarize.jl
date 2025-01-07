

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
    write(io, end_summary(EP;kwargs...))
    close(io)
end