export
    short


##############################################################
######################## Citation ############################
##############################################################

"""
A `Citation` is a type that represents a citation for an algorithm.
It contains the authors, title, journal, and year of the publication.
"""
@kwdef struct Citation
    authors::Vector{String} = [""]
    title::String = ""
    journal::String = ""
    year::Int = 1234
    edition::String = ""
    publisher::String = ""
    url::String = ""
    pages::String = ""
    doi::String = ""
    note::String = ""
end

const NULL_CITATION = Citation()

const Dummit04 = Citation(
    authors = ["Dummit", "Foote"],
    title = "Abstract Algebra",
    publisher = "John Wiley & Sons",
    year = 2004,
    edition = "3rd",
    pages = "pp. 22"
)

const Jordan1870 = Citation(
    authors = ["C. Jordan"],
    title = "Trait´e des substitutions",
    publisher = "Gauthier-Villars",
    year = 1870,
    edition = "1st"
)

function short(c::Citation)
    l=(first(c.authors))
    for i in (reverse(map(x->Base.string(x),digits(c.year)[[1,2]])))
        l*=i
    end
    return("["*l*"]")
end

function Base.show(io::IO, c::Citation)
    print(io, "Citation: ", c.authors[1], ", \"", c.title, "\", ", c.journal, ", ", c.year)
    if !isempty(c.edition)
        print(io, ", ", c.edition)
    end
    if !isempty(c.publisher)
        print(io, ", ", c.publisher)
    end
    if !isempty(c.url)
        print(io, ", URL: ", c.url)
    end
    if !isempty(c.pages)
        print(io, ", Pages: ", c.pages)
    end
    if !isempty(c.doi)
        print(io, ", DOI: ", c.doi)
    end
end