struct Citation
    authors :: Vector{String}
    title :: String
    journal :: String
    year :: Int
end


global NullCitation = Citation([""],"","",0)

