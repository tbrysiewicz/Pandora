struct Citation
    authors :: Vector{String}
    title :: String
    journal :: String
    year :: Int
end


global NULL_CITATION = Citation([""],"","",0)

