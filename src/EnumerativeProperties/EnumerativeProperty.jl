
struct EnumerativeProperty{T} 
    name :: String
end


function get_type(EProp::EnumerativeProperty{T}) where {T}
    T
end

function name(EProp::EnumerativeProperty)
    EProp.name
end


function Base.show(io::IO, EProp::EnumerativeProperty)
    print(io,name(EProp))
end

