

function justify(EP::EnumerativeProblem,fact::Symbol,justification::String)
    (knowledge(EP).explanations)[fact] = justification
end