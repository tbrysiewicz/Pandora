#Make an EnumerativeProperty for the property you want to compute

function my_actual_computation(I1::Type1,I2::Type2) :: Type3
    #Do computations and return result
end

const my_algorithm = EnumerativeAlgorithm(
    name = "my algorithm",
    input_properties = [EProp1,EProp2],
        #EnumerativeProperties which have Type1 and Type2 as values
    core_function = my_actual_computation,
    output_property = property_computed,
    epistemic_status = :numerical
)