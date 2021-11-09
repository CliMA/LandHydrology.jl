module Model2

import .Models: make_rhs, make_tendency_terms, default_initial_conditions, make_update_aux, initialize_states

struct Model2Type <: AbstractModel
    name::model2
end

function Models.make_update_aux(model::Model2Type)
    #makes the update aux function for the auxilary variables Ya.model2.....
end


function Models.make_tendency_terms(model::Model2Type)
    #makes the tendency term functions which update dY.model2 in place
end


# so you can run this model by itself
function Models.make_rhs(model::Model2Type)
    update_aux! = Models.make_update_aux(model)
    tendency_terms = Models.make_tendency_terms!(model)
    function rhs!(dY,Y,Ya,t)
        update_aux!(Ya,t)
        tendency_terms!(dY,Y,Ya,t)
    end
    return rhs!
end

function Models.default_initial_conditions(model::Model2Type)
    # returns a fieldvector for Y, Ya, such that Y.model2 = the state, Ya.model2 = the aux state
end

function Models.initialize_states(model::Model2Type,f::Function, t0::Real)
    # returns a fieldvector for Y, Ya, such that Y.model2 = the state, Ya.model2 = the aux state
end

end
