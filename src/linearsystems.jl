struct StateSpace
    inputs
    outputs
    A
    B
    C
    D
end

struct RationalMatrix
    inputs
    outputs
    mat
end

function tfs(SS,s)
    M = SS.C*inv(s*I - SS.A)*SS.B + SS.D
    return RationalMatrix(SS.inputs,SS.outputs,M)
end

function Base.getindex(SS::RationalMatrix,iout,iin)
    in = findfirst(isequal(iin),SS.inputs)
    out = findfirst(isequal(iout),SS.outputs)
    return SS.mat[out,in]
end