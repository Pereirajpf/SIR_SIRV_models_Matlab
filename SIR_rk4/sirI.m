%funcao para Infetados para SIR
function f=sirI(ys,yi,beta,gama)
    f=(beta*ys*yi)-(gama*yi);
end