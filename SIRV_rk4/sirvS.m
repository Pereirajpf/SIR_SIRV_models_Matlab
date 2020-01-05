%Função Suscetiveis para SIRV
function f=sirvS(ys,yi,beta,ro,yr,mi)
    f=(-beta*ys*yi)-(ro*ys)+(mi*yr);
end